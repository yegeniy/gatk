package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Code for estimating the reference confidence
 *
 * This code can estimate the probability that the data for a single sample is consistent with a
 * well-determined REF/REF diploid genotype.
 *
 */
public final class ReferenceConfidenceModel {

    private final GenomeLocParser genomeLocParser;

    private final SampleList samples;
    private final int indelInformativeDepthIndelSize;

    private static final boolean WRITE_DEBUGGING_BAM = false;
    private final SAMFileWriter debuggingWriter;

    private static final byte REF_MODEL_DELETION_QUAL = (byte) 30;

    /**
     * Create a new ReferenceConfidenceModel
     *
     * @param genomeLocParser how we create genome locs
     * @param samples the list of all samples we'll be considering with this model
     * @param header the SAMFileHeader describing the read information (used for debugging)
     * @param indelInformativeDepthIndelSize the max size of indels to consider when calculating indel informative depths
     */
    public ReferenceConfidenceModel(final GenomeLocParser genomeLocParser,
                                    final SampleList samples,
                                    final SAMFileHeader header,
                                    final int indelInformativeDepthIndelSize) {
        if ( genomeLocParser == null ) throw new IllegalArgumentException("genomeLocParser cannot be null");
        if ( samples == null ) throw new IllegalArgumentException("samples cannot be null");
        if ( samples.numberOfSamples() == 0) throw new IllegalArgumentException("samples cannot be empty");
        if ( header == null ) throw new IllegalArgumentException("header cannot be empty");
        if ( indelInformativeDepthIndelSize < 0) throw new IllegalArgumentException("indelInformativeDepthIndelSize must be >= 1 but got " + indelInformativeDepthIndelSize);

        this.genomeLocParser = genomeLocParser;
        this.samples = samples;
        this.indelInformativeDepthIndelSize = indelInformativeDepthIndelSize;

        if ( WRITE_DEBUGGING_BAM ) {
            final SAMFileWriterFactory factory = new SAMFileWriterFactory();
            factory.setCreateIndex(true);
            debuggingWriter = factory.makeBAMWriter(header, false, new File("refCalc.bam"));
        } else {
            debuggingWriter = null;
        }
    }

    /**
     * Get the VCF header lines to include when emitting reference confidence values via calculateRefConfidence
     * @return a non-null set of VCFHeaderLines
     */
    public Set<VCFHeaderLine> getVCFHeaderLines() {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>();
        headerLines.add(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME, "Represents any possible alternative allele at this location"));
        return headerLines;
    }

    /**
     * Close down this reference model, closing down any debugging information opened during execution
     */
    public void close() {
        if ( debuggingWriter != null ) debuggingWriter.close();
    }


    /**
     * Calculate the reference confidence for a single sample given the its read data
     *
     * Returns a list of variant contexts, one for each position in the {@code activeRegion.getLoc()}, each containing
     * detailed information about the certainty that the sample is hom-ref for each base in the region.
     *
     *
     *
     * @param refHaplotype the reference haplotype, used to get the reference bases across activeRegion.getLoc()
     * @param calledHaplotypes a list of haplotypes that segregate in this region, for realignment of the reads in the
     *                         readLikelihoods, corresponding to each reads best haplotype.  Must contain the refHaplotype.
     * @param paddedReferenceLoc the location of refHaplotype (which might be larger than activeRegion.getLoc())
     * @param activeRegion the active region we want to get the reference confidence over
     * @param readLikelihoods a map from a single sample to its PerReadAlleleLikelihoodMap for each haplotype in calledHaplotypes
     * @param ploidyModel indicate the ploidy of each sample in {@code stratifiedReadMap}.
     * @param model genotyping model.
     * @param variantCalls calls made in this region.  The return result will contain any variant call in this list in the
     *                     correct order by genomic position, and any variant in this list will stop us emitting a ref confidence
     *                     under any position it covers (for snps and insertions that is 1 bp, but for deletions its the entire ref span)
     * @return an ordered list of variant contexts that spans activeRegion.getLoc() and includes both reference confidence
     *         contexts as well as calls from variantCalls if any were provided
     */
    public List<VariantContext> calculateRefConfidence(final Haplotype refHaplotype,
                                                       final Collection<Haplotype> calledHaplotypes,
                                                       final GenomeLoc paddedReferenceLoc,
                                                       final AssemblyRegion activeRegion,
                                                       final ReadLikelihoods<Haplotype> readLikelihoods,
                                                       final PloidyModel ploidyModel,
                                                       final GenotypingModel model,
                                                       final List<VariantContext> variantCalls) {
        if ( refHaplotype == null ) throw new IllegalArgumentException("refHaplotype cannot be null");
        if ( calledHaplotypes == null ) throw new IllegalArgumentException("calledHaplotypes cannot be null");
        if ( !calledHaplotypes.contains(refHaplotype)) throw new IllegalArgumentException("calledHaplotypes must contain the refHaplotype");
        if ( paddedReferenceLoc == null ) throw new IllegalArgumentException("paddedReferenceLoc cannot be null");
        if ( activeRegion == null ) throw new IllegalArgumentException("activeRegion cannot be null");
        if ( readLikelihoods == null ) throw new IllegalArgumentException("readLikelihoods cannot be null");
        if ( readLikelihoods.numberOfSamples() != 1 ) throw new IllegalArgumentException("readLikelihoods must contain exactly one sample but it contained " + readLikelihoods.sampleCount());
        if ( refHaplotype.length() != activeRegion.getExtendedSpan().size() ) throw new IllegalArgumentException("refHaplotype " + refHaplotype.length() + " and activeRegion location size " + activeRegion.getLocation().size() + " are different");
        if ( ploidyModel == null) throw new IllegalArgumentException("the ploidy model cannot be null");
        if ( model == null) throw new IllegalArgumentException("the genotyping model cannot be null");
        final int ploidy = ploidyModel.samplePloidy(0); // the first sample = the only sample in reference-confidence mode.

        final GenomeLoc refSpan = activeRegion.getSpan();
        final List<ReadPileup> refPileups = getPileupsOverReference(refHaplotype, calledHaplotypes, paddedReferenceLoc, activeRegion, refSpan, readLikelihoods);
        final byte[] ref = refHaplotype.getBases();
        final List<VariantContext> results = new ArrayList<>(refSpan.size());
        final String sampleName = readLikelihoods.getSample(0);

        final int globalRefOffset = refSpan.getStart() - activeRegion.getExtendedSpan().getStart();
        for ( final ReadPileup pileup : refPileups ) {
            final Locatable curPos = pileup.getLocation();
            final int offset = curPos.getStart() - refSpan.getStart();

            final VariantContext overlappingSite = getOverlappingVariantContext(curPos, variantCalls);
            if ( overlappingSite != null && overlappingSite.getStart() == curPos.getStart() ) {
                    results.add(overlappingSite);
            } else {
                // otherwise emit a reference confidence variant context
                final int refOffset = offset + globalRefOffset;
                final byte refBase = ref[refOffset];
                final RefVsAnyResult homRefCalc = calcGenotypeLikelihoodsOfRefVsAny(sampleName,ploidy,model,pileup, refBase, (byte)6, null);
                homRefCalc.capByHomRefLikelihood();

                final Allele refAllele = Allele.create(refBase, true);
                final List<Allele> refSiteAlleles = Arrays.asList(refAllele, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
                final VariantContextBuilder vcb = new VariantContextBuilder("HC", curPos.getContig(), curPos.getStart(), curPos.getStart(), refSiteAlleles);
                final GenotypeBuilder gb = new GenotypeBuilder(sampleName, GATKVariantContextUtils.homozygousAlleleList(refAllele, ploidy));
                gb.AD(homRefCalc.AD_Ref_Any);
                gb.DP(homRefCalc.getDP());

                // genotype likelihood calculation
                final GenotypeLikelihoods snpGLs = GenotypeLikelihoods.fromLog10Likelihoods(homRefCalc.genotypeLikelihoods);
                final int nIndelInformativeReads = calcNIndelInformativeReads(pileup, refOffset, ref, indelInformativeDepthIndelSize);
                final GenotypeLikelihoods indelGLs = getIndelPLs(ploidy,nIndelInformativeReads);

                // now that we have the SNP and indel GLs, we take the one with the least confidence,
                // as this is the most conservative estimate of our certainty that we are hom-ref.
                // For example, if the SNP PLs are 0,10,100 and the indel PLs are 0,100,1000
                // we are very certain that there's no indel here, but the SNP confidence imply that we are
                // far less confident that the ref base is actually the only thing here.  So we take 0,10,100
                // as our GLs for the site.
                final GenotypeLikelihoods leastConfidenceGLs = getGLwithWorstGQ(indelGLs, snpGLs);

                gb.GQ((int) (-10 * leastConfidenceGLs.getLog10GQ(GenotypeType.HOM_REF)));
                gb.PL(leastConfidenceGLs.getAsPLs());
                //gb.attribute(INDEL_INFORMATIVE_DEPTH, nIndelInformativeReads);

                vcb.genotypes(gb.make());
                results.add(vcb.make());
//                logger.info("  => VariantContext " + vcb.make());
            }
        }

        return results;
    }

    /**
     * Get the GenotypeLikelihoods with the least strong corresponding GQ value
     * @param gl1 first to consider (cannot be null)
     * @param gl2 second to consider (cannot be null)
     * @return gl1 or gl2, whichever has the worst GQ
     */
    protected final GenotypeLikelihoods getGLwithWorstGQ(final GenotypeLikelihoods gl1, final GenotypeLikelihoods gl2) {
        return gl1.getLog10GQ(GenotypeType.HOM_REF) > gl2.getLog10GQ(GenotypeType.HOM_REF) ? gl1 : gl2;
    }

    /**
     * Get indel PLs corresponding to seeing N nIndelInformativeReads at this site
     *
     * @param nInformativeReads the number of reads that inform us about being ref without an indel at this site
     * @param ploidy the requested ploidy.
     * @return non-null GenotypeLikelihoods given N
     */
    protected final GenotypeLikelihoods getIndelPLs(final int ploidy, final int nInformativeReads) {
        return indelPLCache(ploidy, nInformativeReads > MAX_N_INDEL_INFORMATIVE_READS ? MAX_N_INDEL_INFORMATIVE_READS : nInformativeReads);
    }

    protected static final int MAX_N_INDEL_INFORMATIVE_READS = 40; // more than this is overkill because GQs are capped at 99 anyway
    private static final int INITIAL_INDEL_LK_CACHE_PLOIDY_CAPACITY = 20;
    private static GenotypeLikelihoods[][] indelPLCache = new GenotypeLikelihoods[INITIAL_INDEL_LK_CACHE_PLOIDY_CAPACITY + 1][];
    private static final double INDEL_ERROR_RATE = -4.5; // 10^-4.5 indel errors per bp

    private final GenotypeLikelihoods indelPLCache(final int ploidy, final int nInformativeReads) {
        return initializeIndelPLCache(ploidy)[nInformativeReads];
    }

    private synchronized GenotypeLikelihoods[] initializeIndelPLCache(final int ploidy) {

        if (indelPLCache.length <= ploidy)
            indelPLCache = Arrays.copyOf(indelPLCache, ploidy << 1);

        if (indelPLCache[ploidy] != null)
            return indelPLCache[ploidy];

        final double denominator =  - MathUtils.Log10Cache.get(ploidy);
        final GenotypeLikelihoods[] result = new GenotypeLikelihoods[MAX_N_INDEL_INFORMATIVE_READS + 1];
        result[0] = GenotypeLikelihoods.fromLog10Likelihoods(new double[ploidy + 1]);
        for( int nInformativeReads = 1; nInformativeReads <= MAX_N_INDEL_INFORMATIVE_READS; nInformativeReads++ ) {
            final byte indelQual = (byte) Math.round((INDEL_ERROR_RATE * -10));
            final double refLikelihood = QualityUtils.qualToProbLog10(indelQual);
            final double altLikelihood = QualityUtils.qualToErrorProbLog10(indelQual);
            double[] PLs = new double[ploidy + 1];
            PLs[0] = nInformativeReads * refLikelihood;
            for (int altCount = 1; altCount <= ploidy; altCount++) {
                final double refLikelihoodAccum = refLikelihood + MathUtils.Log10Cache.get(ploidy - altCount);
                final double altLikelihoodAccum = altLikelihood + MathUtils.Log10Cache.get(altCount);
                PLs[altCount] = nInformativeReads * (MathUtils.approximateLog10SumLog10(refLikelihoodAccum ,altLikelihoodAccum) + denominator);
            }
            result[nInformativeReads] = GenotypeLikelihoods.fromLog10Likelihoods(PLs);
        }
        indelPLCache[ploidy] = result;
        return result;
    }

    /**
     * Calculate the genotype likelihoods for the sample in pileup for being hom-ref contrasted with being ref vs. alt
     *
     * @param sampleName target sample name.
     * @param ploidy target sample ploidy.
     * @param genotypingModel model to calculate likelihoods and genotypes.
     * @param pileup the read backed pileup containing the data we want to evaluate
     * @param refBase the reference base at this pileup position
     * @param minBaseQual the min base quality for a read in the pileup at the pileup position to be included in the calculation
     * @param hqSoftClips running average data structure (can be null) to collect information about the number of high quality soft clips
     * @return a RefVsAnyResult genotype call.
     */
    public RefVsAnyResult calcGenotypeLikelihoodsOfRefVsAny(final String sampleName, final int ploidy,
                                                        final GenotypingModel genotypingModel,
                                                        final ReadPileup pileup, final byte refBase, final byte minBaseQual, final MathUtils.RunningAverage hqSoftClips) {
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Allele.create(refBase, true), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        // Notice that the sample name is rather irrelevant as this information is never used, just need to be the same in both lines bellow.

        final int maximumReadCount = pileup.getReads().size();

        final List<GATKRead> reads = new ArrayList<>(maximumReadCount);
        final double[][] likelihoods = new double[2][maximumReadCount];
        final int[] adCounts = new int[2];
        int nextIndex = 0;
        for (final PileupElement p : pileup) {
            final byte qual = p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual();
            if (!p.isDeletion() && qual <= minBaseQual)
                continue;
            final GATKRead read = p.getRead();
            reads.add(read);
            final boolean isAlt = p.getBase() != refBase || p.isDeletion() || p.isBeforeDeletionStart()
                    || p.isAfterDeletionEnd() || p.isBeforeInsertion() || p.isAfterInsertion() || p.isNextToSoftClip();
            final int bestAllele;
            final int worstAllele;
            if (isAlt) {
                bestAllele = 1;
                worstAllele = 0;
            } else {
                bestAllele = 0;
                worstAllele = 1;
            }

            likelihoods[bestAllele][nextIndex] = QualityUtils.qualToProbLog10(qual);
            likelihoods[worstAllele][nextIndex++] = QualityUtils.qualToErrorProbLog10(qual) + MathUtils.LOG_ONE_THIRD;
            adCounts[bestAllele]++;
            if (isAlt && hqSoftClips != null && p.isNextToSoftClip())
                hqSoftClips.add(AlignmentUtils.calcNumHighQualitySoftClips(read, (byte) 28));
        }

        final Map<String,List<GATKRead>> sampleToReads = Collections.singletonMap(sampleName, reads);
        final ReadLikelihoods<Allele> readLikelihoods = new ReadLikelihoods<>(new IndexedSampleList(sampleName),alleleList,sampleToReads);
        final LikelihoodMatrix<Allele> sampleLikelihoods = readLikelihoods.sampleMatrix(0);
        final int readCount = sampleLikelihoods.numberOfReads();
        for (int i = 0; i < readCount; i++) {
            sampleLikelihoods.set(0,i,likelihoods[0][i]);
            sampleLikelihoods.set(1,i,likelihoods[1][i]);
        }

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(new IndexedSampleList(sampleName),ploidy);
        final GenotypingLikelihoods<Allele> genotypingLikelihoods = genotypingModel.calculateLikelihoods(alleleList, new GenotypingData<>(ploidyModel, readLikelihoods));
        final double[] genotypeLikelihoodArray = genotypingLikelihoods.sampleLikelihoods(0).getAsVector();
        final RefVsAnyResult result = new RefVsAnyResult(genotypeLikelihoodArray.length);
        System.arraycopy(genotypeLikelihoodArray, 0, result.genotypeLikelihoods, 0, genotypeLikelihoodArray.length);
        System.arraycopy(adCounts, 0, result.AD_Ref_Any, 0, 2);
        return result;
    }

    /**
     * Get a list of pileups that span the entire active region span, in order, one for each position
     */
    private List<ReadPileup> getPileupsOverReference(final Haplotype refHaplotype,
                                                           final Collection<Haplotype> calledHaplotypes,
                                                           final GenomeLoc paddedReferenceLoc,
                                                           final AssemblyRegion activeRegion,
                                                           final GenomeLoc activeRegionSpan,
                                                           final ReadLikelihoods<Haplotype> readLikelihoods) {

        if ( refHaplotype == null ) throw new IllegalArgumentException("refHaplotype cannot be null");
        if ( calledHaplotypes == null ) throw new IllegalArgumentException("calledHaplotypes cannot be null");
        if ( !calledHaplotypes.contains(refHaplotype)) throw new IllegalArgumentException("calledHaplotypes must contain the refHaplotype");
        if ( paddedReferenceLoc == null ) throw new IllegalArgumentException("paddedReferenceLoc cannot be null");
        if ( activeRegion == null ) throw new IllegalArgumentException("activeRegion cannot be null");
        if ( readLikelihoods == null ) throw new IllegalArgumentException("readLikelihoods cannot be null");
        if ( readLikelihoods.numberOfSamples() != 1 ) throw new IllegalArgumentException("readLikelihoods must contain exactly one sample but it contained " + readLikelihoods.sampleCount());

        final List<GATKRead> reads = activeRegion.getReads();

        if ( debuggingWriter != null )
            for ( final GATKRead read : reads )
                debuggingWriter.addAlignment(read);

        final LocusIteratorByState libs = new LocusIteratorByState(reads.iterator(), LocusIteratorByState.NO_DOWNSAMPLING,
                true, genomeLocParser, SampleListUtils.asSet(samples), false);

        final List<ReadPileup> pileups = new LinkedList<>();
        final int startPos = activeRegionSpan.getStart();
        AlignmentContext next = libs.advanceToLocus(startPos, true);
        for ( int curPos = startPos; curPos <= activeRegionSpan.getStop(); curPos++ ) {
            if ( next != null && next.getLocation().getStart() == curPos ) {
                pileups.add(next.getBasePileup());
                next = libs.hasNext() ? libs.next() : null;
            } else {
                // no data, so we create empty pileups
                pileups.add(new ReadPileup(genomeLocParser.createGenomeLoc(activeRegionSpan.getContig(), curPos)));
            }
        }

        return pileups;
    }

    /**
     * Return the rightmost variant context in maybeOverlapping that overlaps curPos
     *
     * @param curPos non-null genome loc
     * @param maybeOverlapping a collection of variant contexts that might overlap curPos
     * @return a VariantContext, or null if none overlaps
     */
    protected final VariantContext getOverlappingVariantContext(final GenomeLoc curPos, final Collection<VariantContext> maybeOverlapping) {
        VariantContext overlaps = null;
        for ( final VariantContext vc : maybeOverlapping ) {
            if ( genomeLocParser.createGenomeLoc(vc).overlapsP(curPos) ) {
                if ( overlaps == null || vc.getStart() > overlaps.getStart() ) {
                    overlaps = vc;
                }
            }
        }
        return overlaps;
    }

    /**
     * Compute the sum of mismatching base qualities for readBases aligned to refBases at readStart / refStart
     * assuming no insertions or deletions in the read w.r.t. the reference
     *
     * @param readBases non-null bases of the read
     * @param readQuals non-null quals of the read
     * @param readStart the starting position of the read (i.e., that aligns it to a position in the reference)
     * @param refBases the reference bases
     * @param refStart the offset into refBases that aligns to the readStart position in readBases
     * @param maxSum if the sum goes over this value, return immediately
     * @return the sum of quality scores for readBases that mismatch their corresponding ref bases
     */
    protected final int sumMismatchingQualities(final byte[] readBases,
                                                final byte[] readQuals,
                                                final int readStart,
                                                final byte[] refBases,
                                                final int refStart,
                                                final int maxSum) {
        final int n = Math.min(readBases.length - readStart, refBases.length - refStart);
        int sum = 0;

        for ( int i = 0; i < n; i++ ) {
            final byte readBase = readBases[readStart + i];
            final byte refBase  = refBases[refStart + i];
            if ( readBase != refBase ) {
                sum += readQuals[readStart + i];
                if ( sum > maxSum ) // abort early
                    return sum;
            }
        }

        return sum;
    }

    /**
     * Compute whether a read is informative to eliminate an indel of size <= maxIndelSize segregating at readStart/refStart
     *
     * @param readBases non-null bases of the read
     * @param readQuals non-null quals of the read
     * @param readStart the starting position of the read (i.e., that aligns it to a position in the reference)
     * @param refBases the reference bases
     * @param refStart the offset into refBases that aligns to the readStart position in readBases
     * @param maxIndelSize the max indel size to consider for the read to be informative
     * @return true if read can eliminate the possibility that there's an indel of size <= maxIndelSize segregating at refStart
     */
    protected boolean isReadInformativeAboutIndelsOfSize(final byte[] readBases,
                                                         final byte[] readQuals,
                                                         final int readStart,
                                                         final byte[] refBases,
                                                         final int refStart,
                                                         final int maxIndelSize) {
        // fast exit when n bases left < maxIndelSize
        if( readBases.length - readStart < maxIndelSize || refBases.length - refStart < maxIndelSize ) {
            return false;
        }

        final int baselineMMSum = sumMismatchingQualities(readBases, readQuals, readStart, refBases, refStart, Integer.MAX_VALUE);

        // consider each indel size up to max in term, checking if an indel that deletes either the ref bases (deletion
        // or read bases (insertion) would fit as well as the origin baseline sum of mismatching quality scores
        for ( int indelSize = 1; indelSize <= maxIndelSize; indelSize++ ) {
            for ( final boolean checkInsertion : Arrays.asList(true, false) ) {
                final int readI, refI;
                if ( checkInsertion ) {
                    readI = readStart + indelSize;
                    refI = refStart;
                } else {
                    readI = readStart;
                    refI = refStart + indelSize;
                }

                final int score = sumMismatchingQualities(readBases, readQuals, readI, refBases, refI, baselineMMSum);
                if ( score <= baselineMMSum )
                    return false;
            }
        }

        return true;
    }

    /**
     * Calculate the number of indel informative reads at pileup
     *
     * @param pileup a pileup
     * @param pileupOffsetIntoRef the position of the pileup in the reference
     * @param ref the ref bases
     * @param maxIndelSize maximum indel size to consider in the informativeness calculation
     * @return an integer >= 0
     */
    protected final int calcNIndelInformativeReads(final ReadPileup pileup, final int pileupOffsetIntoRef, final byte[] ref, final int maxIndelSize) {
        int nInformative = 0;
        for ( final PileupElement p : pileup ) {
            final GATKRead read = p.getRead();
            final int offset = p.getOffset();

            // doesn't count as evidence
            if ( p.isBeforeDeletionStart() || p.isBeforeInsertion() || p.isDeletion() )
                continue;

            // todo -- this code really should handle CIGARs directly instead of relying on the above tests
            if ( isReadInformativeAboutIndelsOfSize(read.getBases(), read.getBaseQualities(), offset, ref, pileupOffsetIntoRef, maxIndelSize) ) {
                nInformative++;
                if( nInformative > MAX_N_INDEL_INFORMATIVE_READS ) {
                    return MAX_N_INDEL_INFORMATIVE_READS;
                }
            }
        }
        return nInformative;
    }

    /**
     * Create a reference haplotype for an active region
     *
     * @param activeRegion the active region
     * @param refBases the ref bases
     * @param paddedReferenceLoc the location spanning of the refBases -- can be longer than activeRegion.getLocation()
     * @return a reference haplotype
     */
    public static Haplotype createReferenceHaplotype(final AssemblyRegion activeRegion, final byte[] refBases, final GenomeLoc paddedReferenceLoc) {
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final int alignmentStart = activeRegion.getExtendedSpan().getStart() - paddedReferenceLoc.getStart();
        if ( alignmentStart < 0 ) throw new IllegalStateException("Bad alignment start in createReferenceHaplotype " + alignmentStart);
        refHaplotype.setAlignmentStartHapwrtRef(alignmentStart);
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);
        return refHaplotype;
    }
}
