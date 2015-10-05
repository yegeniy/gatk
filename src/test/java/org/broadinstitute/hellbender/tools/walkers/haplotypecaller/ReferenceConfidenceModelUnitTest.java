package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.api.client.repackaged.com.google.common.base.Strings;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.UnvalidatingGenomeLoc;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class ReferenceConfidenceModelUnitTest extends BaseTest {
    GenomeLocParser parser;
    final String RGID = "ID1";
    SAMReadGroupRecord rg;
    final String sample = "NA12878";
    final SampleList samples = SampleList.singletonSampleList(sample);
    SAMFileHeader header;
    ReferenceConfidenceModel model;

    @BeforeClass
    public void setUp() throws Exception {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        rg = new SAMReadGroupRecord(RGID);
        rg.setSample(sample);
        header.addReadGroup(rg);
        parser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @BeforeMethod
    public void setupModel() {
        model = new ReferenceConfidenceModel(parser, samples, header, 10);
    }

    @DataProvider(name = "CalcNIndelInformativeReadsData")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        { // very basic testing
            final String ref  = "ACGT";
            final String read = "ACGT";
            tests.add(new Object[]{read, ref, 1, Arrays.asList(1, 1, 1, 0)});
            tests.add(new Object[]{read, ref, 2, Arrays.asList(1, 1, 0, 0)});
            tests.add(new Object[]{read, ref, 3, Arrays.asList(1, 0, 0, 0)});
            tests.add(new Object[]{read, ref, 4, Arrays.asList(0, 0, 0, 0)});
        }

        { // actually interesting case where some sites aren't informative
            final String ref   = "NNAAAANN";
            final String read1 = "NNA";
            final String read2 = "NNAA";
            final String read3 = "NNAAA";
            final String read4 = "NNAAAA";
            final String read5 = "NNAAAAN";
            tests.add(new Object[]{read1, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read2, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read3, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read4, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read5, ref, 1, Arrays.asList(1, 1, 1, 1, 1, 1, 0, 0)});
        }

        {
            for ( final String repeatUnit : Arrays.asList("A", "CA", "TAG", "TAGC", "TCAGA")) {
                final String anchor = Strings.repeat("N", repeatUnit.length());
                for ( int nUnits = 1; nUnits < 10; nUnits++ ) {
                    final String repeat = Strings.repeat(repeatUnit, nUnits);
                    final String ref = anchor + repeat + anchor;
                    for ( int readLen = repeatUnit.length(); readLen < repeat.length(); readLen++ ) {
                        final String read = anchor + repeat.substring(0, readLen);
                        final List<Integer> expected = new LinkedList<>();
                        for ( int i = 0; i < anchor.length(); i++ ) expected.add(1);
                        for ( int i = 0; i < repeat.length(); i++ ) expected.add(readLen == repeat.length() ? 1 : 0);
                        for ( int i = 0; i < anchor.length(); i++ ) expected.add(0);
                        tests.add(new Object[]{read, ref, repeatUnit.length(), expected});

                        final List<Integer> result = new ArrayList<>(Collections.nCopies(ref.length() - anchor.length(), 1));
                        result.addAll(Collections.nCopies(anchor.length(), 0));
                        tests.add(new Object[]{ref, ref, repeatUnit.length(), result});
                    }
                }

            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CalcNIndelInformativeReadsData")
    public void testCalcNIndelInformativeReads(final String readBases, final String ref, final int maxIndelSize, final List<Integer> expected ) {
        final byte qual = (byte)30;
        final byte[] quals = Utils.dupBytes(qual, readBases.length());

        for ( int i = 0; i < readBases.getBytes().length; i++ ) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(readBases.getBytes(), quals, readBases.length() + "M");
            final GenomeLoc loc = new UnvalidatingGenomeLoc("20", 0, i, i);
            final ReadPileup pileup = new ReadPileup(loc, Collections.singletonList(read), i);
            final int actual = model.calcNIndelInformativeReads(pileup, i, ref.getBytes(), maxIndelSize);
            Assert.assertEquals(actual, (int) expected.get(i), "failed at position " + i);
        }
    }

    @Test
    public void testClose() {
        model.close();
    }

    @Test
    public void testWorstGL() {
        final GenotypeLikelihoods gq10 = GenotypeLikelihoods.fromPLField("0,10,100");
        final GenotypeLikelihoods gq20 = GenotypeLikelihoods.fromPLField("0,20,200");
        final GenotypeLikelihoods gq0 = GenotypeLikelihoods.fromPLField("20,0,200");

        Assert.assertSame(model.getGLwithWorstGQ(gq10, gq20), gq10);
        Assert.assertSame(model.getGLwithWorstGQ(gq20, gq10), gq10);
        Assert.assertSame(model.getGLwithWorstGQ(gq10, gq0), gq0);
        Assert.assertSame(model.getGLwithWorstGQ(gq0, gq10), gq0);
    }

    @Test
    public void testIndelLikelihoods() {
        GenotypeLikelihoods prev = model.getIndelPLs(HomoSapiensConstants.DEFAULT_PLOIDY,0);
        Assert.assertEquals(prev.getAsPLs(), new int[]{0, 0, 0});
        Assert.assertEquals(-10 * prev.getLog10GQ(GenotypeType.HOM_REF), 0.0);

        for ( int i = 1; i <= ReferenceConfidenceModel.MAX_N_INDEL_INFORMATIVE_READS; i++ ) {
            final GenotypeLikelihoods current = model.getIndelPLs(HomoSapiensConstants.DEFAULT_PLOIDY,i);
            final double prevGQ = -10 * prev.getLog10GQ(GenotypeType.HOM_REF);
            final double currGQ = -10 * current.getLog10GQ(GenotypeType.HOM_REF);
            Assert.assertTrue(prevGQ < currGQ, "GQ Failed with prev " + prev + " curr " + current + " at " + i);
            Assert.assertTrue(prev.getAsPLs()[1] < current.getAsPLs()[1], "het PL failed with prev " + prev + " curr " + current + " at " + i);
            Assert.assertTrue(prev.getAsPLs()[2] < current.getAsPLs()[2], "hom-var PL Failed with prev " + prev + " curr " + current + " at " + i);
//            logger.warn("result at " + i + " is " + current);
            prev = current;
        }
    }

    @Test
    public void testOverlappingVariantContext() {
        final VariantContext vc10 = GATKVariantContextUtils.makeFromAlleles("test", "chr1", 10, Arrays.asList("A", "C"));
        final VariantContext vc13 = GATKVariantContextUtils.makeFromAlleles("test", "chr1", 13, Arrays.asList("A", "C"));
        final VariantContext vc12_15 = GATKVariantContextUtils.makeFromAlleles("test", "chr1", 12, Arrays.asList("ACAT", "A"));
        final VariantContext vc18 = GATKVariantContextUtils.makeFromAlleles("test", "chr1", 18, Arrays.asList("A", "ACAT"));

        final List<VariantContext> calls = Arrays.asList(vc13, vc12_15, vc18, vc10);

        checkOverlapping(8, calls, null);
        checkOverlapping(9, calls, null);
        checkOverlapping(10, calls, vc10);
        checkOverlapping(11, calls, null);
        checkOverlapping(12, calls, vc12_15);
        checkOverlapping(13, calls, vc13);
        checkOverlapping(14, calls, vc12_15);
        checkOverlapping(15, calls, vc12_15);
        checkOverlapping(16, calls, null);
        checkOverlapping(17, calls, null);
        checkOverlapping(18, calls, vc18);
        checkOverlapping(19, calls, null);
        checkOverlapping(20, calls, null);
    }

    private void checkOverlapping(final int pos, Collection<VariantContext> calls, final VariantContext expected) {
        final GenomeLoc loc = parser.createGenomeLoc(parser.getContigs().getSequences().get(0).getSequenceName(), pos, pos);
        final VariantContext actual = model.getOverlappingVariantContext(loc, calls);
        Assert.assertEquals(actual, expected);
    }

    //
    // test reference calculation
    //
    private class RefConfData {
        final String ref;
        final int extension;
        final Haplotype refHap;
        final GenomeLoc refLoc, paddedRefLoc;
        final AssemblyRegion region;
        int readCounter = 0;

        private RefConfData(String ref, int extension) {
            this.ref = ref;
            this.extension = extension;

            refLoc = parser.createGenomeLoc("chr1", getStart(), getEnd());
            paddedRefLoc = parser.createGenomeLoc("chr1", getStart() - extension, getEnd() + extension);
            region = new AssemblyRegion(getRefLoc(), parser, extension);
            final String pad = Strings.repeat("N", extension);
            refHap = ReferenceConfidenceModel.createReferenceHaplotype(getActiveRegion(), (pad + ref + pad).getBytes(), getPaddedRefLoc());
        }

        public GenomeLoc getRefLoc() { return refLoc; }
        public GenomeLoc getPaddedRefLoc() { return paddedRefLoc; }
        public AssemblyRegion getActiveRegion() { return region; }
        public Haplotype getRefHap() { return refHap; }
        public int getStart() { return 100; }
        public int getEnd() { return getStart() + getRefLength() - 1; }
        public byte[] getRefBases() { return ref.getBytes(); }
        public int getRefLength() { return ref.length(); }

        public GATKRead makeRead(final int start, final int length) {
            final byte[] quals = Utils.dupBytes((byte)30, length);
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read " + readCounter++, 0, start + getStart(), ref.substring(start, start + length).getBytes(), quals, length + "M");
            read.setReadGroup(rg);
            return read;
        }
    }


    @DataProvider(name = "RefConfidenceData")
    public Object[][] makeRefConfidenceData() {
        List<Object[]> tests = new ArrayList<>();

        for ( int i = 0; i < 10; i++ ) {
            for ( final int extension : Arrays.asList(0, 10) ) {
                tests.add(new Object[]{i, extension});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RefConfidenceData")
    public void testRefConfidenceBasic(final int nReads, final int extension) {
        final RefConfData data = new RefConfData("ACGTAACCGGTT", extension);
        final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
        final List<VariantContext> calls = Collections.emptyList();

        for ( int i = 0; i < nReads; i++ ) {
            data.getActiveRegion().add(data.makeRead(0, data.getRefLength()));
        }

        final ReadLikelihoods<Haplotype> likelihoods = HaplotypeCaller.createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final GenotypingModel genotypingModel = new InfiniteRandomMatingPopulationModel();
        final List<Integer> expectedDPs = Collections.nCopies(data.getActiveRegion().getSpan().size(), nReads);
        final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, genotypingModel, calls);
        checkReferenceModelResult(data, contexts, expectedDPs, calls);
    }

    @Test
    public void testRefConfidencePartialReads() {

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final GenotypingModel genotypingModel = new InfiniteRandomMatingPopulationModel();
        final String ref = "ACGTAACCGGTT";
        for ( int readLen = 3; readLen < ref.length(); readLen++ ) {
            for ( int start = 0; start < ref.length() - readLen; start++ ) {
                final RefConfData data = new RefConfData(ref, 0);
                final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
                final List<VariantContext> calls = Collections.emptyList();

                data.getActiveRegion().add(data.makeRead(start, readLen));
                final ReadLikelihoods<Haplotype> likelihoods = HaplotypeCaller.createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

                final List<Integer> expectedDPs = new ArrayList<>(Collections.nCopies(data.getActiveRegion().getSpan().size(), 0));
                for ( int i = start; i < readLen + start; i++ ) expectedDPs.set(i, 1);
                final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, genotypingModel, calls);
                checkReferenceModelResult(data, contexts, expectedDPs, calls);
            }
        }
    }

    @Test
    public void testRefConfidenceWithCalls() {
        final RefConfData xxxdata = new RefConfData("ACGTAACCGGTT", 0);
        final int start = xxxdata.getStart();
        final int stop = xxxdata.getEnd();

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final GenotypingModel genotypingModel = new InfiniteRandomMatingPopulationModel();

        for ( int nReads = 0; nReads < 2; nReads++ ) {

            final VariantContext vcStart = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start, Arrays.asList("A", "C"));
            final VariantContext vcEnd = GATKVariantContextUtils.makeFromAlleles("test", "chr1", stop, Arrays.asList("A", "C"));
            final VariantContext vcMiddle = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 2, Arrays.asList("A", "C"));
            final VariantContext vcDel = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 4, Arrays.asList("AAC", "A"));
            final VariantContext vcIns = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 8, Arrays.asList("G", "GCG"));

            final List<VariantContext> allCalls = Arrays.asList(vcStart, vcEnd, vcMiddle, vcDel, vcIns);

            for ( int n = 1; n <= allCalls.size(); n++ ) {
                for ( final List<VariantContext> calls : Utils.makePermutations(allCalls, n, false) ) {
//                    logger.warn("Executing " + n + " " + calls.size());
                    final RefConfData data = new RefConfData("ACGTAACCGGTT", 0);
                    final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
                    for ( int i = 0; i < nReads; i++ ) {
                        data.getActiveRegion().add(data.makeRead(0, data.getRefLength()));
                    }

                    final ReadLikelihoods<Haplotype> likelihoods = HaplotypeCaller.createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

                    final List<Integer> expectedDPs = Collections.nCopies(data.getActiveRegion().getSpan().size(), nReads);
                    final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, genotypingModel, calls);
                    checkReferenceModelResult(data, contexts, expectedDPs, calls);
                }
            }
        }
    }

    private void checkReferenceModelResult(final RefConfData data, final List<VariantContext> contexts, final List<Integer> expectedDPs, final List<VariantContext> calls) {
        Assert.assertNotNull(contexts);

        final GenomeLoc loc = data.getActiveRegion().getExtendedSpan();
        final List<Boolean> seenBP = new ArrayList<>(Collections.nCopies(data.getActiveRegion().getSpan().size(), false));

        for ( int i = 0; i < loc.size(); i++ ) {
            final GenomeLoc curPos = parser.createGenomeLoc(loc.getContig(), loc.getStart() + i);
            final VariantContext call = model.getOverlappingVariantContext(curPos, calls);
            final VariantContext refModel = model.getOverlappingVariantContext(curPos, contexts);

            if ( ! data.getActiveRegion().getSpan().containsP(curPos) ) {
                // part of the extended interval, but not the full interval
                Assert.assertNull(refModel);
                continue;
            }

            if ( call != null ) {
                if (call.isVariant() && refModel.getType() ==  VariantContext.Type.SYMBOLIC ) {
                    //Assert.assertEquals(refModel, call, "Should have found call " + call + " but found " + refModel + " instead");
                    Assert.assertTrue(call.getReference().length() > 1); // must be a deletion.
                    Assert.assertTrue(call.getStart() < refModel.getStart()); // the deletion must not start at the same position
                    Assert.assertEquals(call.getReference().getBaseString().substring(refModel.getStart() - call.getStart(),
                            refModel.getStart() - call.getStart() + 1), refModel.getReference().getBaseString(), "" + data.getRefHap()); // the reference must be the same.
                    Assert.assertTrue(refModel.getGenotype(0).getGQ() <= 0); // No confidence in the reference hom-ref call across the deletion
                    Assert.assertEquals(refModel.getAlleles().size(), 2); // the reference and the lonelly <NON_REF>
                    Assert.assertEquals(refModel.getAlleles().get(1), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
                } else {
                    Assert.assertEquals(refModel, call, "Should have found call " + call + " but found " + refModel + " instead");
                }

            } else {
                final int expectedDP = expectedDPs.get(curPos.getStart() - data.getActiveRegion().getSpan().getStart());
                Assert.assertEquals(refModel.getStart(), loc.getStart() + i);
                Assert.assertEquals(refModel.getEnd(), loc.getStart() + i);
                Assert.assertFalse(refModel.hasLog10PError());
                Assert.assertEquals(refModel.getAlternateAlleles().size(), 1);
                Assert.assertEquals(refModel.getAlternateAllele(0), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
                Assert.assertTrue(refModel.hasGenotype(sample));

                final Genotype g = refModel.getGenotype(sample);
                Assert.assertTrue(g.hasAD());
                Assert.assertTrue(g.hasDP());
                Assert.assertEquals(g.getDP(), expectedDP);
                Assert.assertTrue(g.hasGQ());
                Assert.assertTrue(g.hasPL());
            }

            final VariantContext vc = call == null ? refModel : call;
            if ( curPos.getStart() == vc.getStart() ) {
                for ( int pos = vc.getStart(); pos <= vc.getEnd(); pos++ ) {
                    final int j = pos - data.getActiveRegion().getSpan().getStart();
                    Assert.assertFalse(seenBP.get(j));
                    seenBP.set(j, true);
                }
            }
        }

        for ( int i = 0; i < seenBP.size(); i++ ) {
            Assert.assertEquals((boolean) seenBP.get(i), true);
        }
    }
}