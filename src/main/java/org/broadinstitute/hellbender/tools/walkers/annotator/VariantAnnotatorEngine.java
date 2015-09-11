package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException.BadArgumentValue;
import org.broadinstitute.hellbender.tools.walkers.annotator.interfaces.AnnotationGroup;
import org.broadinstitute.hellbender.tools.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.interfaces.VariantAnnotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.walkers.annotator.interfaces.AnnotationGroup.ActiveRegionBased;
import static org.broadinstitute.hellbender.tools.walkers.annotator.interfaces.AnnotationGroup.Standard;

/**
 * The class responsible for computing annotations for variants.
 */
public final class VariantAnnotatorEngine {
    private final List<InfoFieldAnnotation> infoAnnotations;
    private final List<GenotypeAnnotation> genotypeAnnotations;

    private final VariantOverlapAnnotator variantOverlapAnnotator;

    private VariantAnnotatorEngine(final AnnotationManager annots, final FeatureInput<VariantContext> dbSNPBinding){
        infoAnnotations = annots.createInfoFieldAnnotations();
        genotypeAnnotations = annots.createGenotypeAnnotations();
        variantOverlapAnnotator = initializeDBs(dbSNPBinding);
    }

    /**
     * Makes the engine for all known annotation types (minus the excluded ones).
     */
    public static VariantAnnotatorEngine ofAll(final List<String> annotationsToExclude, final FeatureInput<VariantContext> dbSNPBinding) {
        Utils.nonNull(annotationsToExclude, "annotationsToExclude is null");
        return new VariantAnnotatorEngine(AnnotationManager.ofAll(annotationsToExclude), dbSNPBinding);
    }

    /**
     * Makes the engine for given annotation types and groups (minus the excluded ones).
     */
    public static VariantAnnotatorEngine ofSelected(final List<String> annotationGroupsToUse, final List<String> annotationsToUse, final List<String> annotationsToExclude, final FeatureInput<VariantContext> dbSNPBinding) {
        Utils.nonNull(annotationGroupsToUse, "annotationGroupsToUse is null");
        Utils.nonNull(annotationsToUse, "annotationsToUse is null");
        Utils.nonNull(annotationsToExclude, "annotationsToExclude is null");
        return new VariantAnnotatorEngine(AnnotationManager.ofSelected(annotationGroupsToUse, annotationsToUse, annotationsToExclude), dbSNPBinding);
    }

    private VariantOverlapAnnotator initializeDBs(final FeatureInput<VariantContext> dbSNPBinding) {
        final Map<FeatureInput<VariantContext>, String> overlapBindings = new LinkedHashMap<>();
        if ( dbSNPBinding != null ) {
            overlapBindings.put(dbSNPBinding, VCFConstants.DBSNP_KEY); // add overlap detection with DBSNP by default
        }

        return new VariantOverlapAnnotator(dbSNPBinding, overlapBindings);
    }

    /**
     * Returns the list of genotype annotations that will be applied.
     * Note: The returned list is unmodifiable.
     */
    public List<GenotypeAnnotation> getGenotypeAnnotations() {
        return Collections.unmodifiableList(genotypeAnnotations);
    }

    /**
     * Returns the list of info annotations that will be applied.
     * Note: The returned list is unmodifiable.
     */
    public List<InfoFieldAnnotation> getInfoAnnotations() {
        return Collections.unmodifiableList(infoAnnotations);
    }

    /**
     * Returns the set of desciptions to be added to the VCFHeader line (for all annotations in this engine).
     */
    public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {
        final Set<VCFHeaderLine> descriptions = new HashSet<>();

        for ( final InfoFieldAnnotation annotation : infoAnnotations) {
            descriptions.addAll(annotation.getDescriptions());
        }
        for ( final GenotypeAnnotation annotation : genotypeAnnotations) {
            descriptions.addAll(annotation.getDescriptions());
        }
        for ( final String db : variantOverlapAnnotator.getOverlapNames() ) {
            if ( VCFStandardHeaderLines.getInfoLine(db, false) != null ) {
                descriptions.add(VCFStandardHeaderLines.getInfoLine(db));
            } else {
                descriptions.add(new VCFInfoHeaderLine(db, 0, VCFHeaderLineType.Flag, db + " Membership"));
            }
        }

        return descriptions;
    }

    /**
     * Annotates the given variant context - adds all annotations that satisfy the predicate.
     */
    public VariantContext annotateContext(final FeatureContext features,
                                          final ReferenceContext ref,
                                          final VariantContext vc,
                                          final Map<String,PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap,
                                          final Predicate<VariantAnnotation> addAnnot) {
        Utils.nonNull(features, "features cannot be null");
        Utils.nonNull(vc, "vc cannot be null");
        final Map<String, Object> annotValues = new LinkedHashMap<>(vc.getAttributes());

        for ( final InfoFieldAnnotation annotationType : infoAnnotations) {
            if (addAnnot.test(annotationType)){
                final Map<String, Object> annotationsFromCurrentType = annotationType.annotate(ref, vc, perReadAlleleLikelihoodMap);
                if ( annotationsFromCurrentType != null ) {
                    annotValues.putAll(annotationsFromCurrentType);
                }
            }
        }
        final VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(annotValues);

        // annotate genotypes, creating another new VC in the process
        final VariantContext annotated = builder.genotypes(annotateGenotypes(ref, vc, perReadAlleleLikelihoodMap, addAnnot)).make();

        // annotate db occurrences
        return variantOverlapAnnotator.annotateOverlaps(features, variantOverlapAnnotator.annotateRsID(features, annotated));
    }

    private GenotypesContext annotateGenotypes(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final Map<String,PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                               final Predicate<VariantAnnotation> addAnnot) {
        if ( genotypeAnnotations.isEmpty() ) {
            return vc.getGenotypes();
        }

        final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = null;
            if (stratifiedPerReadAlleleLikelihoodMap != null) {
                perReadAlleleLikelihoodMap = stratifiedPerReadAlleleLikelihoodMap.get(genotype.getSampleName());
            }

            final GenotypeBuilder gb = new GenotypeBuilder(genotype);
            for ( final GenotypeAnnotation annotation : genotypeAnnotations) {
                if (addAnnot.test(annotation)) {
                    annotation.annotate(ref, vc, genotype, gb, perReadAlleleLikelihoodMap);
                }
            }
            genotypes.add(gb.make());
        }

        return genotypes;
    }

    private static final class AnnotationManager {

        private final List<String> annotationGroupsToUse;
        private final List<String> annotationsToUse;

        private AnnotationManager(final List<String> annotationGroupsToUse, final List<String> annotationsToUse){
            this.annotationGroupsToUse = annotationGroupsToUse;
            this.annotationsToUse = annotationsToUse;

            final Set<String> unknownAnnots = Sets.difference(new HashSet<>(annotationsToUse), new HashSet<>(AnnotationManager.getAllAnnotationNames()));
            if (!unknownAnnots.isEmpty()){
                throw new BadArgumentValue("annotation", "Annotations " + unknownAnnots + " were not found; please check that you have specified the annotation name correctly");
            }

            final Set<String> unknownGroups =  Sets.difference(new HashSet<>(annotationGroupsToUse), new HashSet<>(AnnotationManager.getAllAnnotationGroupNames()));
            if (!unknownGroups.isEmpty()){
                throw new BadArgumentValue("group", "Unknown annotation group in " + unknownGroups + ". Known groups are " + Arrays.toString(AnnotationGroup.values()) );
            }

        }

        static AnnotationManager ofSelected(final List<String> annotationGroupsToUse, final List<String> annotationsToUse, final List<String> annotationsToExclude){
            final List<String> groups = new ArrayList<>(annotationGroupsToUse);//make copy
            final List<String> annots = new ArrayList<>(annotationsToUse);//make copy
            annots.removeAll(annotationsToExclude);
            return new AnnotationManager(groups, annots);
        }

        static AnnotationManager ofAll(final List<String> annotationsToExclude){
            final List<String> groups = getAllAnnotationGroupNames();
            final List<String> annots = getAllAnnotationNames();
            annots.removeAll(annotationsToExclude);
            return new AnnotationManager(groups, annots);
        }

        private static List<String> getAllAnnotationNames() {
            final Set<VariantAnnotation> union = Sets.union(new HashSet<>(makeAllGenotypeAnnotations()), new HashSet<>(AnnotationManager.makeAllInfoFieldAnnotations()));
            return union.stream().map(a -> a.getClass().getSimpleName()).collect(Collectors.toList());
        }

        public static List<String> getAllAnnotationGroupNames() {
            return Stream.of(AnnotationGroup.values()).map(v -> v.name()).collect(Collectors.toList());
        }

        public List<InfoFieldAnnotation> createInfoFieldAnnotations() {
            final List<InfoFieldAnnotation> all = makeAllInfoFieldAnnotations();
            return filterAnnotations(all);
        }

        private static List<InfoFieldAnnotation> makeAllInfoFieldAnnotations() {
            return Arrays.asList(
                    new InbreedingCoeff(Standard, ActiveRegionBased),
                    new Coverage(Standard, ActiveRegionBased),
                    new ReadPosRankSumTest(Standard, ActiveRegionBased),
                    new BaseQualityRankSumTest(Standard, ActiveRegionBased),
                    new MappingQualityRankSumTest(Standard, ActiveRegionBased),
                    new ChromosomeCounts(Standard, ActiveRegionBased),
                    new StrandOddsRatio(Standard, ActiveRegionBased),
                    new FisherStrand(Standard, ActiveRegionBased),
                    new RMSMappingQuality(Standard, ActiveRegionBased),
                    new QualByDepth(Standard, ActiveRegionBased),
                    new MappingQualityZero(ActiveRegionBased),
                    new SampleList(),
                    new ClippingRankSumTest(ActiveRegionBased),
                    new LikelihoodRankSumTest(ActiveRegionBased),
                    new TandemRepeat(ActiveRegionBased),
                    new GenotypeSummaries(ActiveRegionBased)
                    //        new PossibleDeNovo(),  //TODO enable this
            );
        }

        public List<GenotypeAnnotation> createGenotypeAnnotations() {
            final List<GenotypeAnnotation> all = makeAllGenotypeAnnotations();
            return filterAnnotations(all);
        }

        private static List<GenotypeAnnotation> makeAllGenotypeAnnotations() {
            return Arrays.asList(
                    new BaseQualitySumPerAlleleBySample(Standard, ActiveRegionBased),
                    new DepthPerAlleleBySample(Standard, ActiveRegionBased),
                    new OxoGReadCounts(ActiveRegionBased)
            );
        }

        private <T extends VariantAnnotation> List<T> filterAnnotations(final List<T> all) {
            final SortedSet<T> annotations = new TreeSet<>(Comparator.comparing(t -> t.getClass().getSimpleName()));

            for (final T t : all){
                //if any group matches requested groups, it's in
                if (t.getGroups().stream().anyMatch(g -> annotationGroupsToUse.contains(g.name()))){
                    annotations.add(t);
                }
                if (annotationsToUse.contains(t.getClass().getSimpleName())){
                    annotations.add(t);
                }
            }

            return Collections.unmodifiableList(new ArrayList<>(annotations));
        }
    }

}
