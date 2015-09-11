package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.interfaces.AnnotationGroup;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class VariantAnnotatorEngineUnitTest extends BaseTest {
    @Test
    public void testEmpty(){
        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse= Collections.emptyList();
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelected(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

    @Test
    public void testExclude(){
        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final List<String> annotationsToExclude= annotationsToUse;
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelected(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

    @Test
    public void testAll(){
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofAll(annotationsToExclude, dbSNPBinding);
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions();
        Assert.assertFalse(vcfAnnotationDescriptions.isEmpty());
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadAnnot(){
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<String> annotationGroupsToUse = Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList("fred");
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelected(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding);
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadGroup(){
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<String> annotationGroupsToUse = Arrays.asList("fred");
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelected(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding);
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
    }


    private VariantContext makeVC(final Allele refAllele, final Allele altAllele) {
        final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
        final Genotype g = new GenotypeBuilder("sample1", alleles).make();

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(g).make();
    }

    private Map<String, PerReadAlleleLikelihoodMap> makeReadMap(final int ref, final int alt, final Allele refAllele, final Allele altAllele) {
        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        for (int i = 0; i < alt; i++) {
            final GATKRead r = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
            map.add(r, altAllele, -1.0);
            map.add(r, refAllele, -100.0);
        }
        for (int i = 0; i < ref; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
            map.add(read, altAllele, -100.0);
            map.add(read, refAllele, -1.0);
        }

        return Collections.singletonMap("sample1", map);
    }

    @Test
    public void testCoverageAnnotation() throws Exception {
        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelected(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding);

        final int alt = 5;
        final int ref = 3;
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = makeReadMap(ref, alt, refAllele, altAllele);
        final VariantContext resultVC = vae.annotateContext(new FeatureContext(), null, makeVC(refAllele, altAllele), perReadAlleleLikelihoodMap, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));
    }

    @Test
    public void testAnnotationsAsActiveRegion() throws Exception {
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>("dbsnp", Collections.emptyMap(), file);

        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName(),
                                                            DepthPerAlleleBySample.class.getSimpleName(),
                                                            SampleList.class.getSimpleName());
        final List<String> annotationsToExclude= Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelected(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding);

        final int alt = 5;
        final int ref = 3;
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = makeReadMap(ref, alt, refAllele, altAllele);
        final VariantContext resultVC = vae.annotateContext(new FeatureContext(), null, makeVC(refAllele, altAllele), perReadAlleleLikelihoodMap,
                                            ann->ann.getGroups().contains(AnnotationGroup.ActiveRegionBased));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));

        //skipped because it's not active region annotation
        Assert.assertFalse(resultVC.getCommonInfo().hasAttribute(GATKVCFConstants.SAMPLE_LIST_KEY));
        Assert.assertEquals(resultVC.getGenotype(0).getAD(), new int[]{ref, alt});
    }
}
