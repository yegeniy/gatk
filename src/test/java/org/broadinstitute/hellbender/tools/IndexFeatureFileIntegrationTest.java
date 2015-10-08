package org.broadinstitute.hellbender.tools;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.linear.LinearIndex;
import htsjdk.tribble.index.tabix.TabixIndex;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public final class IndexFeatureFileIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testVCFIndex() {
        final File ORIG_FILE = new File(getTestDataDir(), "test_variants_for_index.vcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.", ".idx");

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        for (int chr = 1; chr <= 4; chr++) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)));
        }
        for (int chr = 5; chr <= 22; chr++) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr));
        }

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
    }

    @Test
    public void testVCFIndex_inferredName() {
        final File ORIG_FILE = new File(getTestDataDir(), "test_variants_for_index.vcf");

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
        final File tribbleIndex = Tribble.indexFile(ORIG_FILE);
        Assert.assertEquals(res, tribbleIndex.getAbsolutePath());
        tribbleIndex.deleteOnExit();

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        for (int chr = 1; chr <= 4; chr++) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)));
        }
        for (int chr = 5; chr <= 22; chr++) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr));
        }

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
    }

    @Test
    public void testVCFGZIndex() {
        final File ORIG_FILE = new File(getTestDataDir(), "test_variants_for_index.vcf.gz");
        final File outName = BaseTest.createTempFile("test_variants_for_index.", ".idx");

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof TabixIndex);
        for (int chr = 1; chr <= 4; chr++) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)));
        }
        for (int chr = 5; chr <= 22; chr++) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr));
        }

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "3", "4"));
    }

    @Test(enabled = false)
    public void testVCFGZIndex_inferredName() throws IOException {
        final File ORIG_FILE = new File(getTestDataDir(), "test_variants_for_index.vcf.gz");
        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
        };
        final Object res = this.runCommandLine(args);
        final File tabixIndex = new File(ORIG_FILE.getAbsolutePath() + IndexFeatureFile.TABIX_INDEX_EXTENSION);;
        Assert.assertEquals(res, tabixIndex.getAbsolutePath());
        tabixIndex.deleteOnExit();


        final Index index0 = IndexFactory.loadIndex(res.toString()); //BUG <- this blows up
        final Index index1 = new TabixIndex(tabixIndex);             //BUG <- this blows up
    }

    @Test
    public void testBCFIndex() {
        final File ORIG_FILE = new File(getTestDataDir(), "test_variants_for_index.bcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.", ".idx");

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1"));
        for (int chr = 1; chr <= 1; chr++) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)));
        }
        for (int chr = 2; chr <= 22; chr++) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr));
        }
    }

    @Test
    public void testGVCF_VCFIndex() {
        final File ORIG_FILE = new File(getTestDataDir(), "test_variants_for_index.gvcf.vcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.gvcf.", ".idx");

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        for (int chr = 1; chr <= 1; chr++) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)));
        }
        for (int chr = 2; chr <= 22; chr++) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr));
        }

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1"));
    }

    @Test
    public void testGVCFIndex() {
        final File ORIG_FILE = new File(getTestDataDir(), "test_variants_for_index.gvcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.gvcf.", ".idx");

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        for (int chr = 1; chr <= 1; chr++) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)));
        }
        for (int chr = 2; chr <= 22; chr++) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr));
        }

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1"));
    }

    @Test
    public void testBedIndex() {
        final File ORIG_FILE = new File(getTestDataDir(), "test_bed_for_index.bed");
        final File outName = BaseTest.createTempFile("test_bed_for_index", ".idx");

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
        Assert.assertEquals(res, outName.getAbsolutePath());

        final Index index = IndexFactory.loadIndex(res.toString());
        Assert.assertTrue(index instanceof LinearIndex);
        for (int chr: new int[]{1,2,4}) {  //note: unusual loop
            Assert.assertTrue(index.containsChromosome(String.valueOf(chr)), String.valueOf(chr));
        }
        for (int chr: new int[]{3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}) { //note: unusual loop
            Assert.assertFalse(index.containsChromosome(String.valueOf(chr)), String.valueOf(chr));
        }
        for (final String chr : Arrays.asList("X", "Y", "MT")){
            Assert.assertFalse(index.containsChromosome(chr), String.valueOf(chr));
        }

        Assert.assertEquals(index.getSequenceNames(), Arrays.asList("1", "2", "4"));
    }


    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testVCFIndex_missingFile() {
        final File ORIG_FILE = new File(getTestDataDir(), "missing_file.vcf");
        final File outName = BaseTest.createTempFile("test_variants_for_index.", ".idx");

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }

    @Test(expectedExceptions = UserException.CouldNotCreateOutputFile.class)
    public void testVCFIndex_cannotWrite() {
        final File ORIG_FILE = new File(getTestDataDir(), "test_variants_for_index.vcf");
        final File outName = new File("/etc/fred.txt");  //we can't write to this

        final String[] args = new String[]{
                "--feature_file" ,  ORIG_FILE.getAbsolutePath(),
                "-O" ,  outName.getAbsolutePath()
        };
        final Object res = this.runCommandLine(args);
    }
}
