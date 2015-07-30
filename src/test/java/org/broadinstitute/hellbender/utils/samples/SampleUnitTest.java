package org.broadinstitute.hellbender.utils.samples;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import org.broadinstitute.hellbender.utils.test.BaseTest;

public class SampleUnitTest extends BaseTest {

    @DataProvider(name="basicSamples")
    public Object[][] basicSamples() {
        return new Object[][] {
                { new Sample("1C", "fam1", "1M", "1F", Sex.UNKNOWN), "1C", "fam1", "1M", "1F", Sex.UNKNOWN, Affection.UNKNOWN },
                { new Sample("1F", "fam1", null, null, Sex.MALE), "1F", "fam1", null, null, Sex.MALE, Affection.UNKNOWN },
                { new Sample("1M", "fam1", null, null, Sex.FEMALE), "1M", "fam1", null, null, Sex.FEMALE, Affection.UNKNOWN },

                // Samples with Affection
                { new Sample("1C", "fam1", "1M", "1F", Sex.UNKNOWN, Affection.AFFECTED), "1C", "fam1", "1M", "1F", Sex.UNKNOWN, Affection.AFFECTED},
                { new Sample("1F", null, null, null, Sex.MALE, Affection.UNAFFECTED), "1F", null, null, null, Sex.MALE, Affection.UNAFFECTED },
                { new Sample("1M", null, null, null, Sex.FEMALE, Affection.OTHER), "1M", null, null, null, Sex.FEMALE, Affection.OTHER }
        };
    }

    /**
     * Basic getters
     */
    @Test(dataProvider="basicSamples")
    public void basicSampleTest(Sample sample, String id, String famID, String paternalID, String maternalID, Sex gender, Affection affection) {
        Assert.assertTrue(id.equals(sample.getID()));
        Assert.assertTrue(famID == null || famID.equals(sample.getFamilyID()));
        Assert.assertTrue(maternalID == null || maternalID.equals(sample.getMaternalID()));
        Assert.assertTrue(paternalID == null || paternalID.equals(sample.getPaternalID()));
        Assert.assertEquals(gender, sample.getSex());
        Assert.assertEquals(affection, sample.getAffection());
    }

    @Test(dataProvider="basicSamples")
    public void testMergeSamples(Sample sample, String id, String famID, String paternalID, String maternalID, Sex gender, Affection affection) {

        Sample newSample = new Sample("newSample", null, null, null, Sex.UNKNOWN, Affection.UNKNOWN);
        Sample mergedSample1 = newSample.mergeSamples(sample);
        Assert.assertTrue(mergedSample1.getID().equals("newSample"));

        if (famID == null) {
            Assert.assertEquals(null, mergedSample1.getFamilyID());
        }
        else {
            Assert.assertTrue(famID.equals(mergedSample1.getFamilyID()));
        }

        if (maternalID == null) {
            Assert.assertEquals(null, mergedSample1.getMaternalID());
        }
        else {
            Assert.assertTrue(maternalID.equals(mergedSample1.getMaternalID()));
        }

        if (paternalID == null) {
            Assert.assertEquals(null, mergedSample1.getPaternalID());
        }
        else {
            Assert.assertTrue(paternalID.equals(mergedSample1.getPaternalID()));
        }

        Assert.assertEquals(mergedSample1.getSex(), gender);
        Assert.assertEquals(mergedSample1.getAffection(), affection);
    }
}
