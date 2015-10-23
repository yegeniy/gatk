package org.broadinstitute.hellbender.tools.walkers.genotyper;

/*
* Copyright (c) 2015 The Broad Institute
*
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/**
 * Created by davidben on 10/22/15.
 * This class is a lightweight wrapper to htsjdk's genotypeLikelihhods class.  The latter is log10-based, while gatk is natural log-based.
 * We have only wrapped public methods that are used in the gatk.  We did not wrap static methods that don't involve likelihoods.
 *
 * We only need to change:
 * 1) the constructor, which takes a log-likelihoods argument
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.util.List;



public class GenotypeLikelihoodsWrapper {
    private final GenotypeLikelihoods genotypeLikelihoods;


    /**
     * Construct from logLikelihoods.  Note that we don't have to convert to log10 before passing on the argument
     * to htsjdk's GenotypeLikelihoods because the htsjdk doen't do anything other than store them.  Thus this wrapper class
     * store's a copy of GenotypeLikelihoods that thinks it's storing log10Likelihoods but is actually storing natural log-likelihoods.
     * @param logLikelihoods
     */
    public GenotypeLikelihoodsWrapper(double[] logLikelihoods) {
        genotypeLikelihoods = GenotypeLikelihoods.fromLog10Likelihoods(logLikelihoods);
    }

    public final static GenotypeLikelihoodsWrapper fromLogLikelihoods(double[] logLikelihoods) {
        return new GenotypeLikelihoodsWrapper(logLikelihoods);
    }

    public final static GenotypeLikelihoodsWrapper fromPLs(final int[] pls) {
        return new GenotypeLikelihoodsWrapper(PLsToGLs(pls));
    }

    public final static GenotypeLikelihoodsWrapper fromGenotype(final Genotype genotype) {
        return genotype.hasLikelihoods() ? fromPLs(genotype.getPL()) : null;
    }


    public double[] getAsVector() {
        return genotypeLikelihoods.getAsVector();
    }

    public String toString() {
        return getAsString();
    }

    public String getAsString() {
        return genotypeLikelihoods.getAsString();
    }

    @Override public boolean equals(Object aThat) {
        //check for self-comparison
        if ( this == aThat ) return true;

        if ( !(aThat instanceof GenotypeLikelihoodsWrapper) ) return false;
        GenotypeLikelihoodsWrapper that = (GenotypeLikelihoodsWrapper) aThat;

        return this.genotypeLikelihoods.equals(that.genotypeLikelihoods);
    }

    public double getLogGQ(Genotype genotype, List<Allele> vcAlleles ) {
        return MathUtils.log10ToLog(genotypeLikelihoods.getLog10GQ(genotype, vcAlleles));
    }

    public double getLogGQ(Genotype genotype, VariantContext context) {
        return MathUtils.log10ToLog(genotypeLikelihoods.getLog10GQ(genotype, context));
    }

    public static double getLogGQFromLikelihoods(int iOfChoosenGenotype, double[] likelihoods){
       return GenotypeLikelihoods.getGQLog10FromLikelihoods(iOfChoosenGenotype, likelihoods);
    }

    //Same as htsjdk GenotypeLikelihoods method, but converts PL to natural log GL, not log10
    private final static double[] PLsToGLs(final int pls[]) {
        double[] likelihoodsAsVector = new double[pls.length];
        for ( int i = 0; i < pls.length; i++ ) {
            likelihoodsAsVector[i] = pls[i] * QualityUtils.PHRED_TO_LOG_PROB_MULTIPLIER;
        }
        return likelihoodsAsVector;
    }
}
