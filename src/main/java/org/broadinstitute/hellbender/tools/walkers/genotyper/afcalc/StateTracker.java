package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * TODO this class (+AFCalculator) is a bit messy... it seems that it combines "debugging" (unnecessarily adding CPU cost in production)
 * TODO but it also contains important part of the AF calculation state... why mix both!!!? It seems that the second
 * TODO part could be just blend into AFCalculator ... one one hand you want to reduce classes code size ... but these
 * TODO two classes code seems to be quite intertwine and makes it difficult to understand what is going on.
 * in the production setting without much need
 *
 * Keeps track of the state information during the exact model AF calculation.
 *
 * Tracks things like the MLE and MAP AC values, their corresponding likelihood and posterior
 * values, the likelihood of the AF == 0 state, and the number of evaluations needed
 * by the calculation to compute the P(AF == 0)
 */
final class StateTracker {
    private static final double VALUE_NOT_CALCULATED = Double.NEGATIVE_INFINITY;
    static final double MAX_LOG_ERROR_TO_STOP_EARLY = MathUtils.log10ToLog(6); // we want the calculation to be accurate to 1 / 10^6

    /**
     * These variables are intended to contain the MLE and MAP (and their corresponding allele counts)
     * of the site over all alternate alleles
     */
    private double logMLE;
    private double logMAP;

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     */
    private int[] alleleCountsOfMLE;
    private int[] alleleCountsOfMAP;

    /**
     * A vector of log10 likelihood values seen, for future summation.  When the size of the
     * vector is exceeed -- because we've pushed more posteriors than there's space to hold
     * -- we simply sum up the existing values, make that the first value, and continue.
     */
    private final double[] logLikelihoodsForAFGt0 = new double[LIKELIHOODS_CACHE_SIZE];
    private static final int LIKELIHOODS_CACHE_SIZE = 5000;
    private int logLikelihoodsForAFGt0CacheIndex = 0;

    /**
     * The actual sum of the likelihoods.  Null if the sum hasn't been computed yet
     */
    private Double logLikelihoodsForAFGt0Sum = null;

    /**
     * Contains the likelihood for the site's being monomorphic (i.e. AF=0 for all alternate alleles)
     */
    private double logLikelihoodOfAFzero = 0.0;

    /**
     * The list of alleles actually used in computing the AF
     */
    private List<Allele> allelesUsedInGenotyping = null;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     *
     * @param maxAltAlleles an integer >= 1
     */
    StateTracker(final int maxAltAlleles) {
        if ( maxAltAlleles < 0 ) {
            throw new IllegalArgumentException("maxAltAlleles must be >= 0, saw " + maxAltAlleles);
        }

        alleleCountsOfMLE = new int[maxAltAlleles];
        alleleCountsOfMAP = new int[maxAltAlleles];

        reset();
    }

    /**
     * Is the likelihood of configuration K too low to consider, related to the
     * maximum likelihood seen already?
     *
     * @param logLofK the log10 likelihood of the configuration we're considering analyzing
     * @return true if the configuration cannot meaningfully contribute to our likelihood sum
     */
    private boolean tooLowLikelihood(final double logLofK) {
        return logLofK < logMLE - MAX_LOG_ERROR_TO_STOP_EARLY;
    }

    /**
     * @return true iff all ACs in this object are less than or equal to their corresponding ACs in the provided set
     */
    private boolean isLowerAC(final ExactACcounts otherACs, final boolean otherACsContainsReference) {
        final int[] otherACcounts = otherACs.getCounts();

        final int firstAltAlleleIndex = otherACsContainsReference ? 1 : 0;

        for ( int i = firstAltAlleleIndex; i < otherACcounts.length; i++ ) {
            if ( alleleCountsOfMLE[i - firstAltAlleleIndex] > otherACcounts[i] ) {
                return false;
            }
        }

        return true;
    }

    /**
     * Should we stop exploring paths from ACs, given it's logLofK
     *
     * @param logLofK the logLofK of these ACs
     * @param ACs the ACs of this state
     * @param exactACcountsContainReference whether the {@code ACs} contains the reference allele count (index == 0) beside all other alternative alleles.
     * @return return true if there's no reason to continue with subpaths of AC, or false otherwise
     */
    @VisibleForTesting
    boolean abort(final double logLofK, final ExactACcounts ACs, final boolean enforceLowerACs, final boolean exactACcountsContainReference) {
        return tooLowLikelihood(logLofK) && (!enforceLowerACs || isLowerAC(ACs,exactACcountsContainReference));
    }

    @VisibleForTesting
    int[] getAlleleCountsOfMAP() {
        return alleleCountsOfMAP;
    }

    /**
     * @return the likelihoods summed across all AC values for AC > 0
     */
    private double getLogLikelihoodOfAFNotZero() {
        if ( logLikelihoodsForAFGt0Sum == null ) {
            if ( logLikelihoodsForAFGt0CacheIndex == 0 ){ // there's nothing to sum up, so make the sum equal to the smallest thing we have
                logLikelihoodsForAFGt0Sum = MathUtils.LOG_P_OF_ZERO;
            } else {
                logLikelihoodsForAFGt0Sum = MathUtils.logSumLog(logLikelihoodsForAFGt0, 0, logLikelihoodsForAFGt0CacheIndex);
            }
        }
        return logLikelihoodsForAFGt0Sum;
    }

    /**
     * @return the log likelihood of AF == 0
     */
    private double getLogLikelihoodOfAFzero() {
        return logLikelihoodOfAFzero;
    }

    /**
     * Convert this state to an corresponding AFCalcResult.
     *
     * Assumes that the values in this state have been filled in with meaningful values during the calculation.
     * For example, that the allelesUsedInGenotyping has been set, that the alleleCountsOfMLE contains meaningful
     * values, etc.
     *
     * @param logPriorsByAC the priors by AC
     *
     * @return an AFCalcResult summarizing the final results of this calculation
     */
    AFCalculationResult toAFCalculationResult(final double[] logPriorsByAC) {
        final int [] subACOfMLE = Arrays.copyOf(alleleCountsOfMLE, allelesUsedInGenotyping.size() - 1);
        final double[] logLikelihoods = MathUtils.normalizeFromLog(new double[]{getLogLikelihoodOfAFzero(), getLogLikelihoodOfAFNotZero()}, true);
        final double[] logPriors = MathUtils.normalizeFromLog(new double[]{logPriorsByAC[0], MathUtils.logSumLog(logPriorsByAC, 1)}, true);

        final Map<Allele, Double> logpRefByAllele = new HashMap<>(allelesUsedInGenotyping.size());
        for ( int i = 0; i < subACOfMLE.length; i++ ) {
            final Allele allele = allelesUsedInGenotyping.get(i+1);
            final double logPRef = alleleCountsOfMAP[i] > 0 ? -10000 : 0; // TODO -- a total hack but in effect what the old behavior was
            logpRefByAllele.put(allele, logPRef);
        }

        return new AFCalculationResult(subACOfMLE, allelesUsedInGenotyping, logLikelihoods, logPriors, logpRefByAllele);
    }

    // --------------------------------------------------------------------------------
    //
    // Protected mutational methods only for use within the calculation models themselves
    //
    // --------------------------------------------------------------------------------

    /**
     * Reset the data in this results object, so that it can be used in a subsequent AF calculation
     *
     * Resetting of the data is done by the calculation model itself, so shouldn't be done by callers any longer
     *
     * @param ensureAltAlleleCapacity indicate the minimum number of alt-alleles that should be supported by the
     *                                tracker.
     */
    void reset(final int ensureAltAlleleCapacity) {
        logMLE = logMAP = logLikelihoodOfAFzero = VALUE_NOT_CALCULATED;
        logLikelihoodsForAFGt0CacheIndex = 0;
        logLikelihoodsForAFGt0Sum = null;
        allelesUsedInGenotyping = null;
        if (alleleCountsOfMAP.length < ensureAltAlleleCapacity) {
            final int newCapacity = Math.max(ensureAltAlleleCapacity, alleleCountsOfMAP.length << 1);
            alleleCountsOfMAP = new int[newCapacity];
            alleleCountsOfMLE = new int[newCapacity];
        } else {
            Arrays.fill(alleleCountsOfMLE, 0);
            Arrays.fill(alleleCountsOfMAP, 0);
        }
        Arrays.fill(logLikelihoodsForAFGt0, Double.POSITIVE_INFINITY);
    }

    /**
     * Reset the data in this results object, so that it can be used in a subsequent AF calculation
     *
     * Resetting of the data is done by the calculation model itself, so shouldn't be done by callers any longer
     */
    void reset() {
        logMLE = logMAP = logLikelihoodOfAFzero = VALUE_NOT_CALCULATED;
        logLikelihoodsForAFGt0CacheIndex = 0;
        logLikelihoodsForAFGt0Sum = null;
        allelesUsedInGenotyping = null;
        Arrays.fill(alleleCountsOfMLE, 0);
        Arrays.fill(alleleCountsOfMAP, 0);
        Arrays.fill(logLikelihoodsForAFGt0, Double.POSITIVE_INFINITY);
    }

    /**
     * Update the maximum log likelihoods seen, if logLofKs is higher, and the corresponding ACs of this state
     *
     * @param logLofK the likelihood of our current configuration state, cannot be the 0 state
     * @param alleleCountsForK the allele counts for this state
     */
    void updateMLEifNeeded(final double logLofK, final int[] alleleCountsForK) {
        addToLikelihoodsCache(logLofK);

        if ( logLofK > logMLE) {
            logMLE = logLofK;
            System.arraycopy(alleleCountsForK, 0, alleleCountsOfMLE, 0, alleleCountsForK.length);
        }
    }

    /**
     * Update the maximum log10 posterior seen, if log10PofKs is higher, and the corresponding ACs of this state
     *
     * @param logPofK the posterior of our current configuration state
     * @param alleleCountsForK the allele counts for this state
     */
    void updateMAPifNeeded(final double logPofK, final int[] alleleCountsForK) {
        if ( logPofK > logMAP) {
            logMAP = logPofK;
            System.arraycopy(alleleCountsForK, 0, alleleCountsOfMAP, 0, alleleCountsForK.length);
        }
    }

    private void addToLikelihoodsCache(final double logLofK) {
        // add to the cache
        logLikelihoodsForAFGt0[logLikelihoodsForAFGt0CacheIndex++] = logLofK;

        // if we've filled up the cache, then condense by summing up all of the values and placing the sum back into the first cell
        if ( logLikelihoodsForAFGt0CacheIndex == LIKELIHOODS_CACHE_SIZE) {
            final double temporarySum = MathUtils.logSumLog(logLikelihoodsForAFGt0, 0, logLikelihoodsForAFGt0CacheIndex);
            Arrays.fill(logLikelihoodsForAFGt0, Double.POSITIVE_INFINITY);
            logLikelihoodsForAFGt0[0] = temporarySum;
            logLikelihoodsForAFGt0CacheIndex = 1;
        }
    }

    void setLogLikelihoodOfAFzero(final double logLikelihoodOfAFzero) {
        this.logLikelihoodOfAFzero = logLikelihoodOfAFzero;
        if ( logLikelihoodOfAFzero > logMLE) {
            logMLE = logLikelihoodOfAFzero;
            Arrays.fill(alleleCountsOfMLE, 0);
        }
    }

    void setLogPosteriorOfAFzero(final double logPosteriorOfAFzero) {
        if ( logPosteriorOfAFzero > logMAP) {
            logMAP = logPosteriorOfAFzero;
            Arrays.fill(alleleCountsOfMAP, 0);
        }
    }

    /**
     * Set the list of alleles used in genotyping
     *
     * @param allelesUsedInGenotyping the list of alleles, where the first allele is reference
     */
    void setAllelesUsedInGenotyping(final List<Allele> allelesUsedInGenotyping) {
        if ( allelesUsedInGenotyping == null || allelesUsedInGenotyping.isEmpty() ) {
            throw new IllegalArgumentException("allelesUsedInGenotyping cannot be null or empty");
        }
        if ( allelesUsedInGenotyping.get(0).isNonReference() ) {
            throw new IllegalArgumentException("The first element of allelesUsedInGenotyping must be the reference allele");
        }

        this.allelesUsedInGenotyping = allelesUsedInGenotyping;
    }

    public void ensureMaximumAlleleCapacity(final int capacity) {
        if (this.alleleCountsOfMAP.length < capacity) {
            reset(capacity);
        }
    }
}
