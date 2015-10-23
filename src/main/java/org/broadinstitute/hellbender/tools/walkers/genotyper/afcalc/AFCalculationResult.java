package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Describes the results of the AFCalc
 *
 * Only the bare essentials are represented here, as all AFCalc models must return meaningful results for
 * all of these fields.
 *
 * Note that all of the values -- i.e. priors -- are checked now that they are meaningful, which means
 * that users of this code can rely on the values coming out of these functions.
 */
public final class AFCalculationResult {
    private static final int AF0 = 0;
    private static final int AF1p = 1;
    private static final int LOG_ARRAY_SIZES = 2;

    private final double[] logLikelihoodsOfAC;
    private final double[] logPriorsOfAC;
    private final double[] logPosteriorsOfAC;

    private final Map<Allele, Double> logpRefByAllele;

    /**
     * The AC values for all ALT alleles at the MLE
     */
    private final int[] alleleCountsOfMLE;

    /**
     * The list of alleles actually used in computing the AF
     */
    private final List<Allele> allelesUsedInGenotyping;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     */
    public AFCalculationResult(final int[] alleleCountsOfMLE,
                               final List<Allele> allelesUsedInGenotyping,
                               final double[] logLikelihoodsOfAC,
                               final double[] logPriorsOfAC,
                               final Map<Allele, Double> logpRefByAllele) {
        Utils.nonNull(alleleCountsOfMLE, "alleleCountsOfMLE cannot be null");
        Utils.nonNull(logPriorsOfAC, "logPriorsOfAC cannot be null");
        Utils.nonNull(logLikelihoodsOfAC, "logLikelihoodsOfAC cannot be null");
        Utils.nonNull(logLikelihoodsOfAC, "logLikelihoodsOfAC cannot be null");
        Utils.nonNull(logpRefByAllele, "logpRefByAllele cannot be null");
        Utils.nonNull(allelesUsedInGenotyping, "allelesUsedInGenotyping cannot be null");
        if ( allelesUsedInGenotyping.isEmpty() ) {
            throw new IllegalArgumentException("allelesUsedInGenotyping must be non-null list of at least 1 value " + allelesUsedInGenotyping);
        }
        if ( alleleCountsOfMLE.length != allelesUsedInGenotyping.size() - 1) {
            throw new IllegalArgumentException("alleleCountsOfMLE.length " + alleleCountsOfMLE.length + " != allelesUsedInGenotyping.size() " + allelesUsedInGenotyping.size());
        }
        if ( logLikelihoodsOfAC.length != 2 ) {
            throw new IllegalArgumentException("logLikelihoodsOfAC must have length equal 2");
        }
        if ( logPriorsOfAC.length != 2 ) {
            throw new IllegalArgumentException("logPriorsOfAC must have length equal 2");
        }
        if ( logpRefByAllele.size() != allelesUsedInGenotyping.size() - 1 ) {
            throw new IllegalArgumentException("logpRefByAllele has the wrong number of elements: logpRefByAllele " + logpRefByAllele + " but allelesUsedInGenotyping " + allelesUsedInGenotyping);
        }
        if ( ! allelesUsedInGenotyping.containsAll(logpRefByAllele.keySet()) ) {
            throw new IllegalArgumentException("logpRefByAllele doesn't contain all of the alleles used in genotyping: logpRefByAllele " + logpRefByAllele + " but allelesUsedInGenotyping " + allelesUsedInGenotyping);
        }
        if ( ! MathUtils.goodLogProbVector(logLikelihoodsOfAC, LOG_ARRAY_SIZES, false) ) {
            throw new IllegalArgumentException("logLikelihoodsOfAC are bad " + Utils.join(",", logLikelihoodsOfAC));
        }
        if ( ! MathUtils.goodLogProbVector(logPriorsOfAC, LOG_ARRAY_SIZES, true) ) {
            throw new IllegalArgumentException("log10priors are bad " + Utils.join(",", logPriorsOfAC));
        }

        //make defensive copies of all arguments
        this.alleleCountsOfMLE = alleleCountsOfMLE.clone();
        this.allelesUsedInGenotyping = Collections.unmodifiableList(new ArrayList<>(allelesUsedInGenotyping));

        this.logLikelihoodsOfAC = Arrays.copyOf(logLikelihoodsOfAC, LOG_ARRAY_SIZES);
        this.logPriorsOfAC = Arrays.copyOf(logPriorsOfAC, LOG_ARRAY_SIZES);
        this.logPosteriorsOfAC = computePosteriors(logLikelihoodsOfAC, logPriorsOfAC);
        this.logpRefByAllele = Collections.unmodifiableMap(new HashMap<>(logpRefByAllele));
    }

    /**
     * Return a new AFCalcResult with a new prior probability
     *
     * @param logPriorsOfAC
     * @return
     */
    public AFCalculationResult copyWithNewPriors(final double[] logPriorsOfAC) {
        Utils.nonNull(logPriorsOfAC);
        return new AFCalculationResult(alleleCountsOfMLE, allelesUsedInGenotyping, logLikelihoodsOfAC, logPriorsOfAC, logpRefByAllele);
    }

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     *
     * This method returns a copy of the internally-stored array.
     *
     * @return a vector with allele counts, not all of which may be meaningful
     */
    public int[] getAlleleCountsOfMLE() {
        return alleleCountsOfMLE.clone();
    }

    /**
     * Returns the AC of allele a la #getAlleleCountsOfMLE
     *
     * @param allele the allele whose AC we want to know.  Error if its not in allelesUsedInGenotyping
     * @throws IllegalStateException if allele isn't in allelesUsedInGenotyping
     * @return the AC of allele
     */
    public int getAlleleCountAtMLE(final Allele allele) {
        Utils.nonNull(allele);
        return alleleCountsOfMLE[altAlleleIndex(allele)];
    }

    /**
     * Get the list of alleles actually used in genotyping.
     *
     * Due to computational / implementation constraints this may be smaller than
     * the actual list of alleles requested
     *
     * @return a non-empty list of alleles used during genotyping, the first of which is the reference allele
     */
    public List<Allele> getAllelesUsedInGenotyping() {
        return allelesUsedInGenotyping;
    }

    /**
     * Get the log normalized -- across all ACs -- posterior probability of AC == 0 for all alleles
     */
    public double getLogPosteriorOfAFEq0() {
        return logPosteriorsOfAC[AF0];
    }

    /**
     * Get the log normalized -- across all ACs -- posterior probability of AC > 0 for any alleles
     */
    public double getLogPosteriorOfAFGT0() {
        return logPosteriorsOfAC[AF1p];
    }

    /**
     * Get the log unnormalized -- across all ACs -- likelihood of AC == 0 for all alleles
     */
    public double getLogLikelihoodOfAFEq0() {
        return logLikelihoodsOfAC[AF0];
    }

    /**
     * Get the log unnormalized -- across all ACs -- likelihood of AC > 0 for any alleles
     */
    public double getLogLikelihoodOfAFGT0() {
        return logLikelihoodsOfAC[AF1p];
    }

    /**
     * Get the log unnormalized -- across all ACs -- prior probability of AC == 0 for all alleles
     */
    public double getLogPriorOfAFEq0() {
        return logPriorsOfAC[AF0];
    }

    /**
     * Get the log unnormalized -- across all ACs -- prior probability of AC > 0
     */
    public double getLogPriorOfAFGT0() {
        return logPriorsOfAC[AF1p];
    }

    @Override
    public String toString() {
        final List<String> byAllele = new LinkedList<>();
        for ( final Allele a : allelesUsedInGenotyping) {
            if (a.isNonReference()) {
                byAllele.add(String.format("%s => MLE %d / posterior %.2f", a, getAlleleCountAtMLE(a), getLogPosteriorOfAFEq0ForAllele(a)));
            }
        }
        return String.format("AFCalc%n\t\tlogPosteriorOfAFGT0=%.2f%n\t\t%s", getLogLikelihoodOfAFGT0(), Utils.join("\n\t\t", byAllele));
    }

    /**
     * Are we sufficiently confident in being non-ref that the site is considered polymorphic?
     *
     * We are non-ref if the probability of being non-ref > the emit confidence (often an argument).
     * Suppose posterior AF > 0 is exp(-5) (-5 in log space)
     * And that logMinPNonRef is -3.
     * We are considered polymorphic since exp(-5) < exp(-3) => -5 < -3
     *
     * Note that logMinPNonRef is really the minimum confidence, scaled as an error rate, so
     * if you want to be 99% confidence, then log10PNonRef should be log10(0.01) = -2.
     *
     * @param logMinPNonRef the log10 scaled min pr of being non-ref to be considered polymorphic
     *
     * @return true if there's enough confidence (relative to logMinPNonRef) to reject AF == 0
     */
    public boolean isPolymorphic(final Allele allele, final double logMinPNonRef) {
        Utils.nonNull(allele);
        return getLogPosteriorOfAFEq0ForAllele(allele) < logMinPNonRef;
    }

    /**
     * Same as #isPolymorphic but takes a phred-scaled quality score as input
     */
    public boolean isPolymorphicPhredScaledQual(final Allele allele, final double minPNonRefPhredScaledQual) {
        Utils.nonNull(allele);
        if ( minPNonRefPhredScaledQual < 0 ) {
            throw new IllegalArgumentException("phredScaledQual " + minPNonRefPhredScaledQual + " < 0 ");
        }
        final double logThreshold = MathUtils.LOG10_TO_LOG_CONVERSION * minPNonRefPhredScaledQual / -10;
        return isPolymorphic(allele, logThreshold);
    }

    /**
     * Returns the log probability that allele is not segregating
     *
     * Note that this function is p not segregating so that we can store
     * internally the log10 value of AF == 0, which grows very quickly
     * negative and yet has sufficient resolution for high confidence tests.
     * For example, if log10pRef == -100, not an unreasonably high number,
     * if we tried to store log10pNonRef we'd be looking at 1 - 10^-100, which
     * quickly underflows to 1.  So the logic here is backward from what
     * you really want (the p of segregating) but we do that for numerical
     * reasons
     *
     * Unlike the sites-level annotation, this calculation is specific to allele, and can be
     * used to separately determine how much evidence there is that allele is independently
     * segregating as opposed to the site being polymorphic with any allele.  In the bi-allelic
     * case these are obviously the same but for multiple alt alleles there can be lots of
     * evidence for one allele but not so much for any other allele
     *
     * @param allele the allele we're interested in, must be in getAllelesUsedInGenotyping
     * @return the log probability that allele is not segregating at this site
     */
    public double getLogPosteriorOfAFEq0ForAllele(final Allele allele) {
        Utils.nonNull(allele);
        final Double logpNonRef = logpRefByAllele.get(allele);
        Utils.nonNull(logpNonRef, "Unknown allele " + allele);
        return logpNonRef;
    }

    /**
     * Returns the log normalized posteriors given the log likelihoods and priors
     *
     * @param logLikelihoodsOfAC
     * @param logPriorsOfAC
     *
     * @return freshly allocated log10 normalized posteriors vector
     */
    private static double[] computePosteriors(final double[] logLikelihoodsOfAC, final double[] logPriorsOfAC) {
        final double[] logUnnormalizedPosteriors = new double[logLikelihoodsOfAC.length];
        for ( int i = 0; i < logLikelihoodsOfAC.length; i++ ) {
            logUnnormalizedPosteriors[i] = logLikelihoodsOfAC[i] + logPriorsOfAC[i];
        }
        return MathUtils.normalizeFromLog(logUnnormalizedPosteriors, true, false);
    }

    /**
     * Computes the offset into linear vectors indexed by alt allele for allele
     *
     * Things like our MLE allele count vector are indexed by alt allele index, with
     * the first alt allele being 0, the second 1, etc.  This function computes the index
     * associated with allele.
     *
     * @param allele the allele whose alt index we'd like to know
     * @throws IllegalArgumentException if allele isn't in allelesUsedInGenotyping
     * @return an index value greater than 0 suitable for indexing into the MLE and other alt allele indexed arrays
     */
    private int altAlleleIndex(final Allele allele) {
        if ( allele.isReference() ) {
            throw new IllegalArgumentException("Cannot get the alt allele index for reference allele " + allele);
        }
        final int index = allelesUsedInGenotyping.indexOf(allele);
        if ( index == -1 ) {
            throw new IllegalArgumentException("could not find allele " + allele + " in " + allelesUsedInGenotyping);
        } else {
            return index - 1;
        }
    }
}