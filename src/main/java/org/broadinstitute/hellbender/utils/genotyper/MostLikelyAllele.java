package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Stores the most likely and second most likely alleles, along with a threshold
 * for assuming computing that a read is informative.
 *
 * If the difference between the most-likely allele and the next-most-likely allele is < INFORMATIVE_LOG_LIKELIHOOD_THRESHOLD
 * then the most likely allele is set to "no call", and isInformative will return false.  This constant can be
 * overridden simply by using one of the version of these calls that accepts informative threshold as an argument.
 *
 * For convenience, there are functions called getAlleleIfInformative that return either the most likely allele, or
 * NO_CALL if two or more alleles have likelihoods within INFORMATIVE_LOG_LIKELIHOOD_THRESHOLD of one another.
 *
 * By default empty allele maps will return NO_CALL, and allele maps with a single entry will return the
 * corresponding key
 */
public final class MostLikelyAllele {
    public static final double INFORMATIVE_LOG_LIKELIHOOD_THRESHOLD = MathUtils.log10ToLog(0.2);

    private final Allele mostLikely;
    private final Allele secondMostLikely;
    private final double logLikelihoodOfMostLikely;
    private final double logLikelihoodOfSecondBest;

    /**
     * Create a new MostLikelyAllele
     *
     * If there's a meaningful most likely allele, allele should be a real allele.  If none can be determined,
     * mostLikely should be a NO_CALL allele.
     *
     * @param mostLikely the most likely allele
     * @param secondMostLikely the most likely allele after mostLikely (may be null)
     * @param logLikelihoodOfMostLikely the log10 likelihood of the most likely allele      (or {@link Double.NEGATIVE_INFINITY} if none is available)
     * @param logLikelihoodOfSecondBest the log10 likelihood of the next most likely allele (or {@link Double.NEGATIVE_INFINITY} if none is available)
     */
    public MostLikelyAllele(final Allele mostLikely, final Allele secondMostLikely, final double logLikelihoodOfMostLikely, final double logLikelihoodOfSecondBest) {
        Utils.nonNull( mostLikely, "mostLikely allele cannot be null");
        if ( logLikelihoodOfMostLikely != Double.NEGATIVE_INFINITY && ! MathUtils.goodLogProbability(logLikelihoodOfMostLikely) ) {
            throw new IllegalArgumentException("logLikelihoodOfMostLikely must be either -Infinity or a good log10 prob but got " + logLikelihoodOfMostLikely);
        }
        if ( logLikelihoodOfSecondBest != Double.NEGATIVE_INFINITY && ! MathUtils.goodLogProbability(logLikelihoodOfSecondBest) ) {
            throw new IllegalArgumentException("logLikelihoodOfSecondBest must be either -Infinity or a good log10 prob but got " + logLikelihoodOfSecondBest);
        }
        if ( logLikelihoodOfMostLikely < logLikelihoodOfSecondBest ) {
            throw new IllegalArgumentException("logLikelihoodOfMostLikely must be <= logLikelihoodOfSecondBest but got " + logLikelihoodOfMostLikely + " vs 2nd " + logLikelihoodOfSecondBest);
        }

        this.mostLikely = mostLikely;
        this.secondMostLikely = secondMostLikely;
        this.logLikelihoodOfMostLikely = logLikelihoodOfMostLikely;
        this.logLikelihoodOfSecondBest = logLikelihoodOfSecondBest;
    }

    /**
     * Retruns the most likely allele.
     */
    public Allele getMostLikelyAllele() {
        return mostLikely;
    }

    /**
     * Retruns the second most likely allele or null if there is none.
     */
    public Allele getSecondMostLikelyAllele() {
        return secondMostLikely;
    }

    /**
     * Retruns the log10 likelihood of the most likely allele or {@link Double.NEGATIVE_INFINITY} if none is available.
     */
    public double getLogLikelihoodOfMostLikely() {
        return logLikelihoodOfMostLikely;
    }

    /**
     * Retruns the log10 likelihood of the second most likely allele or {@link Double.NEGATIVE_INFINITY} if none is available.
     */
    public double getLogLikelihoodOfSecondBest() {
        return logLikelihoodOfSecondBest;
    }

    /**
     * @see #isInformative(double) with threshold of INFORMATIVE_LOG_LIKELIHOOD_THRESHOLD
     */
    public boolean isInformative() {
        return isInformative(INFORMATIVE_LOG_LIKELIHOOD_THRESHOLD);
    }

    /**
     * Was this allele selected from an object that was specifically informative about the allele?
     *
     * The calculation that implements this is whether the likelihood of the most likely allele is larger
     * than the second most likely by at least the logThresholdForInformative
     *
     * @return true if so, false if not
     */
    public boolean isInformative(final double logThresholdForInformative) {
        return logLikelihoodOfMostLikely - logLikelihoodOfSecondBest > logThresholdForInformative;
    }

    /**
     * @see #getAlleleIfInformative(double) with threshold of INFORMATIVE_LOG_LIKELIHOOD_THRESHOLD
     */
    public Allele getAlleleIfInformative() {
        return getAlleleIfInformative(INFORMATIVE_LOG_LIKELIHOOD_THRESHOLD);
    }

    /**
     * Get the most likely allele if isInformative(logThresholdForInformative) is true, or NO_CALL otherwise
     *
     * @param logThresholdForInformative a log threshold to determine if the most likely allele was informative
     * @return a non-null allele
     */
    public Allele getAlleleIfInformative(final double logThresholdForInformative) {
        return isInformative(logThresholdForInformative) ? mostLikely : Allele.NO_CALL;
    }
}