package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Original bi-allelic ~O(N) implementation.  Kept here for posterity and reference
 */
final class OriginalDiploidExactAFCalculator extends DiploidExactAFCalculator {

    @Override
    protected AFCalculationResult computeLogPNonRef(final VariantContext vc,
                                                    @SuppressWarnings("unused")
                                                    final int defaultPloidy,
                                                    final double[] logAlleleFrequencyPriors,
                                                    final StateTracker stateTracker) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(logAlleleFrequencyPriors, "logAlleleFrequencyPriors is null");
        Utils.nonNull(stateTracker, "stateTracker is null");

        final double[] logAlleleFrequencyLikelihoods = new double[logAlleleFrequencyPriors.length];
        final double[] logAlleleFrequencyPosteriors  = new double[logAlleleFrequencyPriors.length];
        final Pair<Integer, Integer> result = linearExact(vc, logAlleleFrequencyPriors, logAlleleFrequencyLikelihoods, logAlleleFrequencyPosteriors);
        final int lastK = result.getLeft();
        final int mleK = result.getRight();

        final double logLikelihoodAFGt0 = lastK == 0 ? MathUtils.LOG_P_OF_ZERO : MathUtils.logSumLog(logAlleleFrequencyLikelihoods, 1, lastK + 1);
        final double[] logLikelihoods = new double[]{logAlleleFrequencyLikelihoods[0], logLikelihoodAFGt0};
        final double[] logPriors = new double[]{logAlleleFrequencyPriors[0], MathUtils.logSumLog(logAlleleFrequencyPriors, 1)};
        final double[] logPosteriors = MathUtils.vectorSum(logLikelihoods, logPriors);

        final double logPRef = logPosteriors[1] > logPosteriors[0] ? MathUtils.LOG_P_OF_ZERO : 0.0;
        final Map<Allele, Double> logpRefByAllele = Collections.singletonMap(vc.getAlternateAllele(0), logPRef);

        return new AFCalculationResult(new int[]{mleK}, vc.getAlleles(),
                MathUtils.normalizeFromLog(logLikelihoods, true),
                MathUtils.normalizeFromLog(logPriors, true),
                logpRefByAllele);
    }

    /**
     * A simple data structure that holds the current, prev, and prev->prev likelihoods vectors
     * for the exact model calculation
     */
    private static final class ExactACCache {
        double[] kMinus2, kMinus1, kMinus0;

        private static double[] create(final int n) {
            return new double[n];
        }

        ExactACCache(final int n) {
            kMinus2 = create(n);
            kMinus1 = create(n);
            kMinus0 = create(n);
        }

        public void rotate() {
            final double[] tmp = kMinus2;
            kMinus2 = kMinus1;
            kMinus1 = kMinus0;
            kMinus0 = tmp;
        }

        public double[] getkMinus2() {
            return kMinus2;
        }

        public double[] getkMinus1() {
            return kMinus1;
        }

        public double[] getkMinus0() {
            return kMinus0;
        }
    }

    private static Pair<Integer, Integer> linearExact(final VariantContext vc,
                                                      final double[] logAlleleFrequencyPriors,
                                                      final double[] logAlleleFrequencyLikelihoods,
                                                      final double[] logAlleleFrequencyPosteriors) {
        final List<double[]> genotypeLikelihoods = getGLs(vc.getGenotypes(), true);
        final int numSamples = genotypeLikelihoods.size()-1;
        final int numChr = 2*numSamples;

        final ExactACCache logY = new ExactACCache(numSamples+1);
        logY.getkMinus0()[0] = 0.0; // the zero case

        double maxLogL = Double.NEGATIVE_INFINITY;
        boolean done = false;
        int lastK = -1, mleK = -1;

        for (int k=0; k <= numChr && ! done; k++ ) {
            final double[] kMinus0 = logY.getkMinus0();

            if ( k == 0 ) { // special case for k = 0
                for ( int j=1; j <= numSamples; j++ ) {
                    kMinus0[j] = kMinus0[j-1] + genotypeLikelihoods.get(j)[0];
                }
            } else { // k > 0
                final double[] kMinus1 = logY.getkMinus1();
                final double[] kMinus2 = logY.getkMinus2();

                for ( int j=1; j <= numSamples; j++ ) {
                    final double[] gl = genotypeLikelihoods.get(j);
                    final double logDenominator = MathUtils.LogCache.get(2*j) + MathUtils.LogCache.get(2*j-1);

                    double aa = Double.NEGATIVE_INFINITY;
                    double ab = Double.NEGATIVE_INFINITY;
                    if (k < 2*j-1) {
                        aa = MathUtils.LogCache.get(2 * j - k) + MathUtils.LogCache.get(2 * j - k - 1) + kMinus0[j - 1] + gl[0];
                    }

                    if (k < 2*j) {
                        ab = MathUtils.LogCache.get(2 * k) + MathUtils.LogCache.get(2 * j - k) + kMinus1[j - 1] + gl[1];
                    }

                    final double logMax;
                    if (k > 1) {
                        final double bb = MathUtils.LogCache.get(k) + MathUtils.LogCache.get(k-1) + kMinus2[j-1] + gl[2];
                        logMax = MathUtils.approximateLogSumLog(aa, ab, bb);
                    } else {
                        // we know we aren't considering the BB case, so we can use an optimized log10 function
                        logMax = MathUtils.approximateLogSumLog(aa, ab);
                    }

                    // finally, update the L(j,k) value
                    kMinus0[j] = logMax - logDenominator;
                }
            }

            // update the posteriors vector
            final double logLofK = kMinus0[numSamples];
            logAlleleFrequencyLikelihoods[k] = logLofK;
            logAlleleFrequencyPosteriors[k] = logLofK + logAlleleFrequencyPriors[k];

            // can we abort early?
            lastK = k;
            if ( logLofK > maxLogL ) {
                maxLogL = logLofK;
                mleK = k;
            }

            if ( logLofK < maxLogL - StateTracker.MAX_LOG_ERROR_TO_STOP_EARLY) {
                //if ( DEBUG ) System.out.printf("  *** breaking early k=%d log10L=%.2f maxLogL=%.2f%n", k, logLofK, maxLogL);
                done = true;
            }

            logY.rotate();
        }

        return Pair.of(lastK, mleK);
    }
}
