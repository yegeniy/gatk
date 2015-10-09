package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;

/**
 * Allele frequency prior provider based on heterozygosity.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HeterozygosityAFPriorProvider extends AFPriorProvider {

    private final double heterozygosity;
    private final double logHeterozygosity;

    /**
     * Construct a new provider given the heterozygosity value.
     * @param heterozygosity must be a valid heterozygosity between larger than 0 and smaller than 1.
     * @throws IllegalArgumentException if {@code heterozygosity} is not valid one in the interval (0,1).
     */
    public HeterozygosityAFPriorProvider(final double heterozygosity) {
        if (heterozygosity <= 0) {
            throw new IllegalArgumentException("the heterozygosity must be greater than 0");
        }
        if (heterozygosity >= 1) {
            throw new IllegalArgumentException("the heterozygosity must be less than 1");
        }
        if (Double.isNaN(heterozygosity)) {
            throw new IllegalArgumentException("the heterozygosity cannot be a NaN");
        }
        this.heterozygosity = heterozygosity;
        this.logHeterozygosity = Math.log(heterozygosity);
    }

    @Override
    protected double[] buildPriors(final int totalPloidy) {
        final double[] result = new double [totalPloidy + 1];
        Arrays.fill(result, logHeterozygosity);
        result[0] = Double.NEGATIVE_INFINITY;
        MathUtils.LogCache.ensureCacheContains(totalPloidy);
        for (int i = 1; i <= totalPloidy; i++) {
            result[i] -= MathUtils.LogCache.get(i);
        }
        final double logSum = MathUtils.approximateLogSumLog(result);
        if (logSum >= 0) {
            throw new IllegalArgumentException("heterozygosity " + heterozygosity + " is too large of total ploidy " + totalPloidy);
        }
        result[0] = MathUtils.log1mexp(logSum);
        return result;
    }
}
