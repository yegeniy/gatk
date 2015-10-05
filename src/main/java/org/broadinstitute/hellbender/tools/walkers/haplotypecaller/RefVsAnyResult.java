package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

/**
 * Holds information about a genotype call of a single sample reference vs. any non-ref event
 */
final class RefVsAnyResult {
    /**
     * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
     */
    final double[] genotypeLikelihoods;

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     * @param likelihoodCapacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihoodCapacity} is negative.
     */
    public RefVsAnyResult(final int likelihoodCapacity) {
        if (likelihoodCapacity < 0)
            throw new IllegalArgumentException("likelihood capacity is negative");
        genotypeLikelihoods = new double[likelihoodCapacity];
    }

    /**
     * AD field value for ref / non-ref
     */
    final int[] AD_Ref_Any = new int[2];

    /**
     * @return Get the DP (sum of AD values)
     */
    int getDP() { return AD_Ref_Any[0] + AD_Ref_Any[1]; }

    /**
     * Cap the het and hom var likelihood values by the hom ref likelihood.
     */
    void capByHomRefLikelihood() {
        final int likelihoodCount = genotypeLikelihoods.length;
        for (int i = 1; i < likelihoodCount; i++) {
            genotypeLikelihoods[i] = Math.min(genotypeLikelihoods[0], genotypeLikelihoods[i]);
        }
    }
}
