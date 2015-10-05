package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

public final class HaplotypeCaller {
    /**
     * Create a context that maps each read to the reference haplotype with log10 L of 0
     * @param refHaplotype a non-null reference haplotype
     * @param samples a list of all samples
     * @param region the active region containing reads
     * @return a map from sample -> PerReadAlleleLikelihoodMap that maps each read to ref
     */
    public static ReadLikelihoods<Haplotype> createDummyStratifiedReadMap(final Haplotype refHaplotype,
                                                                          final SampleList samples,
                                                                          final AssemblyRegion region) {
        return new ReadLikelihoods<>(samples, new IndexedAlleleList<>(refHaplotype),
                splitReadsBySample(samples, region.getReads()));
    }
    public static Map<String, List<GATKRead>> splitReadsBySample( final SampleList samplesList, final Collection<GATKRead> reads ) {
        final Map<String, List<GATKRead>> returnMap = new HashMap<>();
        final int sampleCount = samplesList.numberOfSamples();
        for (int i = 0; i < sampleCount; i++)
            returnMap.put(samplesList.getSample(i), new ArrayList<>());

        for( final GATKRead read : reads )
            returnMap.get(read.getReadGroup().getSample()).add(read);

        return returnMap;
    }
}
