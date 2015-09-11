package org.broadinstitute.hellbender.tools.walkers.annotator.interfaces;

/**
 * An annotation group is a set of annotation that have something in common and should be added at the same time.
 * This enum lists the known groups.
 */
public enum AnnotationGroup{
    Standard,
    ActiveRegionBased
}