package org.broadinstitute.hellbender.tools.walkers.annotator.interfaces;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

/**
 * Superclass of all variant annotations.
 */
public abstract class VariantAnnotation {

    private final EnumSet<AnnotationGroup> groups;

    /**
     * Creates this annotation and marks is as a member of the given set of groups.
     * @param groups set of groups to add this annotation to, or null to not add this annotation to any group.
     */
    protected VariantAnnotation(final AnnotationGroup... groups){
        this.groups = groups == null || groups.length == 0 ? EnumSet.noneOf(AnnotationGroup.class) : EnumSet.copyOf(Arrays.asList(groups));
    }

    /**
     * Return the annotation groups this annotation is part of.
     */
    public EnumSet<AnnotationGroup> getGroups() {
        return groups;
    }

    /**
     * Return the keys
     */
    public abstract List<String> getKeyNames();
}
