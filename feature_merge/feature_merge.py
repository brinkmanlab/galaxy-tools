#!/usr/bin/env python
"""
Feature Merge - Merge annotations in GFF files.
Can merge overlapping annotations based on criteria.

Run without arguments for options help.

Nolan Woods (nolan_w@sfu.ca) 2019
Brinkman Lab, SFU
"""
import shutil

import gffutils
import sys
import os
import getopt

usage = """
Usage: feature_merge.py [-i] [-e] [-x] [-v] [-m merge|append|error|skip|replace] [-f type[,type..]].. <input1> [<input_n>..]
Accepts GFF or GTF format.
-v Print version and exit
-f Comma seperated types of features to merge. Must be terms or accessions from the SOFA sequence ontology, \"ALL\", or \"NONE\". (Can be provided more than once to specify multiple merge groups)
-i Ignore strand, merge feature regardless of strand
-x Only merge features with identical coordinates
-e Exclude component features from output
-m Merge strategy used to deal with id collisions between input files.
    merge: attributes of all features with the same primary key will be merged
    append: entry will have a unique, autoincremented primary key assigned to it (default)
    error: exception will be raised. This means you will have to edit the file yourself to fix the duplicated IDs
    skip: ignore duplicates, emitting a warning
    replace: keep last duplicate
"""[1:-1]

merge_strategies = {"merge": "merge", "append": "create_unique", "error": "error", "skip": "warning", "replace": "replace"}

# The following two functions are waiting to be integrated into the gffutils library.
# PR https://github.com/daler/gffutils/pull/130
# PR https://github.com/daler/gffutils/pull/131

def merge(self, features, exact_only=False, ignore_strand=False, ignore_featuretype=False):
    """
    Merge overlapping features together.

    Parameters
    ----------

    features : iterator of Feature instances

    exact_only : bool
        If True, will only merge features with identical coordinates;
        otherwise will merge any overlapping or contiguous features.

    ignore_strand : bool
        If True, features on multiple strands will be merged, and the final
        strand will be set to '.'

    ignore_feauretype : bool
        If True, features of multiple types will be merged, and the final
        type will be set to 'sequence_feature'.

    Returns
    -------
    A generator object that yields :class:`Feature` objects representing
    the newly merged features and a list of the features that compose it.
    """

    # To start, we create a merged feature of just the first feature.
    features = iter(features)
    try:
        feature = next(features)
    except StopIteration:
        return
    current_merged_start = feature.start
    current_merged_stop = feature.stop
    current_merged_seqid = feature.seqid
    strand = feature.strand
    frame = feature.frame
    featuretype = feature.featuretype
    feature_components = [feature]
    if exact_only:
        coordinate_criteria = lambda seqid, start, stop: start == current_merged_start and stop == current_merged_stop and seqid == current_merged_seqid
    else:
        coordinate_criteria = lambda seqid, start, stop: current_merged_start <= start and start <= current_merged_stop + 1 and seqid == current_merged_seqid

    for feature in features:
        # Does this feature start within the currently merged feature?...
        if coordinate_criteria(feature.seqid, feature.start, feature.stop) and (ignore_strand or feature.strand == strand) and (ignore_featuretype or feature.featuretype == featuretype):
            feature_components.append(feature)
            if feature.strand != strand: strand = '.'
            if feature.frame != frame: frame = '.'
            if feature.featuretype != featuretype: featuretype = "sequence_feature"
            # ...It starts within, so leave current_merged_start where it
            # is.  Does it extend any farther?
            if feature.stop >= current_merged_stop:
                # Extends further, so set a new stop position
                current_merged_stop = feature.stop
            else:
                # If feature.stop < current_merged_stop, it's completely
                # within the previous feature.  Nothing more to do.
                continue
        else:
            # The start position is outside the merged feature, so we're
            # done with the current merged feature.  Prepare for output...
            attributes={}
            for component in feature_components: attributes = gffutils.helpers.merge_attributes(component.attributes, attributes)
            yield self._feature_returner(
                seqid=current_merged_seqid,
                source=",".join(set(component.source for component in feature_components)),
                featuretype=featuretype,
                start=current_merged_start,
                end=current_merged_stop,
                score='.',
                strand=strand,
                frame=frame,
                attributes=attributes), feature_components

            # and we start a new one, initializing with this feature's
            # start and stop.
            current_merged_start = feature.start
            current_merged_stop = feature.stop
            current_merged_seqid = feature.seqid
            strand = feature.strand
            frame = feature.frame
            featuretype = feature.featuretype
            feature_components = [feature]

    attributes = {}
    for component in feature_components: attributes = gffutils.helpers.merge_attributes(component.attributes, attributes)
    yield self._feature_returner(
        seqid=current_merged_seqid,
        source = ",".join(set(component.source for component in feature_components)),
        featuretype=featuretype,
        start=current_merged_start,
        end=current_merged_stop,
        score='.',
        strand=strand,
        frame=frame,
        attributes=attributes), feature_components


def update(self, data, **kwargs):
    """
    Ripped this out of FeatureDB.update() to deal with a bug.
    Update database with features in `data`.
    If the file is empty, return rather than throw exception.

    data : str, iterable, FeatureDB instance
        If FeatureDB, all data will be used. If string, assume it's
        a filename of a GFF or GTF file.  Otherwise, assume it's an
        iterable of Feature objects.  The classes in gffutils.iterators may
        be helpful in this case.
    """
    from gffutils import create
    from gffutils import iterators

    # Handle all sorts of input
    data = iterators.DataIterator(data)

    if self.dialect['fmt'] == 'gtf':
        if 'id_spec' not in kwargs:
            kwargs['id_spec'] = {'gene': 'gene_id', 'transcript': 'transcript_id'}
        db = create._GTFDBCreator(
            data=data, dbfn=self.dbfn, dialect=self.dialect, **kwargs)
    elif self.dialect['fmt'] == 'gff3':
        if 'id_spec' not in kwargs:
            kwargs['id_spec'] = 'ID'
        db = create._GFFDBCreator(
            data=data, dbfn=self.dbfn, dialect=self.dialect, **kwargs)
    else:
        raise ValueError

    peek, data._iter = iterators.peek(data._iter, 1)
    if len(peek) == 0: return # If the file is empty then do nothing

    db._autoincrements.update(self._autoincrements)
    db._populate_from_lines(data)
    db._update_relations()
    db._finalize()
    self._autoincrements.update(db._autoincrements)

# -- Begin feature_merge specific code --

if __name__ == '__main__':
    ignore_strand = False
    ignore_featuretypes = False
    exact_only = False
    exclude_components = False
    featuretypes_groups = []
    merge_strategy = "create_unique"
    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'viexf:m:')
        for opt, val in opts:
            if opt == '-v':
                import __version
                print(__version.__versionstr__)
                exit(0)
            elif opt == '-i':
                ignore_strand = True
            elif opt == '-e':
                exclude_components = True
            elif opt == '-f':
                if val != "NONE": ignore_featuretypes = True
                if val != "ALL": featuretypes_groups.append(set(filter(None, val.split(','))))
            elif opt == '-x':
                exact_only = True
            elif opt == '-m':
                if val in merge_strategies:
                    merge_strategy = merge_strategies[val]
                else:
                    raise getopt.GetoptError("Invalid merge strategy", opt)

    except getopt.GetoptError as err:
        print("Argument error(", err.opt, "): ", err.msg, file=sys.stderr)
        args = []

    if len(args) < 1:
        print(usage, file=sys.stderr)
        exit(1)

    # Remove any empty files as GFFutils gets angry
    args = list(filter(os.path.getsize, args))
    if not len(args): exit(0)

    if not featuretypes_groups:
        featuretypes_groups.append(None)

    if ignore_featuretypes:
        if ignore_strand:
            merge_order = ('seqid', 'start', 'featuretype')
        else:
            merge_order = ('seqid', 'strand', 'start', 'featuretype')
    else:
        if ignore_strand:
            merge_order = ('seqid', 'featuretype', 'start')
        else:
            merge_order = ('seqid', 'featuretype', 'strand', 'start')


    # Load input data
    try:
        input = args[0]
        db = gffutils.create_db(input, ":memory:", merge_strategy=merge_strategy)

        for input in args[1:]:
            update(db, input, merge_strategy=merge_strategy)
    except Exception as e:
        print("Error while parsing ", input, e, file=sys.stderr)
        raise

    remaining_featuretypes = set(db.featuretypes())

    # Output header
    print("##gff-version 3")

    # Merge features per featuregroup
    for featuregroup in featuretypes_groups:
        if featuregroup:
            remaining_featuretypes -= featuregroup
        else:
            remaining_featuretypes = set()
        for merged, components in merge(db, db.all_features(featuretype=featuregroup, order_by=merge_order), exact_only, ignore_strand, ignore_featuretypes):
            # Store merged record components as hierarchy
            merged.id = "merged"
            if len(components) > 1:
                # Build a unique id for merged record
                for component in components:
                    if not component.id:
                        component.id = hex(hash(component) + 2**63)[2:]
                    merged.id += "-" + component.id

                # If id is too long, hash and encode it
                if len(merged.id) > 32:
                    merged.id = hex(hash(merged.id)+ 2**63)[2:]

                merged.attributes["ID"] = [merged.id]

                # Set merged record as component parent
                for component in components:
                    if "Parent" not in component.attributes:
                        component.attributes["Parent"] = []
                    component.attributes["Parent"].append(merged.id)

                # Set source
                merged.attributes["sources"] = merged.source.split(',')
                merged.source = "feature_merge"

            # Output components
            if len(components) == 1 or not exclude_components:
                for component in components:
                    component.attributes["ID"] = component.attributes.get('ID', [component.id])
                    print(component)

            # Output merged record if more than one component
            if len(components) > 1:
                print(merged)

    # Output any features that may have not been in the -f arguments
    if remaining_featuretypes:
        for feature in db.all_features(featuretype=remaining_featuretypes, order_by=merge_order):
            print(feature)

