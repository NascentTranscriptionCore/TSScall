#!/usr/bin/env python

# CREATED BY CHRISTOPHER LAVENDER
# BASED ON WORK BY ADAM BURKHOLDER
# INTEGRATIVE BIOINFORMATICS, NIEHS
# WORKING OBJECT ORIENTED VERSION
#
# ---------------------------------------------------------------------------
# PERFORMANCE REFACTOR (single-node). Statistical approach UNCHANGED.
# Verified to produce byte-identical BED and detail output to the original
# across bin_winner / global / FDR / no-annotation configurations, plus unit
# fuzz tests for the two rewritten hot paths.
#
# Changes, each result-preserving:
#   1. bin_winner TSS selection: O(H^2) per window -> O(H) via a two-pointer +
#      prefix sum. Same OVERLAPPING, hit-anchored sliding windows as the
#      original (NOT fixed tiling); identical totals and identical tie-breaking.
#   2. bedGraph storage: one Python dict per single-nt position (~350 B) ->
#      compact per-(strand,chromosome) int64 numpy arrays (~16-24 B). This is
#      what removes the memory wall that forced pre-filtering.
#   3. countLoci / total read count: vectorized (identical values).
#   4. Window merges in createSearchWindowsFromAnnotation and
#      createUnannotatedSearchWindowsFromBedgraph: list.pop(0) O(n^2) -> index
#      cursor O(n) (identical output order).
#
# Requires numpy. Python 2.7 or 3.
#
# TWO PRE-EXISTING BEHAVIORS ARE PRESERVED VERBATIM so that output matches the
# original exactly. Both silently discard data and are almost certainly bugs;
# they are left intact and clearly marked (search "PRESERVED LEGACY BEHAVIOR")
# so the decision to keep or fix them is explicit, not accidental:
#   (a) filterBedGraphListByWindows drops every bedGraph position ordered after
#       the last filter window (uTSS candidates past the last known TSS).
#   (b) createUnannotatedSearchWindowsFromBedgraph drops the final (highest-
#       sorted) uTSS search window.
# ---------------------------------------------------------------------------

import os
import math
import argparse
import sys
from operator import itemgetter
from array import array

import numpy as np


def writeBedHeader(file_name, description, OUTPUT):
    OUTPUT.write('track name="{}" description="{}"\n'.format(
        file_name,
        description,
    ))


# STRAND_STATUS IS USED TO DETERMINE IF STRAND IS USED IN SORT
def sortList(input_list, strand_status):
    if strand_status == 'sort_by_strand':
        return sorted(input_list, key=lambda k: (
            k['strand'],
            k['chromosome'],
            k['start']
            ))
    elif strand_status == 'ignore_strand':
        return sorted(input_list, key=lambda k: (
            k['chromosome'],
            k['start']
            ))


# ENTRY 1 IS LESS THAN ENTRY 2?
def isLessThan(entry_1, entry_2):
    for val in ['strand', 'chromosome', 'start']:
        if entry_1[val] < entry_2[val]:
            return True
        elif entry_1[val] > entry_2[val]:
            return False
    return False


# ENTRY 1 IS WITHIN ENTRY 2?
def isWithin(entry_1, entry_2):
    if entry_1['strand'] == entry_2['strand'] and\
            entry_1['chromosome'] == entry_2['chromosome']:
        if entry_1['start'] >= entry_2['start'] and\
                entry_1['end'] <= entry_2['end']:
            return True
    return False


def getID(base_name, count):
    max_entries = 999999
    feature_name = base_name + '_'
    for i in range(len(str(count)), len(str(max_entries))):
        feature_name += '0'
    feature_name += str(count)
    return feature_name


def readInReferenceAnnotation(annotation_file):
    reference_annotation = dict()
    all_gtf_keys = []
    with open(annotation_file) as f:
        for line in f:

            if not line.startswith('#'):  # Check for headers

                chromosome, source, feature, start, end, score, strand, \
                    frame, attributes = line.strip().split('\t')

                if feature == 'transcript' or feature == 'exon':

                    keys = []
                    values = []
                    gtf_fields = dict()

                    for entry in attributes.split(';')[:-1]:
                        # Check for key-value pair
                        if len(entry.split('\"')) > 1:
                            keys.append(entry.split('\"')[0].strip())
                            values.append(entry.split('\"')[1].strip())
                    for key, value in zip(keys, values):
                        gtf_fields[key] = [value]
                    for key in keys:
                        if key not in all_gtf_keys:
                            all_gtf_keys.append(key)

                    tr_id = gtf_fields.pop('transcript_id')[0]
                    gene_id = gtf_fields.pop('gene_id')[0]
                    for val in ('transcript_id', 'gene_id'):
                        all_gtf_keys.remove(val)

                    if feature == 'exon':

                        ref_id = (tr_id, chromosome)
                        if ref_id not in reference_annotation:
                            reference_annotation[ref_id] = {
                                'chromosome': chromosome,
                                'strand': strand,
                                'exons': [],
                                'gene_id': gene_id,
                                'gtf_fields': gtf_fields,
                                }
                        reference_annotation[ref_id]['exons'].append(
                            [int(start), int(end)]
                            )

    for ref_id in reference_annotation:
        t = reference_annotation[ref_id]
        # TAKE ADDITIONAL INFORMATION FROM EXON LISTS
        t['exons'].sort(key=lambda x: x[0])
        t['tr_start'] = t['exons'][0][0]
        t['tr_end'] = t['exons'][len(t['exons'])-1][1]
        if t['strand'] == '+':
            t['tss'] = t['tr_start']
        if t['strand'] == '-':
            t['tss'] = t['tr_end']
        t['gene_length'] = t['tr_end'] - t['tr_start']
        # POPULATE MISSING GTF FIELD ENTRIES
        for key in all_gtf_keys:
            if key not in t['gtf_fields']:
                t['gtf_fields'][key] = [None]
    return reference_annotation, all_gtf_keys


def _group_cmp(a_strand, a_chrom, b_strand, b_chrom):
    # Mirrors the (strand, chromosome) key ordering used by isLessThan and
    # sortList: plain Python string comparison, strand first then chromosome.
    if a_strand != b_strand:
        return -1 if a_strand < b_strand else 1
    if a_chrom != b_chrom:
        return -1 if a_chrom < b_chrom else 1
    return 0


class BedGraph(object):
    # Compact, behavior-preserving replacement for the original
    # list-of-single-nucleotide-dicts. In the original, combineAndSortBedGraphs
    # produced one dict per single-nt position
    #     {'chromosome', 'start', 'end', 'reads', 'strand'}   (start == end)
    # sorted by (strand, chromosome, start). That cost ~340 bytes/position.
    #
    # Here, positions are grouped by (strand, chromosome); each group holds two
    # parallel numpy arrays (pos, reads) in ascending-position order (~8 bytes
    # per int). Iterating self.order (sorted keys) and, within each group, the
    # pos array ascending, reproduces the EXACT global order that
    # sortList(..., 'sort_by_strand') produced on the original list, because:
    #   - the original sort key is (strand_str, chromosome_str, start_int);
    #   - sorted() over (strand, chromosome) tuple keys uses the same string
    #     comparison for the first two fields;
    #   - positions are unique within a (strand, chromosome) group (bedtools
    #     genomecov emits non-overlapping intervals per strand), so the third
    #     field is a total order with no ties.
    # 'end' is dropped because it always equalled 'start'; every original read
    # of entry['end'] is reproduced by using the position.
    def __init__(self):
        self.groups = {}   # (strand, chrom) -> {'pos': ndarray, 'reads': ndarray}
        self.order = []    # sorted list of (strand, chrom) keys

    def add_group(self, key, pos, reads):
        self.groups[key] = {'pos': pos, 'reads': reads}

    def finalize_order(self):
        self.order = sorted(self.groups.keys())

    def total_reads(self):
        # == sum(entry['reads'] for entry in bedgraph_list)
        return int(sum(int(g['reads'].sum()) for g in self.groups.values()))

    def count_loci(self, value):
        # == countLoci(bedgraph_list, value): number of single-nt positions
        # whose read count is >= value.
        loci = 0
        for g in self.groups.values():
            loci += int(np.count_nonzero(g['reads'] >= value))
        return loci


class TSSCalling(object):

    def __init__(self, **kwargs):

        self.forward_bedgraph = kwargs['forward_bedgraph']
        self.reverse_bedgraph = kwargs['reverse_bedgraph']
        self.chrom_sizes = kwargs['chrom_sizes']
        self.annotation_file = kwargs['annotation_file']
        self.output_bed = kwargs['output_bed']

        assert os.path.exists(self.forward_bedgraph)
        assert os.path.exists(self.reverse_bedgraph)
        assert os.path.exists(self.chrom_sizes)
        if self.annotation_file:
            assert os.path.exists(self.annotation_file)

        self.fdr_threshold = kwargs['fdr']
        self.false_positives = kwargs['false_positives']
        self.utss_filter_size = kwargs['utss_filter_size']
        self.utss_search_window = kwargs['utss_search_window']
        self.bidirectional_threshold = kwargs['bidirectional_threshold']
        self.cluster_threshold = kwargs['cluster_threshold']
        self.detail_file = kwargs['detail_file']
        self.cluster_bed = kwargs['cluster_bed']
        self.call_method = kwargs['call_method']
        self.annotation_join_distance = kwargs['annotation_join_distance']
        self.annotation_search_window = kwargs['annotation_search_window']
        self.bin_winner_size = kwargs['bin_winner_size']

        self.set_read_threshold = kwargs['set_read_threshold']
        try:
            int(self.set_read_threshold)
        except:
            pass
        else:
            self.set_read_threshold = int(self.set_read_threshold)

        # EVALUATE THRESHOLD METHOD ARGUMENTS; IF NONE, SET FDR_THRESHOLD
        # AT 0.001
        implied_threshold_methods = 0
        for val in [
                self.fdr_threshold,
                self.false_positives,
                self.set_read_threshold]:
            implied_threshold_methods += int(bool(val))
        if implied_threshold_methods == 1:
            pass
        elif implied_threshold_methods > 1:
            raise ValueError('More than 1 read threshold method implied!!')
        elif implied_threshold_methods == 0:
            self.fdr_threshold = 0.001

        self.tss_list = []
        self.reference_annotation = None
        self.gtf_attribute_fields = []
        self.annotated_tss_count = 0
        self.unannotated_tss_count = 0
        self.tss_cluster_count = 0
        self.unobserved_ref_count = 0

        self.execute()

    def createSearchWindowsFromAnnotation(self):
        # VALUE USED TO MERGE SEARCH WINDOWS BY PROXIMITY
        join_window = self.annotation_join_distance
        window_size = self.annotation_search_window

        current_entry = sorted(self.reference_annotation, key=lambda k: (
            self.reference_annotation[k]['strand'],
            self.reference_annotation[k]['chromosome'],
            self.reference_annotation[k]['tss'],
            # self.reference_annotation[k]['gene'],
            k,
            ))

        # POPULATE TRANSCRIPT LIST FROM SORTED LIST;
        # ADD SEARCH WINDOW EDGES TO ENTRIES
        transcript_list = []

        for ref in current_entry:
            transcript_list.append({
                'transcript_id': [ref[0]],
                'chromosome':
                    self.reference_annotation[ref]['chromosome'],
                'tss': [self.reference_annotation[ref]['tss']],
                'strand': self.reference_annotation[ref]['strand'],
                'gene_id': [self.reference_annotation[ref]['gene_id']],
                'hits': [],
                'gtf_fields': self.reference_annotation[ref]['gtf_fields'],
            })
            if self.reference_annotation[ref]['strand'] == '+':
                transcript_list[-1]['start'] = \
                    transcript_list[-1]['tss'][0] - window_size
                # MAKE SURE WINDOW END DOES NOT GO PAST TRANSCRIPT END
                end = transcript_list[-1]['tss'][0] + window_size
                if end > self.reference_annotation[ref]['tr_end']:
                    transcript_list[-1]['end'] = \
                        self.reference_annotation[ref]['tr_end']
                else:
                    transcript_list[-1]['end'] = end
            elif self.reference_annotation[ref]['strand'] == '-':
                # MAKE SURE WINDOW START DOES NOT GO PAST TRANSCRIPT START
                start = transcript_list[-1]['tss'][0] - window_size
                if start < self.reference_annotation[ref]['tr_start']:
                    transcript_list[-1]['start'] = \
                        self.reference_annotation[ref]['tr_end']
                else:
                    transcript_list[-1]['start'] = start
                transcript_list[-1]['end'] = \
                    transcript_list[-1]['tss'][0] + window_size

        merged_windows = []

        # MERGE WINDOWS BASED PROXIMITY;
        # IF WINDOWS ARE WITHIN JOIN THRESHOLD, THEY ARE MERGED;
        # IF NOT, BUT STILL OVERLAPPING, MIDPOINT BECOMES BOUNDARY
        # NOTE: index cursor instead of transcript_list.pop(0). pop(0) is O(W)
        # per call -> O(W^2) overall for W transcripts; the cursor is O(W).
        # Output order and merge logic are unchanged.
        _tl_idx = 0
        working_entry = transcript_list[_tl_idx]
        _tl_idx += 1
        while _tl_idx < len(transcript_list):
            next_entry = transcript_list[_tl_idx]
            _tl_idx += 1
            if (working_entry['strand'] == next_entry['strand']) and \
                    (working_entry['chromosome'] == next_entry['chromosome']):
                if working_entry['tss'][-1] + join_window >= \
                        next_entry['tss'][0]:
                    working_entry['transcript_id'].append(
                        next_entry['transcript_id'][0]
                        )
                    working_entry['gene_id'].append(
                        next_entry['gene_id'][0]
                        )
                    for key in working_entry['gtf_fields']:
                        working_entry['gtf_fields'][key].append(
                            next_entry['gtf_fields'][key][0]
                            )
                    # working_entry['genes'].append(next_entry['genes'][0])
                    working_entry['end'] = next_entry['end']
                    working_entry['tss'].append(next_entry['tss'][0])
                elif working_entry['end'] >= next_entry['start']:
                    working_entry['end'] = int(math.floor(
                        (working_entry['end']+next_entry['start'])/2
                        ))
                    next_entry['start'] = working_entry['end'] + 1
                    merged_windows.append(working_entry)
                    working_entry = next_entry
                else:
                    merged_windows.append(working_entry)
                    working_entry = next_entry
            else:
                merged_windows.append(working_entry)
                working_entry = next_entry
        merged_windows.append(working_entry)

        return merged_windows

    def combineAndSortBedGraphs(self, forward_bedgraph, reverse_bedgraph):
        # Behavior-preserving compact reader. Same line filtering and same
        # single-nucleotide expansion as the original readBedGraph (a line
        # spanning [start, end) contributes one position per integer in
        # range(start+1, end+1), each carrying the line's read count), but
        # positions accumulate into C-backed int64 arrays grouped by
        # (strand, chromosome) instead of one Python dict each.

        def readBedGraph(acc, bedgraph_fn, strand):
            with open(bedgraph_fn) as f:
                for line in f:
                    if not ('track' in line or line == '\n'):
                        chromosome, start, end, reads = line.strip().split()
                        start = int(start)
                        end = int(end)
                        reads = int(reads)
                        key = (strand, chromosome)
                        pair = acc.get(key)
                        if pair is None:
                            pair = [array('q'), array('q')]
                            acc[key] = pair
                        pos_arr, read_arr = pair
                        for i in range(start + 1, end + 1):
                            pos_arr.append(i)
                            read_arr.append(reads)

        acc = {}
        readBedGraph(acc, forward_bedgraph, '+')
        readBedGraph(acc, reverse_bedgraph, '-')

        bg = BedGraph()
        for key in list(acc.keys()):
            pos_arr, read_arr = acc.pop(key)
            pos_v = np.frombuffer(pos_arr, dtype=np.int64)
            read_v = np.frombuffer(read_arr, dtype=np.int64)
            # Sort within group by position (== sortList's 'start' key; strand
            # and chromosome are constant within a group). Positions are unique
            # per group, so ordering is fully determined.
            order = np.argsort(pos_v, kind='stable')
            bg.add_group(key, np.array(pos_v[order]), np.array(read_v[order]))
            del pos_arr, read_arr, pos_v, read_v, order
        bg.finalize_order()
        return bg

    # CONSIDERS TAB-DELIMITED CHROM_SIZES FILE (UCSC)
    def findGenomeSize(self, chrom_sizes):
        genome_size = 0
        with open(chrom_sizes) as f:
            for line in f:
                genome_size += int(line.strip().split()[1])
        return genome_size

    # FIND THRESHOLD FOR TSS CALLING, BASED ON
    # JOTHI ET AL. (2008) NUCLEIC ACIDS RES 36: 5221-5231.
    def findReadThreshold(self, bedgraph_list, genome_size):

        def countLoci(bedgraph_list, value):
            # Identical result to the original per-entry loop; vectorized.
            return bedgraph_list.count_loci(value)

        if self.fdr_threshold or self.false_positives:
            self.false_positives = 1
            mappable_size = 0.8 * 2 * float(genome_size)
            read_count = bedgraph_list.total_reads()
            expected_count = float(read_count)/mappable_size

            cume_probability = ((expected_count**0)/math.factorial(0)) * \
                math.exp(-expected_count)
            threshold = 1
            while True:
                probability = 1 - cume_probability
                expected_loci = probability * mappable_size
                if self.fdr_threshold:
                    observed_loci = countLoci(bedgraph_list, threshold)
                    fdr = float(expected_loci)/float(observed_loci)
                    if fdr < self.fdr_threshold:
                        return threshold
                else:
                    if expected_loci < self.false_positives:
                        return threshold
                cume_probability += \
                    ((expected_count**threshold)/math.factorial(threshold)) * \
                    math.exp(-expected_count)
                threshold += 1
        else:
            return self.set_read_threshold

    # FIND INTERSECTION WITH SEARCH_WINDOWS, BEDGRAPH_LIST;
    # HITS ARE ADDED TO WINDOW_LIST, REQUIRES SORTED LIST
    def findIntersectionWithBedGraph(self, search_windows, bedgraph_list):
        # In the original global merge-walk, a hit is recorded only when the
        # bedGraph position and the window share strand AND chromosome (isWithin
        # requires it); across (strand, chromosome) boundaries the walk merely
        # advances indices. Grouping the (already globally-sorted) windows by
        # (strand, chromosome) and walking each group's positions against its
        # windows therefore produces identical hits and identical window order.
        # search_windows is globally sorted, so per-key window lists stay sorted
        # by start.
        windows_by_key = {}
        for w in search_windows:
            windows_by_key.setdefault(
                (w['strand'], w['chromosome']), []).append(w)

        for key, win_list in windows_by_key.items():
            group = bedgraph_list.groups.get(key)
            if group is None:
                continue
            pos = group['pos']
            reads = group['reads']
            n = len(pos)
            nw = len(win_list)
            bedgraph_index = 0
            search_index = 0
            while (search_index < nw) and (bedgraph_index < n):
                w = win_list[search_index]
                p = pos[bedgraph_index]
                if w['start'] <= p <= w['end']:          # isWithin (same group)
                    w['hits'].append(
                        [int(p), int(reads[bedgraph_index])])
                    bedgraph_index += 1
                elif p < w['start']:                     # isLessThan (same group)
                    bedgraph_index += 1
                else:
                    search_index += 1

    # CREATE WINDOWS ABOUT KNOWN TSS FOR UNANNOTATED TSSs CALLING;
    # CONSIDERS ANNOTATED AND CALLED TSSs IN INSTANCE LISTS
    def createFilterWindowsFromAnnotationAndCalledTSSs(self):

        filter_windows = []

        if self.reference_annotation:
            for transcript in self.reference_annotation:
                filter_windows.append({
                    'strand': self.reference_annotation[transcript]['strand'],
                    'chromosome':
                        self.reference_annotation[transcript]['chromosome'],
                    'start': self.reference_annotation[transcript]['tss'] -
                        self.utss_filter_size,
                    'end': self.reference_annotation[transcript]['tss'] +
                        self.utss_filter_size
                    })

        if self.tss_list != []:
            for tss in self.tss_list:
                filter_windows.append({
                    'strand': tss['strand'],
                    'chromosome': tss['chromosome'],
                    'start': tss['start'] - self.utss_filter_size,
                    'end': tss['start'] + self.utss_filter_size
                    })

        return sortList(filter_windows, 'sort_by_strand')

    def filterBedGraphListByWindows(self, bedgraph_list, filter_windows):
        # FILTER BY OVERLAP WITH FILTER WINDOWS.
        # Faithful reproduction of the original single global merge-walk over
        # the sorted bedGraph vs. the sorted filter_windows, operating on the
        # grouped arrays. A single global filter-window cursor (filter_index) is
        # carried across (strand, chromosome) groups so that cross-group
        # advancement matches isLessThan exactly.
        #
        # PRESERVED LEGACY BEHAVIOR (do not "fix" without a decision):
        #   The original while-loop terminates as soon as filter_index reaches
        #   the end of filter_windows, so any bedGraph positions ordered AFTER
        #   the last filter window are dropped (they never reach working_list) --
        #   both the remainder of the current group and every later group. This
        #   silently discards uTSS candidates past the last filter window. It is
        #   reproduced verbatim here; see the '_tail_dropped' break below.
        if filter_windows == []:
            return bedgraph_list

        FW = filter_windows
        nfw = len(FW)
        filter_index = 0
        new_bg = BedGraph()

        for key in bedgraph_list.order:
            if filter_index >= nfw:
                break  # legacy tail-drop: remaining groups discarded
            strand, chrom = key
            group = bedgraph_list.groups[key]
            pos = group['pos']
            reads = group['reads']
            n = len(pos)
            keep = np.zeros(n, dtype=bool)
            bedgraph_index = 0
            while (filter_index < nfw) and (bedgraph_index < n):
                w = FW[filter_index]
                gc = _group_cmp(strand, chrom, w['strand'], w['chromosome'])
                if gc == 0:
                    p = pos[bedgraph_index]
                    if w['start'] <= p <= w['end']:      # within -> drop
                        bedgraph_index += 1
                    elif p < w['start']:                 # before window -> keep
                        keep[bedgraph_index] = True
                        bedgraph_index += 1
                    else:                                # past window -> advance
                        filter_index += 1
                elif gc < 0:                             # bg group < window group
                    keep[bedgraph_index] = True          #   isLessThan -> keep
                    bedgraph_index += 1
                else:                                    # window group < bg group
                    filter_index += 1                    #   advance window
            # positions not reached before filter windows ran out are dropped
            idx = np.nonzero(keep)[0]
            if len(idx):
                new_bg.add_group(key, pos[idx], reads[idx])

        new_bg.finalize_order()
        return new_bg

    # CREATES WINDOWS FOR UNANNOTATED TSS CALLING
    def createUnannotatedSearchWindowsFromBedgraph(self,
                                                   bedgraph_list,
                                                   read_threshold):
        # Seed one window per qualifying position (reads strictly > threshold),
        # in the same global sorted order the original produced by iterating the
        # sorted bedgraph_list. entry['end'] == entry['start'] == position, so
        # end + utss_search_window becomes position + utss_search_window.
        windows = []
        for key in bedgraph_list.order:
            strand, chromosome = key
            group = bedgraph_list.groups[key]
            pos = group['pos']
            reads = group['reads']
            sel = np.nonzero(reads > read_threshold)[0]
            for i in sel:
                p = int(pos[i])
                windows.append({
                    'strand': strand,
                    'chromosome': chromosome,
                    'start': p - self.utss_search_window,
                    'end': p + self.utss_search_window,
                    'hits': []
                    })

        # MERGE OVERLAPPING WINDOWS.
        # index cursor instead of windows.pop(0): pop(0) is O(len(windows)) per
        # call -> O(M^2) for M qualifying positions (can be very large on deep
        # data); the cursor is O(M). Merge logic and output order are unchanged.
        #
        # PRESERVED LEGACY BEHAVIOR (do not "fix" without a decision): the
        # original does NOT append the final working_entry after the loop, so
        # the last (highest-sorted) uTSS window is dropped. Reproduced verbatim.
        merged_windows = []
        if len(windows) == 0:
            # The original would raise IndexError on windows.pop(0) here; on real
            # data there is always >=1 qualifying position, so this guard only
            # changes the degenerate empty case (returns [] instead of crashing).
            return merged_windows
        _w_idx = 0
        working_entry = windows[_w_idx]
        _w_idx += 1
        while _w_idx < len(windows):
            next_entry = windows[_w_idx]
            _w_idx += 1
            if (working_entry['strand'] == next_entry['strand']) and\
                    (working_entry['chromosome'] == next_entry['chromosome']):
                if working_entry['end'] >= next_entry['start']:
                    working_entry['end'] = next_entry['end']
                else:
                    merged_windows.append(working_entry)
                    working_entry = next_entry
            else:
                merged_windows.append(working_entry)
                working_entry = next_entry
        return merged_windows

    # SORT CALLED TSSs AND ASSOCIATE INTO BIDIRECTIONAL PAIRS
    def associateBidirectionalTSSs(self):
        self.tss_list = sortList(self.tss_list, 'ignore_strand')
        for i in range(len(self.tss_list)-1):
            if self.tss_list[i]['chromosome'] == \
                    self.tss_list[i+1]['chromosome']:
                if self.tss_list[i]['strand'] == '-' and \
                        self.tss_list[i+1]['strand'] == '+':
                    if self.tss_list[i]['start'] + \
                            self.bidirectional_threshold >= \
                            self.tss_list[i+1]['start']:
                        distance = abs(self.tss_list[i]['start'] -
                                       self.tss_list[i+1]['start'])
                        self.tss_list[i]['divergent partner'] = \
                            self.tss_list[i+1]['id']
                        self.tss_list[i+1]['divergent partner'] = \
                            self.tss_list[i]['id']
                        self.tss_list[i]['divergent distance'] = distance
                        self.tss_list[i+1]['divergent distance'] = distance
                if self.tss_list[i]['strand'] == '+' and \
                        self.tss_list[i+1]['strand'] == '-':
                    if self.tss_list[i]['start'] + \
                            self.bidirectional_threshold >= \
                            self.tss_list[i+1]['start']:
                        distance = abs(self.tss_list[i]['start'] -
                                       self.tss_list[i+1]['start'])
                        self.tss_list[i]['convergent partner'] = \
                            self.tss_list[i+1]['id']
                        self.tss_list[i+1]['convergent partner'] = \
                            self.tss_list[i]['id']
                        self.tss_list[i]['convergent distance'] = distance
                        self.tss_list[i+1]['convergent distance'] = distance

    def findTSSExonIntronOverlap(self):
        exons = []
        introns = []

        if self.reference_annotation:
            for transcript in self.reference_annotation:
                for i in range(len(
                        self.reference_annotation[transcript]['exons'])):
                    strand = self.reference_annotation[transcript]['strand']
                    chromosome =\
                        self.reference_annotation[transcript]['chromosome']
                    start =\
                        self.reference_annotation[transcript]['exons'][i][0]
                    end = self.reference_annotation[transcript]['exons'][i][1]
                    exons.append({
                        'strand': strand,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end
                        })

                for i in range(
                        len(self.reference_annotation[transcript]['exons'])-1):
                    strand = self.reference_annotation[transcript]['strand']
                    chromosome =\
                        self.reference_annotation[transcript]['chromosome']
                    start = \
                        self.reference_annotation[transcript]['exons'][i][1]+1
                    end = \
                        self.reference_annotation[transcript]['exons'][i+1][0]\
                        - 1
                    introns.append({
                        'strand': strand,
                        'chromosome': chromosome,
                        'start': start,
                        'end': end
                        })

            exons = sortList(exons, 'sort_by_strand')
            introns = sortList(introns, 'sort_by_strand')
            self.tss_list = sortList(self.tss_list, 'sort_by_strand')

        def findFeatureOverlap(tss_list, feature_list, feature_key):
            if feature_list == []:
                for tss in tss_list:
                    tss[feature_key] = False
            else:
                feature_index = 0
                tss_index = 0
                while (feature_index < len(feature_list)) and\
                        (tss_index < len(tss_list)):
                    if isWithin(tss_list[tss_index],
                                feature_list[feature_index]):
                        tss_list[tss_index][feature_key] = True
                        tss_index += 1
                    else:
                        if isLessThan(tss_list[tss_index],
                                      feature_list[feature_index]):
                            tss_list[tss_index][feature_key] = False
                            tss_index += 1
                        else:
                            feature_index += 1

        findFeatureOverlap(self.tss_list, exons, 'exon_overlap')
        findFeatureOverlap(self.tss_list, introns, 'intron_overlap')

    # ASSOCIATE TSSs INTO CLUSTERS BY PROXIMITY;
    # ADD TSS CLUSTER AND NUMBER OF TSSs IN ASSOCIATED CLUSTER IN TSS ENTRY
    def associateTSSsIntoClusters(self):
        cluster_count = dict()
        self.tss_list = sortList(self.tss_list, 'ignore_strand')

        current_cluster = getID('cluster', self.tss_cluster_count)
        self.tss_cluster_count += 1
        self.tss_list[0]['cluster'] = current_cluster
        cluster_count[current_cluster] = 1

        for i in range(1, len(self.tss_list)):
            if not (self.tss_list[i-1]['chromosome'] ==
                    self.tss_list[i]['chromosome'] and
                    self.tss_list[i-1]['start'] + self.cluster_threshold >=
                    self.tss_list[i]['start']):
                current_cluster = getID('cluster', self.tss_cluster_count)
                self.tss_cluster_count += 1
            self.tss_list[i]['cluster'] = current_cluster
            if current_cluster not in cluster_count:
                cluster_count[current_cluster] = 1
            else:
                cluster_count[current_cluster] += 1

        for tss in self.tss_list:
            tss['cluster_count'] = cluster_count[tss['cluster']]

    def createDetailFile(self):
        def checkHits(window):
            for hit in window['hits']:
                if hit[1] >= self.read_threshold:
                    return True
            return False

        def writeUnobservedEntry(OUTPUT, tss, tr_ids, gene_ids, window):
            tss_id = getID('annoTSS', self.unobserved_ref_count)
            self.unobserved_ref_count += 1

            transcripts = tr_ids[0]
            genes = gene_ids[0]
            for i in range(1, len(tr_ids)):
                transcripts += ';' + tr_ids[i]
                genes += ';' + gene_ids[i]

            reads = 0
            for hit in window['hits']:
                if int(tss) == int(hit[0]):
                    reads = hit[1]

            OUTPUT.write(('{}' + '\t{}' * 15)
                         .format(
                            tss_id,
                            'unobserved reference TSS',
                            transcripts,
                            genes,
                            window['strand'],
                            window['chromosome'],
                            str(tss),
                            str(reads),
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                            'NA',
                         ))
            for key in self.gtf_attribute_fields:
                # OUTPUT.write('\t' + ';'.join(window['gtf_fields'][key]))
                OUTPUT.write('\t' + ';'.join(['None' if v is None else v for
                                              v in window['gtf_fields'][key]]))
            OUTPUT.write('\n')

        # self.findTSSExonIntronOverlap()
        # self.associateTSSsIntoClusters()
        # Remove GTF fields 'exon_number' and 'exon_id' if present
        skip_fields = ['exon_number', 'exon_id']
        for entry in skip_fields:
            if entry in self.gtf_attribute_fields:
                self.gtf_attribute_fields.remove(entry)

        with open(self.detail_file, 'w') as OUTPUT:
            OUTPUT.write(
                ('{}' + '\t{}' * 15)
                .format(
                    'TSS ID',
                    'Type',
                    'Transcripts',
                    'Gene ID',
                    'Strand',
                    'Chromosome',
                    'Position',
                    'Reads',
                    'Divergent?',
                    'Divergent partner',
                    'Divergent distance',
                    'Convergent?',
                    'Convergent partner',
                    'Convergent distance',
                    'TSS cluster',
                    'TSSs in associated cluster',
                    ))
            for field in self.gtf_attribute_fields:
                OUTPUT.write('\t' + field)
            OUTPUT.write('\n')

            for tss in self.tss_list:
                OUTPUT.write(tss['id'])
                OUTPUT.write('\t' + tss['type'])
                for key in ('transcript_id', 'gene_id'):
                    if key in tss:
                        OUTPUT.write('\t' + ';'.join(tss[key]))
                    else:
                        OUTPUT.write('\tNA')
                for entry in ['strand', 'chromosome', 'start', 'reads']:
                    OUTPUT.write('\t' + str(tss[entry]))
                if 'divergent partner' in tss:
                    OUTPUT.write('\t{}\t{}\t{}'.format(
                        'True',
                        tss['divergent partner'],
                        str(tss['divergent distance']),
                    ))
                else:
                    OUTPUT.write('\tFalse\tNA\tNA')
                if 'convergent partner' in tss:
                    OUTPUT.write('\t{}\t{}\t{}'.format(
                        'True',
                        tss['convergent partner'],
                        str(tss['convergent distance']),
                    ))
                else:
                    OUTPUT.write('\tFalse\tNA\tNA')
                # OUTPUT.write('\t' + str(
                #     tss['exon_overlap'] or tss['intron_overlap']))
                for entry in [
                        'cluster',
                        'cluster_count']:
                    OUTPUT.write('\t' + str(tss[entry]))
                if 'gtf_fields' in tss:
                    for key in self.gtf_attribute_fields:
                        # OUTPUT.write('\t' + ';'.join(tss['gtf_fields'][key]))
                        OUTPUT.write('\t' + ';'.join(
                            ['None' if v is None else
                             v for v in tss['gtf_fields'][key]]
                            ))
                else:
                    for key in self.gtf_attribute_fields:
                        OUTPUT.write('\tNA')
                OUTPUT.write('\n')

            if self.annotation_file:
                for window in self.ref_search_windows:
                    if not checkHits(window):
                        window_tss = []
                        for tr_id, gene_id, tss in zip(window['transcript_id'],
                                                       window['gene_id'],
                                                       window['tss']):
                            window_tss.append({
                                'transcript_id': tr_id,
                                'gene_id': gene_id,
                                'tss': int(tss),
                            })
                        window_tss.sort(key=itemgetter('tss'))

                        current_tss = window_tss[0]['tss']
                        current_tr_ids = [window_tss[0]['transcript_id']]
                        current_genes = [window_tss[0]['gene_id']]
                        window_index = 1
                        while window_index < len(window_tss):
                            if current_tss == window_tss[window_index]['tss']:
                                current_tr_ids.append(
                                    window_tss[window_index]['transcript_id'])
                                current_genes.append(
                                    window_tss[window_index]['gene_id'])
                            else:
                                writeUnobservedEntry(OUTPUT, current_tss,
                                                     current_tr_ids,
                                                     current_genes,
                                                     window)
                                current_tss = window_tss[window_index]['tss']
                                current_tr_ids = \
                                    [window_tss[window_index]['transcript_id']]
                                current_genes = [window_tss[0]['gene_id']]
                            window_index += 1
                        writeUnobservedEntry(OUTPUT, current_tss,
                                             current_tr_ids, current_genes,
                                             window)

    def writeClusterBed(self, tss_list, cluster_bed):
        clusters = dict()
        with open(cluster_bed, 'w') as OUTPUT:
            writeBedHeader(
                cluster_bed.split('.bed')[0],
                'TSScall clusters',
                OUTPUT,
            )
            for tss in tss_list:
                if tss['cluster'] in clusters:
                    clusters[tss['cluster']]['tss'].append(tss['start'])
                else:
                    clusters[tss['cluster']] = {
                        'chromosome': tss['chromosome'],
                        'tss': [tss['start']],
                    }
            for cluster in sorted(clusters):
                tss = sorted(clusters[cluster]['tss'])
                OUTPUT.write('{}\t{}\t{}\t{}\n'.format(
                    clusters[cluster]['chromosome'],
                    str(tss[0] - 1),
                    str(tss[-1]),
                    cluster,
                ))

    def writeBedFile(self, tss_list, output_bed):
        with open(output_bed, 'w') as OUTPUT:
            writeBedHeader(
                output_bed.split('.bed')[0],
                'TSScall TSSs',
                OUTPUT,
            )
            for tss in tss_list:
                OUTPUT.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    tss['chromosome'],
                    str(tss['start'] - 1),
                    str(tss['start']),
                    tss['id'],
                    '0',
                    tss['strand']
                    ))

    # FROM HITS IN SEARCH WINDOWS, CALL TSSs
    # COUNT IS RETURNED IN ORDER TO UPDATE INSTANCE VARIABLES
    def callTSSsFromIntersection(self, intersection, read_threshold, base_name,
                                 count, tss_type, nearest_allowed):
        def callTSS(hits, strand):
            if self.call_method == 'global':
                max_reads = float('-inf')
                max_position = None
                for hit in hits:
                    if hit[1] > max_reads:
                        max_position = hit[0]
                        max_reads = hit[1]
                    elif hit[1] == max_reads:
                        if strand == '+':
                            if hit[0] < max_position:
                                max_position = hit[0]
                        elif strand == '-':
                            if hit[0] > max_position:
                                max_position = hit[0]
                return max_position, max_reads
            if self.call_method == 'bin_winner':
                bin_size = self.bin_winner_size
                # O(H) EQUIVALENT OF THE ORIGINAL O(H^2) BINNING.
                # ORIGINAL: one OVERLAPPING bin anchored at every hit i,
                # covering [pos_i, pos_i + bin_size] (forward span, because
                # hits are position-sorted and j runs from i upward). This
                # slides a hit-anchored window; it does NOT tile the axis.
                # Two-pointer reproduces the identical bins/totals/ties.
                hits.sort(key=itemgetter(0))
                n = len(hits)
                prefix = [0] * (n + 1)
                for k in range(n):
                    prefix[k + 1] = prefix[k] + hits[k][1]
                # SELECT BIN WITH HIGHEST TOTAL READS
                # BECAUSE SORTED, WILL TAKE UPSTREAM BIN IN TIES (first strict max)
                max_bin_reads = float('-inf')
                max_bin_index = None
                win_r = None
                r = 0
                for i in range(n):
                    if r < i:
                        r = i
                    while r + 1 < n and hits[r + 1][0] - hits[i][0] <= bin_size:
                        r += 1
                    total_reads = prefix[r + 1] - prefix[i]
                    if total_reads > max_bin_reads:
                        max_bin_reads = total_reads
                        max_bin_index = i
                        win_r = r
                # GET LOCAL WINNER
                # BECAUSE SORTED, WILL TAKE UPSTREAM TSS IN TIES (first strict max)
                max_reads = float('-inf')
                max_position = None
                for k in range(max_bin_index, win_r + 1):
                    if hits[k][1] > max_reads:
                        max_position = hits[k][0]
                        max_reads = hits[k][1]
                return max_position, max_reads

        # ITERATE THROUGH WINDOWS IN INTERSECTION
        for entry in intersection:
            entry_hits = entry['hits']
            # LOOP WHILE 'HITS' IS POPULATED
            while len(entry_hits) != 0:
                # CALL A TSS
                tss_position, tss_reads = callTSS(entry_hits, entry['strand'])
                if tss_reads >= read_threshold:
                    self.tss_list.append({
                        'id': getID(base_name, count),
                        'type': tss_type,
                        'start': tss_position,
                        'end': tss_position,
                        'reads': tss_reads,
                        })
                    # IF VAL IN ENTRY, ADD TO DICT IN TSS LIST
                    for val in ['transcript_id', 'gene_id', 'strand',
                                'chromosome', 'gtf_fields']:
                        if val in entry:
                            self.tss_list[-1][val] = entry[val]
                    count += 1
                # GO THROUGH HITS, KEEP THOSE WITHIN NEAREST_ALLOWED
                temp = []
                for hit in entry_hits:
                    if abs(hit[0] - tss_position) > nearest_allowed:
                        temp.append(hit)
                entry_hits = temp
        return count

    def callTSSsFromAnnotation(self, bedgraph_list, read_threshold):
        self.ref_search_windows = self.createSearchWindowsFromAnnotation()
        self.findIntersectionWithBedGraph(self.ref_search_windows,
                                          bedgraph_list)
        self.annotated_tss_count = self.callTSSsFromIntersection(
            self.ref_search_windows,
            read_threshold,
            'obsTSS',
            self.annotated_tss_count,
            'called from reference window',
            float('inf')
            )

    def callUnannotatedTSSs(self, bedgraph_list, read_threshold):
        filter_windows = self.createFilterWindowsFromAnnotationAndCalledTSSs()
        filtered_bedgraph = self.filterBedGraphListByWindows(bedgraph_list,
                                                             filter_windows)
        unannotated_search_windows =\
            self.createUnannotatedSearchWindowsFromBedgraph(filtered_bedgraph,
                                                            read_threshold)
        self.findIntersectionWithBedGraph(unannotated_search_windows,
                                          filtered_bedgraph)
        self.unannotated_tss_count = self.callTSSsFromIntersection(
            unannotated_search_windows,
            read_threshold,
            'uTSS',
            self.unannotated_tss_count,
            'unannotated',
            self.utss_search_window
            )

    def execute(self):
        sys.stdout.write('Reading in bedGraph files...\n')
        bedgraph_list = self.combineAndSortBedGraphs(self.forward_bedgraph,
                                                     self.reverse_bedgraph)
        genome_size = self.findGenomeSize(self.chrom_sizes)
        sys.stdout.write('Calculating read threshold...\n')
        self.read_threshold = \
            self.findReadThreshold(bedgraph_list, genome_size)
        sys.stdout.write('Read threshold set to {}\n'.format(
            str(self.read_threshold)))

        if self.annotation_file:
            sys.stdout.write('Reading in annotation file...\n')
            self.reference_annotation, self.gtf_attribute_fields =\
                readInReferenceAnnotation(self.annotation_file)
            sys.stdout.write('Calling TSSs from annotation...\n')
            self.callTSSsFromAnnotation(bedgraph_list, self.read_threshold)
            sys.stdout.write('{} TSSs called from annotation\n'.format(
                str(self.annotated_tss_count)))
        sys.stdout.write('Calling unannotated TSSs...\n')
        self.callUnannotatedTSSs(bedgraph_list, self.read_threshold)
        sys.stdout.write('{} unannotated TSSs called\n'.format(
            str(self.unannotated_tss_count)))
        sys.stdout.write('Associating bidirectional TSSs...\n')
        self.associateBidirectionalTSSs()
        self.associateTSSsIntoClusters()
        if self.detail_file:
            sys.stdout.write('Creating detail file...\n')
            self.createDetailFile()
        if self.cluster_bed:
            sys.stdout.write('Creating cluster bed...\n')
            self.writeClusterBed(self.tss_list, self.cluster_bed)
        sys.stdout.write('Creating output bed...\n')
        self.writeBedFile(self.tss_list, self.output_bed)
        sys.stdout.write('TSS calling complete\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fdr', default=None, type=float,
                        help='set read threshold by FDR (FLOAT) (Default \
                        method: less than 0.001)')
    parser.add_argument('--false_positives', default=None, type=int,
                        help='set read threshold by false positive count')
    parser.add_argument('--utss_filter_size', default=750, type=int,
                        help='set uTSS filter size; any read within INTEGER \
                        of obsTSS/annoTSS is filtered prior to uTSS calling \
                        (Default: 750)')
    parser.add_argument('--utss_search_window', default=250, type=int,
                        help='set uTSS search window size to INTEGER \
                        (Default: 250)')
    parser.add_argument('--bidirectional_threshold', default=1000, type=int,
                        help='INTEGER threshold to associate bidirectional \
                        TSSs (Default: 1000)')
    parser.add_argument('--detail_file', default=None, type=str,
                        help='create a tab-delimited TXT file with details \
                        about TSS calls')
    parser.add_argument('--cluster_threshold', default=1000, type=int,
                        help='INTEGER threshold to associate TSSs into \
                        clusters (Default: 1000)')
    parser.add_argument('--annotation_file', '-a', type=str,
                        help='annotation in GTF format')
    parser.add_argument('--call_method', type=str, default='bin_winner',
                        choices=['global', 'bin_winner'],
                        help='TSS calling method to use (Default: bin_winner)')
    parser.add_argument('--annotation_join_distance', type=int, default=200,
                        help='set INTEGER distace threshold for joining search \
                        windows from annotation (Default: 200)')
    parser.add_argument('--annotation_search_window', type=int, default=1000,
                        help='set annotation search window size to INTEGER \
                        (Default: 1000)')
    parser.add_argument('--set_read_threshold', type=float, default=None,
                        help='set read threshold for TSS calling to FLOAT; do \
                        not determine threshold from data')
    parser.add_argument('--bin_winner_size', type=int, default=200,
                        help='set bin size for call method bin_winner \
                        (Default: 200)')
    parser.add_argument('--cluster_bed', type=str, default=None,
                        help='write clusters to output bed file')
    parser.add_argument('forward_bedgraph', type=str,
                        help='forward strand Start-seq bedgraph file')
    parser.add_argument('reverse_bedgraph', type=str,
                        help='reverse strand Start-seq bedgraph file')
    parser.add_argument('chrom_sizes', type=str,
                        help='standard tab-delimited chromosome sizes file')
    parser.add_argument('output_bed', type=str, help='output TSS BED file')
    args = parser.parse_args()

    TSSCalling(**vars(args))
