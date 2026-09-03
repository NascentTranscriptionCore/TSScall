#!/usr/bin/env python3
"""
analyze_overlap.py -- comparative peak analytics for TSScall BED outputs.

Answers, for any pair of TSScall BED files:
  * how many called TSSs sit at IDENTICAL coordinates,
  * how many are the same peak SHIFTED by a few bp (and by how much),
  * how many are UNIQUE to one side (genuine gain/loss),
broken down by TSS type (obsTSS = annotated, uTSS = unannotated/enhancer-class).

Run over a results_<EXPNAME> tree produced by the benchmark harness. It walks
conditions/<script>_preFilt<PF>_callThresh<CT>/<EXPNAME>_output.bed and runs
three comparison families:
  1. preFilt effect at matched callThresh: reference preFilt (default the
     lowest present, i.e. no pre-filter) vs each other preFilt level.
     -> speaks directly to "is it reasonable to standardize on no pre-filter?"
  2. callThresh effect at fixed preFilt: adjacent calling thresholds.
  3. ORIG vs REFAC at every cell (coordinate-level identity check).

Outputs (under the results dir):
  overlap_summary.tsv     one row per comparison, aggregate metrics
  overlap_by_type.tsv     same, split into obsTSS / uTSS
  shift_histogram.tsv     distance buckets per comparison
and prints a readable digest with a standardization verdict.

Pure stdlib + numpy. No bedtools required.
"""
import argparse
import os
import sys
import glob
from collections import defaultdict, OrderedDict

import numpy as np


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #
def read_bed(path):
    """Return dict[(chrom, strand)] -> (pos_sorted[int64], type_sorted[U8]).
    TSScall BED row: chrom, start0, start1, id, score, strand.
    Peak coordinate = start1 (col 3). Type = id prefix before '_'."""
    groups_pos = defaultdict(list)
    groups_typ = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('track') or line.startswith('browser'):
                continue
            if not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 6:
                continue
            chrom, pos, name, strand = c[0], int(c[2]), c[3], c[5]
            ttype = name.split('_', 1)[0]
            key = (chrom, strand)
            groups_pos[key].append(pos)
            groups_typ[key].append(ttype)
    out = {}
    for key in groups_pos:
        p = np.asarray(groups_pos[key], dtype=np.int64)
        t = np.asarray(groups_typ[key], dtype=object)
        order = np.argsort(p, kind='stable')
        out[key] = (p[order], t[order])
    return out


def total_peaks(bed):
    return int(sum(len(p) for p, _ in bed.values()))


def type_counts(bed):
    d = defaultdict(int)
    for _, t in bed.values():
        for x in t:
            d[x] += 1
    return dict(d)


# --------------------------------------------------------------------------- #
# matching
# --------------------------------------------------------------------------- #
def match_pair(A, B, max_shift):
    """Match peaks between two BEDs within each (chrom, strand) group.

    Strategy: exact matches first (dist 0), then greedy nearest pairing on the
    leftovers, accepting the smallest distances first so each peak is used once.

    Returns:
      matched  : list of (dist, signed, a_type, b_type, chrom, strand, a_pos, b_pos)
      uniqueA  : list of (a_type, chrom, strand, pos)
      uniqueB  : list of (b_type, chrom, strand, pos)"""
    matched = []
    uniqueA = []
    uniqueB = []
    keys = set(A) | set(B)
    for key in keys:
        chrom, strand = key
        ap, at = A.get(key, (np.empty(0, np.int64), np.empty(0, object)))
        bp, bt = B.get(key, (np.empty(0, np.int64), np.empty(0, object)))
        na, nb = len(ap), len(bp)
        a_used = np.zeros(na, dtype=bool)
        b_used = np.zeros(nb, dtype=bool)

        # 1) exact matches via a position dict on B
        if na and nb:
            bpos_index = {}
            for j, pv in enumerate(bp):
                bpos_index.setdefault(int(pv), []).append(j)
            for i in range(na):
                lst = bpos_index.get(int(ap[i]))
                if lst:
                    j = lst.pop()          # one-to-one
                    a_used[i] = True
                    b_used[j] = True
                    matched.append((0, 0, at[i], bt[j], chrom, strand,
                                    int(ap[i]), int(bp[j])))

        # 2) greedy nearest on leftovers within max_shift
        rem_a = np.nonzero(~a_used)[0]
        rem_b = np.nonzero(~b_used)[0]
        if len(rem_a) and len(rem_b) and max_shift > 0:
            bpos = bp[rem_b]
            order_b = np.argsort(bpos, kind='stable')
            bpos_s = bpos[order_b]
            cand = []             # (dist, i_local, j_local)
            for ia, i in enumerate(rem_a):
                p = ap[i]
                k = np.searchsorted(bpos_s, p)
                for kk in (k - 1, k):
                    if 0 <= kk < len(bpos_s):
                        d = int(abs(int(p) - int(bpos_s[kk])))
                        if d <= max_shift:
                            cand.append((d, ia, kk))
            cand.sort(key=lambda x: x[0])
            la_used = np.zeros(len(rem_a), dtype=bool)
            lb_used = np.zeros(len(bpos_s), dtype=bool)
            for d, ia, kk in cand:
                if la_used[ia] or lb_used[kk]:
                    continue
                la_used[ia] = True
                lb_used[kk] = True
                i = rem_a[ia]
                j = rem_b[order_b[kk]]
                signed = int(bp[j]) - int(ap[i])
                matched.append((d, signed, at[i], bt[j], chrom, strand,
                                int(ap[i]), int(bp[j])))
                a_used[i] = True
                b_used[j] = True

        for i in np.nonzero(~a_used)[0]:
            uniqueA.append((at[i], chrom, strand, int(ap[i])))
        for j in np.nonzero(~b_used)[0]:
            uniqueB.append((bt[j], chrom, strand, int(bp[j])))
    return matched, uniqueA, uniqueB


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #
BUCKETS = [(0, 0, 'exact'), (1, 10, '1-10'), (11, 50, '11-50'),
           (51, 100, '51-100'), (101, 500, '101-500'), (501, 10**9, '501+')]


def summarize(name, A, B, matched, uniqueA, uniqueB):
    nA, nB = total_peaks(A), total_peaks(B)
    dists = np.asarray([m[0] for m in matched], dtype=np.int64)
    signed = np.asarray([m[1] for m in matched], dtype=np.int64)
    n_exact = int(np.count_nonzero(dists == 0))
    n_shift = int(np.count_nonzero(dists > 0))
    n_uA, n_uB = len(uniqueA), len(uniqueB)
    union = nA + nB - (n_exact + n_shift)
    jacc_exact = n_exact / (nA + nB - n_exact) if (nA + nB - n_exact) else 1.0
    jacc_w = (n_exact + n_shift) / union if union else 1.0
    row = OrderedDict()
    row['comparison'] = name
    row['n_A'] = nA
    row['n_B'] = nB
    row['exact'] = n_exact
    row['shifted'] = n_shift
    row['unique_A'] = n_uA
    row['unique_B'] = n_uB
    row['jaccard_exact'] = round(jacc_exact, 4)
    row['jaccard_within_window'] = round(jacc_w, 4)
    row['median_shift_bp'] = int(np.median(signed[dists > 0])) if n_shift else 0
    row['mean_abs_shift_bp'] = round(float(dists[dists > 0].mean()), 1) if n_shift else 0.0
    return row, dists


def _write_bed(path, records):
    # records: iterable of (chrom, strand, pos, name)
    with open(path, 'w') as f:
        for chrom, strand, pos, name in records:
            f.write('%s\t%d\t%d\t%s\t0\t%s\n'
                    % (chrom, pos - 1, pos, name, strand))


def dump_comparison(dump_dir, name, matched, uniqueA, uniqueB):
    """Write browser-loadable BEDs of what changed between A and B."""
    safe = name.replace(':', '__').replace('@', '_at_')
    d = os.path.join(dump_dir, safe)
    os.makedirs(d, exist_ok=True)
    # lost = in A (reference) but not B ; gained = in B but not A
    _write_bed(os.path.join(d, 'lost_from_reference.bed'),
               ((c, s, p, t) for (t, c, s, p) in uniqueA))
    _write_bed(os.path.join(d, 'gained_in_other.bed'),
               ((c, s, p, t) for (t, c, s, p) in uniqueB))
    # shifted: reference-side position, labelled with the signed shift
    _write_bed(os.path.join(d, 'shifted_reference_side.bed'),
               ((c, s, ap, '%s_shift%+d' % (at, sg))
                for (dist, sg, at, bt, c, s, ap, bp) in matched if dist > 0))


def bucket_row(name, dists, n_uA, n_uB):
    row = OrderedDict()
    row['comparison'] = name
    for lo, hi, label in BUCKETS:
        row[label] = int(np.count_nonzero((dists >= lo) & (dists <= hi)))
    row['unique_A'] = n_uA
    row['unique_B'] = n_uB
    return row


def type_breakdown(name, matched, uniqueA, uniqueB):
    rows = []
    for ttype in ('obsTSS', 'uTSS'):
        ex = sum(1 for m in matched if m[0] == 0 and m[2] == ttype)
        sh = sum(1 for m in matched if m[0] > 0 and m[2] == ttype)
        uA = sum(1 for t in uniqueA if t[0] == ttype)
        uB = sum(1 for t in uniqueB if t[0] == ttype)
        rows.append(OrderedDict(comparison=name, tss_type=ttype,
                                exact=ex, shifted=sh,
                                unique_A=uA, unique_B=uB))
    return rows


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def bed_path(results, expname, tag, pf, ct):
    return os.path.join(results, 'conditions',
                        '%s_preFilt%s_callThresh%s' % (tag, pf, ct),
                        '%s_output.bed' % expname)


def discover(results, expname, tag):
    pfs, cts = set(), set()
    for d in glob.glob(os.path.join(results, 'conditions',
                                    '%s_preFilt*_callThresh*' % tag)):
        base = os.path.basename(d)
        try:
            rest = base.split('_preFilt', 1)[1]
            pf = int(rest.split('_callThresh')[0])
            ct = int(rest.split('_callThresh')[1])
        except (IndexError, ValueError):
            continue
        bp = bed_path(results, expname, tag, pf, ct)
        if os.path.isfile(bp) and os.path.getsize(bp) > 0:
            pfs.add(pf); cts.add(ct)
    return sorted(pfs), sorted(cts)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-dir', required=True)
    ap.add_argument('--expname', required=True)
    ap.add_argument('--script', default='refac',
                    help='which script\'s outputs to analyze (default refac; '
                         'orig is identical)')
    ap.add_argument('--max-shift', type=int, default=1000,
                    help='max bp to still call two peaks "the same, shifted" '
                         '(default 1000 = the annotation search window)')
    ap.add_argument('--ref-prefilt', type=int, default=None,
                    help='preFilt level treated as reference (default: the '
                         'lowest present, i.e. no pre-filter if 1 is included)')
    ap.add_argument('--out-prefix', default='overlap')
    ap.add_argument('--dump-dir', default=None,
                    help='if set, write lost/gained/shifted peaks as BEDs (for '
                         'the preFilt comparisons) into this directory')
    args = ap.parse_args()

    results, exp, tag = args.results_dir, args.expname, args.script
    pfs, cts = discover(results, exp, tag)
    if not pfs:
        sys.exit('No conditions found under %s/conditions for script "%s"'
                 % (results, tag))
    ref_pf = args.ref_prefilt if args.ref_prefilt is not None else pfs[0]
    if ref_pf not in pfs:
        sys.exit('reference preFilt %d not among discovered levels %s'
                 % (ref_pf, pfs))
    print('script=%s  preFilt=%s  callThresh=%s  ref_preFilt=%d  max_shift=%d bp\n'
          % (tag, pfs, cts, ref_pf, args.max_shift))

    cache = {}

    def load(pf, ct, which=None):
        which = which or tag
        key = (which, pf, ct)
        if key not in cache:
            cache[key] = read_bed(bed_path(results, exp, which, pf, ct))
        return cache[key]

    summary_rows, type_rows, bucket_rows = [], [], []

    def run(name, A, B, dump=False):
        matched, uA, uB = match_pair(A, B, args.max_shift)
        row, dists = summarize(name, A, B, matched, uA, uB)
        summary_rows.append(row)
        bucket_rows.append(bucket_row(name, dists, len(uA), len(uB)))
        type_rows.extend(type_breakdown(name, matched, uA, uB))
        if dump and args.dump_dir:
            dump_comparison(args.dump_dir, name, matched, uA, uB)
        return row

    # 1) PREFILT effect at matched callThresh (reference preFilt vs each other)
    print('=== 1. preFilt effect (reference preFilt%d) ===' % ref_pf)
    for ct in cts:
        ref = load(ref_pf, ct)
        for pf in pfs:
            if pf == ref_pf:
                continue
            name = 'preFilt:%d_vs_%d@callThresh%d' % (ref_pf, pf, ct)
            r = run(name, ref, load(pf, ct), dump=True)
            print('  callThresh%-2d  preFilt%d vs preFilt%d : exact=%d shifted=%d '
                  'lost_from_preFilt%d=%d gained=%d  Jexact=%.3f'
                  % (ct, ref_pf, pf, r['exact'], r['shifted'], ref_pf,
                     r['unique_A'], r['unique_B'], r['jaccard_exact']))

    # 2) callThresh effect at fixed preFilt
    print('\n=== 2. callThresh effect (adjacent levels, per preFilt) ===')
    for pf in pfs:
        for a, b in zip(cts, cts[1:]):
            name = 'callThresh:%d_vs_%d@preFilt%d' % (a, b, pf)
            r = run(name, load(pf, a), load(pf, b))
            print('  preFilt%-2d  callThresh%d vs callThresh%d : exact=%d shifted=%d '
                  'only_in_%d=%d only_in_%d=%d'
                  % (pf, a, b, r['exact'], r['shifted'], a, r['unique_A'],
                     b, r['unique_B']))

    # 3) ORIG vs REFAC identity at each cell (only where both BEDs exist)
    print('\n=== 3. ORIG vs REFAC coordinate identity ===')
    checked = nonident = 0
    for pf in pfs:
        for ct in cts:
            op = bed_path(results, exp, 'orig', pf, ct)
            rp = bed_path(results, exp, 'refac', pf, ct)
            if not (os.path.isfile(op) and os.path.getsize(op) > 0 and
                    os.path.isfile(rp) and os.path.getsize(rp) > 0):
                continue
            checked += 1
            A = load(pf, ct, 'orig')
            B = load(pf, ct, 'refac')
            matched, uA, uB = match_pair(A, B, 0)
            nex = sum(1 for x in matched if x[0] == 0)
            ident = (uA == [] and uB == [] and
                     nex == total_peaks(A) == total_peaks(B))
            if not ident:
                nonident += 1
                print('  DIFF at preFilt%d callThresh%d: exact=%d uniqueA=%d '
                      'uniqueB=%d' % (pf, ct, nex, len(uA), len(uB)))
    if checked:
        print('  %d of %d comparable cells differ (0 = perfect reproduction)'
              % (nonident, checked))
    else:
        print('  (no cells with both orig and refac BEDs found; skipping)')

    # write tables
    def write_tsv(path, rows):
        if not rows:
            return
        with open(path, 'w') as f:
            f.write('\t'.join(rows[0].keys()) + '\n')
            for r in rows:
                f.write('\t'.join(str(v) for v in r.values()) + '\n')

    p1 = os.path.join(results, args.out_prefix + '_summary.tsv')
    p2 = os.path.join(results, args.out_prefix + '_by_type.tsv')
    p3 = os.path.join(results, args.out_prefix + '_shift_histogram.tsv')
    write_tsv(p1, summary_rows)
    write_tsv(p2, type_rows)
    write_tsv(p3, bucket_rows)

    # verdict on standardizing to the reference preFilt level
    print('\n=== STANDARDIZATION READ (reference preFilt%d) ===' % ref_pf)
    for ct in cts:
        rs = [r for r in summary_rows
              if r['comparison'].startswith('preFilt:')
              and r['comparison'].endswith('@callThresh%d' % ct)]
        for r in rs:
            frac_lost = r['unique_A'] / r['n_A'] if r['n_A'] else 0
            frac_gain = r['unique_B'] / r['n_B'] if r['n_B'] else 0
            print('  callThresh%-2d %s: %.1f%% of preFilt%d peaks absent after '
                  'pre-filter, %.1f%% pre-filter-only, %d shifted'
                  % (ct, r['comparison'].split(':')[1], 100 * frac_lost,
                     ref_pf, 100 * frac_gain, r['shifted']))

    print('\nwrote:\n  %s\n  %s\n  %s' % (p1, p2, p3))


if __name__ == '__main__':
    main()
