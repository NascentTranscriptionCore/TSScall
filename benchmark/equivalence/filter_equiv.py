#!/usr/bin/env python3
"""Confirm the grouped filter exactly reproduces the original list-based filter,
including the legacy tail-drop and cross-(strand,chrom) cursor behavior."""
import random
import numpy as np


# ---- original logic, verbatim ----
def isWithin(e1, e2):
    if e1['strand'] == e2['strand'] and e1['chromosome'] == e2['chromosome']:
        if e1['start'] >= e2['start'] and e1['end'] <= e2['end']:
            return True
    return False


def isLessThan(e1, e2):
    for val in ['strand', 'chromosome', 'start']:
        if e1[val] < e2[val]:
            return True
        elif e1[val] > e2[val]:
            return False
    return False


def original_filter(bedgraph_list, filter_windows):
    if filter_windows != []:
        fi = 0
        bi = 0
        working = []
        while (fi < len(filter_windows)) and (bi < len(bedgraph_list)):
            if isWithin(bedgraph_list[bi], filter_windows[fi]):
                bi += 1
            else:
                if isLessThan(bedgraph_list[bi], filter_windows[fi]):
                    working.append(bedgraph_list[bi])
                    bi += 1
                else:
                    fi += 1
        bedgraph_list = working
    return bedgraph_list


# ---- refactored grouped logic (mirrors TSScall_refactored) ----
def _group_cmp(a_s, a_c, b_s, b_c):
    if a_s != b_s:
        return -1 if a_s < b_s else 1
    if a_c != b_c:
        return -1 if a_c < b_c else 1
    return 0


def grouped_filter(order, groups, filter_windows):
    if filter_windows == []:
        return order, groups
    FW = filter_windows
    nfw = len(FW)
    fi = 0
    new_groups = {}
    new_order = []
    for key in order:
        if fi >= nfw:
            break
        strand, chrom = key
        pos = groups[key]['pos']
        reads = groups[key]['reads']
        n = len(pos)
        keep = np.zeros(n, dtype=bool)
        bi = 0
        while (fi < nfw) and (bi < n):
            w = FW[fi]
            gc = _group_cmp(strand, chrom, w['strand'], w['chromosome'])
            if gc == 0:
                p = pos[bi]
                if w['start'] <= p <= w['end']:
                    bi += 1
                elif p < w['start']:
                    keep[bi] = True
                    bi += 1
                else:
                    fi += 1
            elif gc < 0:
                keep[bi] = True
                bi += 1
            else:
                fi += 1
        idx = np.nonzero(keep)[0]
        if len(idx):
            new_groups[key] = {'pos': pos[idx], 'reads': reads[idx]}
            new_order.append(key)
    return new_order, new_groups


# ---- searchsorted grouped logic (mirrors TSScall filterBedGraphListByWindows) ----
def searchsorted_filter(order, groups, filter_windows):
    if filter_windows == []:
        return order, groups
    windows_by_key = {}
    for w in filter_windows:
        windows_by_key.setdefault(
            (w['strand'], w['chromosome']), []).append(w)
    last_w = filter_windows[-1]
    last_key = (last_w['strand'], last_w['chromosome'])
    last_end = last_w['end']
    new_groups = {}
    new_order = []
    for key in order:
        if key > last_key:                     # legacy tail-drop: later groups
            break
        pos = groups[key]['pos']
        reads = groups[key]['reads']
        keep = np.ones(len(pos), dtype=bool)
        wins = windows_by_key.get(key)
        if wins:
            starts = np.fromiter((w['start'] for w in wins), dtype=np.int32,
                                 count=len(wins))
            ends = np.fromiter((w['end'] for w in wins), dtype=np.int32,
                               count=len(wins))
            los = np.searchsorted(pos, starts, side='left')
            his = np.searchsorted(pos, ends, side='right')
            for lo, hi in zip(los.tolist(), his.tolist()):
                keep[lo:hi] = False
        if key == last_key:                    # legacy tail-drop within group
            cut = int(np.searchsorted(pos, np.int32(last_end), side='right'))
            keep[cut:] = False
        idx = np.nonzero(keep)[0]
        if len(idx):
            new_groups[key] = {'pos': pos[idx], 'reads': reads[idx]}
            new_order.append(key)
    return new_order, new_groups


def make_case(rng):
    chroms = ['chr1', 'chr2', 'chr10']       # note string order chr1<chr10<chr2
    strands = ['+', '-']
    # bedgraph: unique (strand,chrom,pos)
    entries = []
    for s in strands:
        for c in chroms:
            k = rng.randint(0, 12)
            ps = rng.sample(range(0, 400), k) if k else []
            for p in ps:
                entries.append({'strand': s, 'chromosome': c, 'start': p,
                                'end': p, 'reads': rng.randint(1, 9)})
    entries.sort(key=lambda e: (e['strand'], e['chromosome'], e['start']))
    # filter windows: random, then sorted like sortList('sort_by_strand')
    fw = []
    for _ in range(rng.randint(0, 6)):
        s = rng.choice(strands)
        c = rng.choice(chroms)
        a = rng.randint(0, 380)
        w = rng.randint(0, 40)
        fw.append({'strand': s, 'chromosome': c, 'start': a, 'end': a + w})
    fw.sort(key=lambda e: (e['strand'], e['chromosome'], e['start']))
    return entries, fw


def to_groups(entries):
    order = []
    groups = {}
    from collections import OrderedDict
    tmp = OrderedDict()
    for e in entries:
        tmp.setdefault((e['strand'], e['chromosome']), []).append(e)
    for key in sorted(tmp.keys()):
        arr = tmp[key]
        groups[key] = {'pos': np.array([a['start'] for a in arr]),
                       'reads': np.array([a['reads'] for a in arr])}
        order.append(key)
    return order, groups


def tail_dropped_count(entries, fw, kept_orig):
    # positions the original dropped that are NOT inside any window == tail-drop
    def in_any(e):
        return any(e['strand'] == w['strand']
                   and e['chromosome'] == w['chromosome']
                   and w['start'] <= e['start'] <= w['end'] for w in fw)
    not_in_window = sum(1 for e in entries if not in_any(e))
    return not_in_window - len(kept_orig)


def keys_of(order, groups):
    out = []
    for key in order:
        for p in groups[key]['pos']:
            out.append((key[0], key[1], int(p)))
    return out


def main():
    rng = random.Random(99)
    trials = 40000
    mism_grouped = 0
    mism_ss = 0
    tail_fired = 0
    for _ in range(trials):
        entries, fw = make_case(rng)
        kept_orig = original_filter(list(entries), fw)
        orig_keys = [(e['strand'], e['chromosome'], e['start'])
                     for e in kept_orig]
        order, groups = to_groups(entries)
        if keys_of(*grouped_filter(order, groups, fw)) != orig_keys:
            mism_grouped += 1
        if keys_of(*searchsorted_filter(order, groups, fw)) != orig_keys:
            mism_ss += 1
            if mism_ss <= 3:
                print("MISMATCH (searchsorted)")
        if fw and tail_dropped_count(entries, fw, kept_orig) > 0:
            tail_fired += 1
    print("trials={}  grouped_mismatches={}  searchsorted_mismatches={}"
          .format(trials, mism_grouped, mism_ss))
    print("cases where the legacy tail-drop actually fired: {}".format(
        tail_fired))
    print("PASS" if (mism_grouped == 0 and mism_ss == 0) else "FAIL")


if __name__ == '__main__':
    main()
