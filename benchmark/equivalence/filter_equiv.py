#!/usr/bin/env python3
"""Confirm the vectorized filter implements the intended behavior -- keep every
covered position that is NOT inside any filter window -- after the legacy
tail-drop was removed.

The original filter shared one global filter-window cursor and silently dropped
every position ordered after the LAST filter window (a terminal-region bug).
That tail-drop has been removed; this test checks the searchsorted filter now
equals a plain "keep if outside every window" reference, and reports how often
the removed tail-drop would have changed the result.
"""
import random

import numpy as np


# ---- old list-based filter, verbatim (WITH the tail-drop) -- kept only to
#      measure how often the removed bug actually changed the output ----
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


# ---- intended (fixed) behavior: keep positions outside every filter window ----
def correct_filter(entries, fw):
    def in_any(e):
        return any(e['strand'] == w['strand']
                   and e['chromosome'] == w['chromosome']
                   and w['start'] <= e['start'] <= w['end'] for w in fw)
    return [e for e in entries if not in_any(e)]


# ---- searchsorted grouped filter (mirrors TSScall filterBedGraphListByWindows,
#      tail-drop removed) ----
def searchsorted_filter(order, groups, filter_windows):
    if filter_windows == []:
        return order, groups
    windows_by_key = {}
    for w in filter_windows:
        windows_by_key.setdefault(
            (w['strand'], w['chromosome']), []).append(w)
    new_groups = {}
    new_order = []
    for key in order:
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
        idx = np.nonzero(keep)[0]
        if len(idx):
            new_groups[key] = {'pos': pos[idx], 'reads': reads[idx]}
            new_order.append(key)
    return new_order, new_groups


def make_case(rng):
    chroms = ['chr1', 'chr2', 'chr10']       # note string order chr1<chr10<chr2
    strands = ['+', '-']
    entries = []
    for s in strands:
        for c in chroms:
            k = rng.randint(0, 12)
            ps = rng.sample(range(0, 400), k) if k else []
            for p in ps:
                entries.append({'strand': s, 'chromosome': c, 'start': p,
                                'end': p, 'reads': rng.randint(1, 9)})
    entries.sort(key=lambda e: (e['strand'], e['chromosome'], e['start']))
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
    from collections import OrderedDict
    tmp = OrderedDict()
    for e in entries:
        tmp.setdefault((e['strand'], e['chromosome']), []).append(e)
    order = []
    groups = {}
    for key in sorted(tmp.keys()):
        arr = tmp[key]
        groups[key] = {'pos': np.array([a['start'] for a in arr], dtype=np.int32),
                       'reads': np.array([a['reads'] for a in arr], dtype=np.int32)}
        order.append(key)
    return order, groups


def keys_of(order, groups):
    out = []
    for key in order:
        for p in groups[key]['pos']:
            out.append((key[0], key[1], int(p)))
    return out


def main():
    rng = random.Random(99)
    trials = 40000
    mism = 0
    bug_footprint = 0
    for _ in range(trials):
        entries, fw = make_case(rng)
        correct = [(e['strand'], e['chromosome'], e['start'])
                   for e in correct_filter(entries, fw)]
        order, groups = to_groups(entries)
        got = keys_of(*searchsorted_filter(order, groups, fw))
        if got != correct:
            mism += 1
            if mism <= 3:
                print("MISMATCH\n correct:", correct, "\n got    :", got)
        old = [(e['strand'], e['chromosome'], e['start'])
               for e in original_filter(list(entries), fw)]
        if old != correct:
            bug_footprint += 1
    print("trials={}  mismatches(fixed filter vs intended)={}".format(
        trials, mism))
    print("cases where the removed tail-drop had changed the result: {}".format(
        bug_footprint))
    print("PASS" if mism == 0 else "FAIL")


if __name__ == '__main__':
    main()
