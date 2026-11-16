#!/usr/bin/env python3
"""Prove the searchsorted intersection == the original per-group merge-walk, for
the non-overlapping, start-sorted windows that both window builders guarantee.

findIntersectionWithBedGraph records, per window, [position, reads] for every
covered position in [start, end] on the same (strand, chromosome). The original
stepped a cursor through every position; the replacement takes two
np.searchsorted bounds per window. These agree ONLY because the windows reaching
the method are non-overlapping (each position belongs to at most one window). An
overlapping-windows negative control confirms the test detects the difference,
documenting why that invariant matters.
"""
import random

import numpy as np


def original_intersection(pos, reads, win_list):
    # verbatim transcription of the per-(strand,chromosome) merge-walk
    for w in win_list:
        w['hits'] = []
    n = len(pos)
    nw = len(win_list)
    bi = 0
    si = 0
    while si < nw and bi < n:
        w = win_list[si]
        p = pos[bi]
        if w['start'] <= p <= w['end']:
            w['hits'].append([int(p), int(reads[bi])])
            bi += 1
        elif p < w['start']:
            bi += 1
        else:
            si += 1
    return [list(w['hits']) for w in win_list]


def fast_intersection(pos, reads, win_list):
    posN = np.asarray(pos, dtype=np.int32)
    readsN = np.asarray(reads, dtype=np.int32)
    starts = np.fromiter((w['start'] for w in win_list), dtype=np.int32,
                         count=len(win_list))
    ends = np.fromiter((w['end'] for w in win_list), dtype=np.int32,
                       count=len(win_list))
    los = np.searchsorted(posN, starts, side='left')
    his = np.searchsorted(posN, ends, side='right')
    out = []
    for i in range(len(win_list)):
        lo = int(los[i])
        hi = int(his[i])
        if hi > lo:
            out.append(list(map(list, zip(posN[lo:hi].tolist(),
                                          readsN[lo:hi].tolist()))))
        else:
            out.append([])
    return out


def make_case(rng, overlapping=False):
    n = rng.randint(0, 60)
    pos = sorted(rng.sample(range(0, 5000), n)) if n else []
    reads = [rng.randint(1, 20) for _ in pos]     # ties common
    wins = []
    x = rng.randint(0, 100)
    for _ in range(rng.randint(0, 12)):
        start = x + rng.randint(0, 50)
        end = start + rng.randint(0, 80)          # width 0 windows included
        wins.append({'start': start, 'end': end})
        if overlapping:
            x = start + rng.randint(0, 40)        # may overlap the previous
        else:
            x = end + 1 + rng.randint(0, 40)      # strictly non-overlapping
    return pos, reads, wins


def main():
    rng = random.Random(31)
    trials = 40000
    mism = 0
    for _ in range(trials):
        pos, reads, wins = make_case(rng, overlapping=False)
        a = original_intersection(pos, reads, [dict(w) for w in wins])
        b = fast_intersection(pos, reads, [dict(w) for w in wins])
        if a != b:
            mism += 1
    ctrl_div = 0
    for _ in range(trials):
        pos, reads, wins = make_case(rng, overlapping=True)
        a = original_intersection(pos, reads, [dict(w) for w in wins])
        b = fast_intersection(pos, reads, [dict(w) for w in wins])
        if a != b:
            ctrl_div += 1
    print("non-overlapping (the real contract): trials={} mismatches={}".format(
        trials, mism))
    print("overlapping control: diverged on {} / {} (test detects the "
          "difference)".format(ctrl_div, trials))
    print("PASS" if mism == 0 else "FAIL")


if __name__ == '__main__':
    main()
