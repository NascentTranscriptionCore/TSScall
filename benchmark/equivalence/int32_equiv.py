#!/usr/bin/env python3
"""Verify the int32 coverage arrays preserve behavior, and -- the one hazard md5
end-to-end tests cannot reach on normal-scale data -- that the genome-wide
read-count accumulation does not overflow int32.

The stored coverage arrays (pos, reads) are int32: positions (<= ~250 M) and
per-position read counts both fit. Output is unchanged because every value is
cast to Python int on the way into hits/output. The FDR path's genome-wide
read_count sum must accumulate in int64, because a deep library exceeds 2^31
total reads; BedGraph.total_reads uses .sum(dtype=np.int64) for exactly this.
"""
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import TSScall  # noqa: E402  (path set above; __main__ guard makes import safe)


def make_bg(groups):
    bg = TSScall.BedGraph()
    for key, (pos, reads) in groups.items():
        bg.add_group(key, np.asarray(pos, dtype=np.int32),
                     np.asarray(reads, dtype=np.int32))
    bg.finalize_order()
    return bg


def ref_total(groups):
    return sum(int(r) for _, (_, reads) in groups.items() for r in reads)


def test_overflow():
    # 200 positions * 20,000,000 reads = 4.0e9, which exceeds int32 max
    # (2,147,483,647). total_reads must return the true total, not a wrapped
    # int32 value. (On 64-bit platforms numpy already widens an int32 .sum()
    # accumulator to int64; the explicit dtype=int64 makes that guarantee
    # platform-independent -- on a 32-bit platform a bare .sum() would wrap.)
    reads = [20_000_000] * 200
    pos = list(range(1, 201))
    bg = make_bg({('+', 'chr1'): (pos, reads)})
    got = bg.total_reads()
    assert got == 4_000_000_000, got
    assert got > 2 ** 31
    assert bg.count_loci(1) == 200
    # a 32-bit accumulation of the same int32 values wraps (intentional
    # overflow here), which is exactly what the int64 accumulation avoids:
    with np.errstate(over='ignore'):
        wrapped = np.int32(0)
        for r in np.asarray(reads, dtype=np.int32):
            wrapped = np.int32(wrapped + r)
    assert int(wrapped) != got
    print("overflow: total_reads={} (> 2^31); int32-accumulated value wraps to "
          "{}".format(got, int(wrapped)))


def main():
    test_overflow()
    rng = random.Random(7)
    trials = 20000
    mism = 0
    for _ in range(trials):
        groups = {}
        for s in ('+', '-'):
            for c in rng.sample(['chr1', 'chr2', 'chrX'], rng.randint(0, 3)):
                n = rng.randint(0, 50)
                pos = sorted(rng.sample(range(1, 250_000_000), n)) if n else []
                reads = [rng.randint(1, 10 ** rng.randint(0, 9))
                         for _ in range(n)]       # values fit int32 (< 2^31)
                if pos:
                    groups[(s, c)] = (pos, reads)
        bg = make_bg(groups)
        if bg.total_reads() != ref_total(groups):
            mism += 1
            continue
        v = rng.randint(1, 10 ** 8)
        exp_loci = sum(1 for _, (_, reads) in groups.items()
                       for r in reads if r >= v)
        if bg.count_loci(v) != exp_loci:
            mism += 1
            continue
        # positions/reads round-trip through int() unchanged (int32 -> Python int)
        ok = True
        for key, (pos, reads) in groups.items():
            g = bg.groups[key]
            if ([int(x) for x in g['pos']] != list(pos)
                    or [int(x) for x in g['reads']] != list(reads)):
                ok = False
                break
        if not ok:
            mism += 1
    print("trials={}  mismatches={}".format(trials, mism))
    print("PASS" if mism == 0 else "FAIL")


if __name__ == '__main__':
    main()
