#!/usr/bin/env python3
"""Prove the two-pointer bin_winner == original, and show tiling != original."""
import random
from operator import itemgetter


def original_bin_winner(hits, bin_size):
    """Verbatim transcription of TSScall.py callTSS bin_winner branch."""
    hits = sorted(hits, key=itemgetter(0))
    bins = []
    for i in range(len(hits)):
        bins.append({'total_reads': 0, 'bin_hits': []})
        for j in range(i, len(hits)):
            if abs(hits[i][0] - hits[j][0]) <= bin_size:
                bins[-1]['total_reads'] += hits[j][1]
                bins[-1]['bin_hits'].append(hits[j])
    max_bin_reads = float('-inf')
    max_bin_index = None
    for i, entry in enumerate(bins):
        if entry['total_reads'] > max_bin_reads:
            max_bin_index = i
            max_bin_reads = entry['total_reads']
    max_reads = float('-inf')
    max_position = None
    for hit in bins[max_bin_index]['bin_hits']:
        if hit[1] > max_reads:
            max_position = hit[0]
            max_reads = hit[1]
    return max_position, max_reads


def fast_bin_winner(hits, bin_size):
    """O(H) two-pointer + prefix sum. Must match original exactly."""
    hits = sorted(hits, key=itemgetter(0))
    n = len(hits)
    pos = [h[0] for h in hits]
    reads = [h[1] for h in hits]
    prefix = [0] * (n + 1)
    for k in range(n):
        prefix[k + 1] = prefix[k] + reads[k]
    best_total = float('-inf')
    best_i = None
    r = 0
    r_of = [0] * n
    for i in range(n):
        if r < i:
            r = i
        while r + 1 < n and pos[r + 1] - pos[i] <= bin_size:
            r += 1
        total = prefix[r + 1] - prefix[i]
        r_of[i] = r
        if total > best_total:      # strict > -> first (lowest-position) bin
            best_total = total
            best_i = i
    max_reads = float('-inf')
    max_position = None
    for k in range(best_i, r_of[best_i] + 1):
        if reads[k] > max_reads:    # strict > -> first (lowest-position) winner
            max_reads = reads[k]
            max_position = pos[k]
    return max_position, max_reads


def tiling_bin_winner(hits, bin_size):
    """The TRAP: fixed non-overlapping tiles (pos // bin_size)."""
    from collections import defaultdict
    tiles = defaultdict(list)
    for h in sorted(hits, key=itemgetter(0)):
        tiles[h[0] // bin_size].append(h)
    best_total = float('-inf')
    best_key = None
    for key in sorted(tiles):
        total = sum(h[1] for h in tiles[key])
        if total > best_total:
            best_total = total
            best_key = key
    max_reads = float('-inf')
    max_position = None
    for h in tiles[best_key]:
        if h[1] > max_reads:
            max_reads = h[1]
            max_position = h[0]
    return max_position, max_reads


def make_hits(rng, mode):
    n = rng.randint(1, 40)
    if mode == 'sparse':
        positions = rng.sample(range(0, 20000), n)
    elif mode == 'dense':                       # tight cluster, straddles tiles
        base = rng.randint(0, 5000)
        positions = [base + rng.randint(0, 600) for _ in range(n)]
        positions = list(dict.fromkeys(positions))  # unique
    else:                                       # mixed
        positions = [rng.randint(0, 3000) for _ in range(n)]
        positions = list(dict.fromkeys(positions))
    hits = [[p, rng.randint(1, 12)] for p in positions]   # ties in reads common
    return hits


def main():
    rng = random.Random(1234)
    bin_size = 200
    trials = 60000
    fast_mismatch = 0
    tiling_mismatch = 0
    examples = []
    for t in range(trials):
        mode = rng.choice(['sparse', 'dense', 'dense', 'mixed'])
        hits = make_hits(rng, mode)
        orig = original_bin_winner(hits, bin_size)
        fast = fast_bin_winner(hits, bin_size)
        tile = tiling_bin_winner(hits, bin_size)
        if orig != fast:
            fast_mismatch += 1
            if len(examples) < 5:
                examples.append(('FAST', mode, hits, orig, fast))
        if orig != tile:
            tiling_mismatch += 1
            if mode == 'dense' and len([e for e in examples if e[0] == 'TILE']) < 3:
                examples.append(('TILE', mode, hits, orig, tile))
    print(f"trials: {trials}, bin_size={bin_size}")
    print(f"two-pointer (proposed) mismatches vs original: {fast_mismatch}")
    print(f"tiling (prior-bug model) mismatches vs original: {tiling_mismatch}")
    print()
    for tag, mode, hits, orig, other in examples:
        if tag == 'TILE':
            print(f"[tiling divergence, {mode}] original={orig} tiling={other}")
            print(f"   hits(sorted)={sorted(hits, key=itemgetter(0))}")
    print()
    print("PASS" if fast_mismatch == 0 else "FAIL: two-pointer diverges")


if __name__ == '__main__':
    main()
