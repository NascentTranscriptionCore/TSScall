#!/usr/bin/env python3
"""Prove the vectorized pandas/numpy bedGraph parse == the original per-line,
per-position parse.

The original combineAndSortBedGraphs.readBedGraph skips any line containing
'track' or a blank line, then expands each line [start, end) into one position
per integer in range(start+1, end+1) carrying the line's read count, grouped by
(strand, chromosome) and sorted by position. The proposed replacement does the
same with pandas (C parser) + numpy vectorized range-expansion.

Covers: width-1 and wider intervals, multiple chromosomes, read-count ties,
'track' header lines, blank lines, shuffled (unsorted) input, and cross-strand
grouping. Must be 0 mismatches.
"""
import io
import random
from array import array

import numpy as np
import pandas as pd


# ---- original logic, verbatim transcription ----
def original_parse(fwd_lines, rev_lines):
    acc = {}

    def readBedGraph(lines, strand):
        for line in lines:
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

    readBedGraph(fwd_lines, '+')
    readBedGraph(rev_lines, '-')
    out = {}
    for key, (pa, ra) in acc.items():
        pv = np.frombuffer(pa, dtype=np.int64)
        rv = np.frombuffer(ra, dtype=np.int64)
        order = np.argsort(pv, kind='stable')
        out[key] = (np.array(pv[order]), np.array(rv[order]))
    return out


# ---- proposed vectorized logic ----
def _isint(s):
    s2 = s[1:] if s[:1] == '-' else s
    return s2.isdigit()


def vectorized_parse(fwd_text, rev_text):
    def read_one(text, strand):
        # Skip leading header/blank lines (a 'track'/'browser' line or blanks at
        # the top of a bedGraph) so the first row pandas sees is a 4-column data
        # row -- otherwise the C parser infers the column count from the header.
        # Mid-file blanks are dropped by skip_blank_lines; a stray mid-file
        # 'track' line is dropped by the numeric coercion + chrom filter below.
        skip = 0
        for line in text.splitlines(keepends=True):
            p = line.split()
            if len(p) >= 4 and _isint(p[1]) and _isint(p[2]) and _isint(p[3]):
                break
            skip += 1
        else:
            return {}   # no data rows
        df = pd.read_csv(
            io.StringIO(text), sep='\t', header=None, usecols=[0, 1, 2, 3],
            names=['chrom', 'start', 'end', 'reads'], skiprows=skip,
            skip_blank_lines=True, dtype={'chrom': 'category'},
            keep_default_na=False, na_values=[''], comment=None,
        )
        for col in ('start', 'end', 'reads'):
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df = df.dropna(subset=['start', 'end', 'reads'])
        df = df[~df['chrom'].str.contains('track', regex=False, na=False)]
        groups = {}
        for chrom, sub in df.groupby('chrom', sort=False, observed=True):
            start = sub['start'].to_numpy(dtype=np.int64)
            end = sub['end'].to_numpy(dtype=np.int64)
            reads = sub['reads'].to_numpy(dtype=np.int64)
            widths = np.maximum(end - start, 0)
            total = int(widths.sum())
            reads_exp = np.repeat(reads, widths)
            starts_exp = np.repeat(start + 1, widths)
            offset = np.arange(total) - np.repeat(np.cumsum(widths) - widths,
                                                  widths)
            pos = starts_exp + offset
            groups[(strand, chrom)] = (pos, reads_exp)
        return groups

    out = {}
    for text, strand in ((fwd_text, '+'), (rev_text, '-')):
        for key, (pos, reads) in read_one(text, strand).items():
            order = np.argsort(pos, kind='stable')
            out[key] = (pos[order], reads[order])
    return out


def make_file(rng):
    chroms = rng.sample(['chr1', 'chr2', 'chr10', 'chrX', 'chrMT'],
                        k=rng.randint(1, 4))
    lines = []
    if rng.random() < 0.5:
        lines.append('track type=bedGraph name="x"\n')
    for c in chroms:
        p = rng.randint(1, 50)
        for _ in range(rng.randint(0, 40)):
            w = rng.randint(1, 5)                 # width-1 and wider
            rd = rng.randint(1, 8)                # ties common
            lines.append(f"{c}\t{p}\t{p + w}\t{rd}\n")
            p += w + rng.randint(0, 3)            # non-overlapping, with gaps
        if rng.random() < 0.3:
            lines.append('\n')                    # blank line
    # Shuffle sometimes so the input is not position-sorted (stresses the sort);
    # 'track'/blank filtering is position-independent in both implementations.
    body = lines[1:] if lines and lines[0].startswith('track') else lines
    head = lines[:1] if lines and lines[0].startswith('track') else []
    if rng.random() < 0.4:
        rng.shuffle(body)
    return "".join(head + body)


def eq(a, b):
    if set(a.keys()) != set(b.keys()):
        return False
    for k in a:
        pa, ra = a[k]
        pb, rb = b[k]
        if not (np.array_equal(pa, pb) and np.array_equal(ra, rb)):
            return False
    return True


def main():
    rng = random.Random(2024)
    trials = 20000
    mism = 0
    with_expansion = 0
    for t in range(trials):
        fwd = make_file(rng)
        rev = make_file(rng)
        o = original_parse(fwd.splitlines(keepends=True),
                           rev.splitlines(keepends=True))
        v = vectorized_parse(fwd, rev)
        if any('\t' in ln and int(ln.split('\t')[2]) - int(ln.split('\t')[1]) > 1
               for ln in (fwd + rev).splitlines() if ln and 'track' not in ln):
            with_expansion += 1
        if not eq(o, v):
            mism += 1
            if mism <= 3:
                print("MISMATCH on trial", t)
                print(" fwd:", repr(fwd[:200]))
    print(f"trials={trials}  mismatches={mism}  "
          f"trials_with_wide_intervals={with_expansion}")
    print("PASS" if mism == 0 else "FAIL")


if __name__ == '__main__':
    main()
