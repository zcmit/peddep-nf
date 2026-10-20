#!/usr/bin/env python3
"""Fast drop-in replacement for ROSE_bamToGFF.py (ROSE 1.3.x).

ROSE_main.py maps every BAM to every region with ROSE_bamToGFF.py, which runs one
`samtools view` per region and counts coverage base by base in Python dicts. On a
deep H3K27ac library that takes hours. This script gives the same table from one
streaming pass over the BAM, with the same arithmetic:

  * each read covers POS .. POS+len(SEQ), extended by -e bp in its own direction
    (reads with an N in the CIGAR are skipped, as in ROSE);
  * positions covered by more than -f reads count;
  * a region's value is the summed coverage strictly inside it divided by
    (length - 1) / bins, per million mapped reads with -r (samtools flagstat).

The ROSE module puts this first on PATH as ROSE_bamToGFF.py; options and output
format are those of ROSE_bamToGFF.py. Needs python3, numpy and samtools.
"""
import argparse
import math
import os
import re
import subprocess
import sys
from array import array

import numpy as np


def mapped_millions(bam):
    out = subprocess.run(['samtools', 'flagstat', bam], stdout=subprocess.PIPE, check=True).stdout.decode()
    for line in out.splitlines():
        if line.count('mapped (') == 1:
            return round(float(int(line.split(' ')[0])) / 1000000, 4)
    sys.exit('ERROR: no mapped read count in samtools flagstat %s' % bam)


def read_regions(gff, strip_chr):
    """GFF lines -> [chrom, start, end, sense, ID] like ROSE_utils.Locus."""
    regions = []
    with open(gff) as handle:
        for raw in handle:
            line = raw.rstrip('\n').split('\t')
            if len(line) < 7:
                continue
            start, end = sorted((int(line[3]), int(line[4])))
            chrom = re.sub('chr', '', line[0]) if strip_chr else line[0]
            regions.append([chrom, start, end, line[6], line[1]])
    return regions


class Coverage(object):
    """Piecewise-constant coverage of one chromosome; F(p) = sum of counted coverage over x < p."""

    def __init__(self, starts, ends, floor):
        first = np.frombuffer(starts, dtype=np.int64)
        after = np.frombuffer(ends, dtype=np.int64) + 1
        events = np.concatenate((first, after))
        delta = np.concatenate((np.ones(len(first), np.int64), -np.ones(len(after), np.int64)))
        self.bounds, inverse = np.unique(events, return_inverse=True)
        cov = np.cumsum(np.bincount(inverse, weights=delta).astype(np.int64))
        self.weight = np.where(cov > floor, cov, 0)
        self.integral = np.concatenate(([0], np.cumsum(self.weight[:-1] * np.diff(self.bounds))))

    def upto(self, pos):
        k = int(np.searchsorted(self.bounds, pos, side='right')) - 1
        if k < 0:
            return 0
        return int(self.integral[k]) + int(self.weight[k]) * (pos - int(self.bounds[k]))

    def between(self, low, high):
        """Counted coverage summed over integers x with low < x < high (float bounds allowed)."""
        first, last = math.floor(low) + 1, math.ceil(high)
        return self.upto(last) - self.upto(first) if last > first else 0


class NoCoverage(object):
    def between(self, low, high):
        return 0


def chromosome_coverage(bam, extension, floor):
    """Yield (chrom, Coverage) for each chromosome of a coordinate-sorted BAM."""
    view = subprocess.Popen(['samtools', 'view', bam], stdout=subprocess.PIPE)
    current, starts, ends, done = None, array('q'), array('q'), set()
    for line in view.stdout:
        fields = line.split(b'\t', 10)
        chrom = fields[2]
        if chrom == b'*':
            continue
        if b'N' in fields[5]:
            continue
        if chrom != current:
            if current is not None:
                yield current.decode(), Coverage(starts, ends, floor)
                done.add(current)
            if chrom in done:
                sys.exit('ERROR: %s is not coordinate-sorted' % bam)
            current, starts, ends = chrom, array('q'), array('q')
        pos = int(fields[3])
        length = len(fields[9])
        if int(fields[1]) & 16:
            starts.append(pos - extension)
            ends.append(pos + length)
        else:
            starts.append(pos)
            ends.append(pos + length + extension)
    if view.wait() != 0:
        sys.exit('ERROR: samtools view %s failed' % bam)
    if current is not None:
        yield current.decode(), Coverage(starts, ends, floor)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--bam', required=True)
    parser.add_argument('-i', '--input', required=True, help='GFF of regions')
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-s', '--sense', default='both', choices=['both', '.'])
    parser.add_argument('-f', '--floor', type=int, default=0)
    parser.add_argument('-e', '--extension', type=int, default=200)
    parser.add_argument('-r', '--rpm', action='store_true')
    parser.add_argument('-m', '--matrix', type=int, required=True, help='bins per region')
    args = parser.parse_args()

    first = subprocess.run('samtools view %s | head -n 1' % args.bam, shell=True,
                           stdout=subprocess.PIPE).stdout.decode().split('\t')
    has_chr = len(first) > 2 and 'chr' in first[2]
    mmr = mapped_millions(args.bam) if args.rpm else 1
    regions = read_regions(args.input, strip_chr=not has_chr)

    by_chrom = {}
    for index, region in enumerate(regions):
        by_chrom.setdefault(region[0], []).append(index)
    values = [None] * len(regions)

    def fill(chrom, coverage):
        for index in by_chrom.pop(chrom, []):
            _, start, end, sense, _ = regions[index]
            bin_size = (end - start) / args.matrix
            if bin_size == 0:
                values[index] = ['NA'] * args.matrix
                continue
            row = []
            if sense == '-':
                pos = end
                for _ in range(args.matrix):
                    row.append(round(float(coverage.between(pos - bin_size, pos)) / bin_size / mmr, 4))
                    pos = pos - bin_size
            else:
                pos = start
                for _ in range(args.matrix):
                    row.append(round(float(coverage.between(pos, pos + bin_size)) / bin_size / mmr, 4))
                    pos = pos + bin_size
            values[index] = row

    for chrom, coverage in chromosome_coverage(args.bam, args.extension, args.floor):
        fill(chrom, coverage)
    for chrom in list(by_chrom):
        fill(chrom, NoCoverage())

    bam_name = args.bam.split('/')[-1]
    partial = args.output + '.partial'
    with open(partial, 'w') as out:
        out.write('\t'.join(['GENE_ID', 'locusLine'] + ['bin_%d_%s' % (n, bam_name) for n in range(1, args.matrix + 1)]) + '\n')
        for (chrom, start, end, sense, ident), row in zip(regions, values):
            locus = '%s%s(%s):%d-%d' % ('' if has_chr else 'chr', chrom, sense, start, end)
            out.write('\t'.join([ident, locus] + [str(v) for v in row]) + '\n')
    # ROSE_main.py polls for the output file, so it must appear only when complete
    os.rename(partial, args.output)


if __name__ == '__main__':
    main()
