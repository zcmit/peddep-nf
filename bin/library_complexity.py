#!/usr/bin/env python3
"""ENCODE library-complexity metrics from a coordinate-sorted, duplicate-marked BAM.

A fragment is identified by (contig, start, end, read-1 strand) of a properly
paired read 1 with MAPQ >= --min-mapq. Duplicates are counted (they are only
marked, not removed), so:
    NRF  = distinct fragments / total fragments
    PBC1 = fragments seen exactly once / distinct fragments
    PBC2 = fragments seen exactly once / fragments seen exactly twice
"""
import argparse
import sys
from collections import Counter

import pysam


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--exclude-contig", action="append", default=[],
                        help="contig to ignore (e.g. the mitochondrial genome); repeatable")
    args = parser.parse_args()

    excluded = set(c for c in args.exclude_contig if c)
    total = distinct = one = two = 0

    def tally(counts):
        nonlocal total, distinct, one, two
        for n in counts.values():
            total += n
            distinct += 1
            if n == 1:
                one += 1
            elif n == 2:
                two += 1

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        counts = Counter()
        current = None
        for read in bam.fetch(until_eof=True):
            if (not read.is_paired or not read.is_proper_pair or not read.is_read1
                    or read.is_unmapped or read.mate_is_unmapped or read.is_secondary
                    or read.is_supplementary or read.is_qcfail
                    or read.mapping_quality < args.min_mapq or read.template_length == 0):
                continue
            contig = read.reference_name
            if contig in excluded:
                continue
            if contig != current:
                tally(counts)
                counts = Counter()
                current = contig
            start = min(read.reference_start, read.next_reference_start)
            counts[(start, start + abs(read.template_length), read.is_reverse)] += 1
        tally(counts)

    nrf = distinct / total if total else 0.0
    pbc1 = one / distinct if distinct else 0.0
    pbc2 = "%.4f" % (one / two) if two else "NA"
    with open(args.out, "w") as out:
        out.write("Sample\tTotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2\n")
        out.write("%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%s\n" % (args.sample, total, distinct, one, two, nrf, pbc1, pbc2))


if __name__ == "__main__":
    main()
