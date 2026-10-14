#!/usr/bin/env python3
"""Print the number of non-N bases in a FASTA (used as the MACS3 genome size)."""
import gzip
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: genome_size.py genome.fa[.gz]")
    path = sys.argv[1]
    opener = gzip.open if path.endswith(".gz") else open
    total = 0
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq = line.strip().upper()
            total += len(seq) - seq.count("N")
    print(total)


if __name__ == "__main__":
    main()
