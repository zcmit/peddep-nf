#!/usr/bin/env python3
"""Keep only complete read pairs from a name-sorted (or collated) BAM.

After region/MAPQ filtering one mate of a pair can be removed while the other
stays; such orphans would be miscounted as fragments, so they are dropped here.
"""
import sys

import pysam


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: rm_orphan_pairs.py name_sorted_in.bam out.bam")
    kept = dropped = 0
    with pysam.AlignmentFile(sys.argv[1], "rb") as bam_in, \
            pysam.AlignmentFile(sys.argv[2], "wb", template=bam_in) as bam_out:
        group = []
        name = None

        def flush(reads):
            r1 = [r for r in reads if r.is_read1]
            r2 = [r for r in reads if r.is_read2]
            if len(r1) == 1 and len(r2) == 1:
                bam_out.write(r1[0])
                bam_out.write(r2[0])
                return 1, len(reads) - 2
            return 0, len(reads)

        for read in bam_in.fetch(until_eof=True):
            if read.query_name != name:
                if group:
                    k, d = flush(group)
                    kept += k
                    dropped += d
                group = []
                name = read.query_name
            group.append(read)
        if group:
            k, d = flush(group)
            kept += k
            dropped += d
    print("rm_orphan_pairs: kept %d pairs, dropped %d reads" % (kept, dropped), file=sys.stderr)


if __name__ == "__main__":
    main()
