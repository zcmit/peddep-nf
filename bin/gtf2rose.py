#!/usr/bin/env python3
"""Convert a GTF to the UCSC refGene-style table ROSE uses for TSS exclusion and gene mapping.

Port of the conversion embedded in the PEDDEP/SEAseq ROSE task: one row per
`transcript` record (or `gene` record when the GTF has no transcripts), with
transcript_id as name and gene_name as name2. Unlike the WDL version, contig
names are kept as they are (no forced "chr" prefix) so non-UCSC genomes work.
"""
import argparse
import gzip
import re
import sys

ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    features = set()
    with open_text(args.gtf) as handle:
        for line in handle:
            if not line.startswith("#"):
                fields = line.split("\t")
                if len(fields) > 2:
                    features.add(fields[2])
    if "transcript" in features:
        feature = "transcript"
    elif "gene" in features:
        feature = "gene"
    elif "exon" in features:
        feature = "exon"
    else:
        sys.exit("ERROR: the GTF has no transcript, gene or exon records")

    spans = {}
    order = []
    with open_text(args.gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != feature:
                continue
            attrs = dict(ATTRIBUTE.findall(fields[8]))
            name = attrs.get("transcript_id") or attrs.get("gene_id") or fields[8].split(";")[0]
            gene = attrs.get("gene_name") or attrs.get("gene_id") or name
            start, end = int(fields[3]), int(fields[4])
            if feature == "exon":
                # GTFs without transcript rows (e.g. Ensembl yeast): span = union of exons
                if name in spans:
                    spans[name][3] = min(spans[name][3], start)
                    spans[name][4] = max(spans[name][4], end)
                    continue
            elif name in spans:
                name = "%s_%d" % (name, len(order))
            spans[name] = [name, fields[0], fields[6], start, end, gene]
            order.append(name)

    with open(args.out, "w") as out:
        out.write("#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\tX\tX\tX\t\tX\tname2\n")
        for key in order:
            name, chrom, strand, start, end, gene = spans[key]
            out.write("0\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t.\t.\t.\t.\t%s\n" % (name, chrom, strand, start, end, start, end, gene))
    print("Wrote %d %s records to %s" % (len(order), feature, args.out), file=sys.stderr)


if __name__ == "__main__":
    main()
