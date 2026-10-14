#!/usr/bin/env python3
"""Write transcription start sites from a GTF as a sorted 6-column BED.

Each transcript's span is taken from its `transcript` row or, when the GTF has
none (e.g. Ensembl yeast), from the union of its exons. Protein-coding
transcripts are kept when biotype attributes exist; otherwise all transcripts
are used. Identical (contig, position, strand) sites are written once.
"""
import argparse
import gzip
import re
import sys

ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')
BIOTYPE_KEYS = ("transcript_type", "transcript_biotype", "gene_type", "gene_biotype")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    transcripts = {}
    with open_text(args.gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon"):
                continue
            attrs = dict(ATTRIBUTE.findall(fields[8]))
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            start, end = int(fields[3]), int(fields[4])
            record = transcripts.get(tid)
            if record is None:
                biotypes = [attrs[k] for k in BIOTYPE_KEYS if k in attrs]
                transcripts[tid] = {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "strand": fields[6],
                    "name": attrs.get("gene_name") or attrs.get("gene_id") or tid,
                    "has_biotype": bool(biotypes),
                    "coding": "protein_coding" in biotypes,
                }
            else:
                record["start"] = min(record["start"], start)
                record["end"] = max(record["end"], end)

    if not transcripts:
        sys.exit("ERROR: no transcript or exon records with transcript_id found in %s" % args.gtf)

    records = list(transcripts.values())
    coding = [r for r in records if r["coding"]]
    if coding:
        records = coding

    sites = {}
    for r in records:
        if r["strand"] == "-":
            pos = r["end"] - 1
        else:
            pos = r["start"] - 1
        key = (r["chrom"], pos, r["strand"] if r["strand"] in "+-" else "+")
        sites.setdefault(key, r["name"])

    with open(args.out, "w") as out:
        for (chrom, pos, strand), name in sorted(sites.items()):
            out.write("%s\t%d\t%d\t%s\t0\t%s\n" % (chrom, pos, pos + 1, name, strand))

    kind = "protein-coding" if coding else "all"
    print("Wrote %d unique TSS (%s transcripts) to %s" % (len(sites), kind, args.out), file=sys.stderr)


if __name__ == "__main__":
    main()
