#!/usr/bin/env python3
"""PEDDEP graded QC report for one HiChIP library (target or INPUT).

Sequencing metrics as in the other PEDDEP reports (raw pairs, base quality,
sequence diversity), contact metrics from pairtools (mapped pairs, duplicates,
unique contacts, cis / long-range cis fractions of the deduplicated contacts) and,
for targets, enrichment metrics (fraction of read ends in peaks, peaks, loops).
Graded -2..+2 with cutoffs from a JSON file; the overall grade is the mean.
INPUT libraries get no peaks or loops, so those metrics are not graded for them.

The HiLOW report divided cis counts of the deduplicated contacts by all valid
pairs, duplicates included; here every contact fraction uses deduplicated
contacts on both sides.

Writes <prefix>.hichip_qc_report.txt, .html and <prefix>.hichip_qc_metrics.json.
"""
import argparse
import json

from qc_common import fastqc_summary, read_table, to_number, build_rows, write_reports

METRICS = [
    ("Raw_Read_Pairs", "Read pairs sequenced, all lanes", "FastQC (raw)"),
    ("Mean_Read_Length", "Mean raw read length (bp)", "FastQC (raw)"),
    ("Base_Quality", "Per-base sequence quality, worst of all raw FASTQs", "FastQC (raw)"),
    ("Sequence_Diversity", "Overrepresented sequences, worst of all raw FASTQs", "FastQC (raw)"),
    ("Pct_Reads_Q30", "% reads with mean base quality >= 30 (worst FASTQ)", "FastQC (raw)"),
    ("Mapped_Pairs_Pct", "% of pairs with both ends mapped uniquely (pair types UU, UR, RU)", "pairtools"),
    ("Duplication_Pct", "% of mapped pairs that are duplicates", "pairtools dedup"),
    ("Unique_Contacts", "Deduplicated contacts (mapped pairs minus duplicates)", "pairtools dedup"),
    ("Cis_Pct", "% of unique contacts within one chromosome", "pairtools dedup"),
    ("Cis_Under_1kb_Pct", "% of unique contacts that are cis and < 1 kb apart", "pairtools dedup"),
    ("Cis_Over_20kb_Pct", "% of unique contacts that are cis and >= 20 kb apart (HiLOW cis_longRange)", "pairtools dedup"),
    ("Trans_Pct", "% of unique contacts between chromosomes", "pairtools dedup"),
    ("FRiP", "Fraction of read ends in peaks", "MACS3 + samtools"),
    ("Total_Peaks", "1D peaks (loop anchors) called vs the INPUT", "MACS3"),
    ("Loops_Q0.01", "FitHiChIP loops at q < 0.01", "FitHiChIP"),
    ("Merged_Loops_Q0.01", "FitHiChIP loops at q < 0.01 after merging nearby loops", "FitHiChIP"),
]
CONTROL_UNGRADED = ("FRiP", "Total_Peaks", "Loops_Q0.01", "Merged_Loops_Q0.01")


def read_pairtools_stats(path):
    stats = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 2:
                stats[fields[0]] = to_number(fields[1])
    return stats


def pct(numerator, denominator):
    if numerator is None or not denominator:
        return None
    return 100.0 * numerator / denominator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--control", action="store_true", help="this library is an INPUT control")
    parser.add_argument("--fastqc", nargs="+", required=True)
    parser.add_argument("--dedup-stats", required=True)
    parser.add_argument("--frip")
    parser.add_argument("--loops")
    args = parser.parse_args()

    with open(args.thresholds) as handle:
        thresholds = json.load(handle)

    fq = fastqc_summary(args.fastqc)
    stats = read_pairtools_stats(args.dedup_stats)
    total = stats.get("total")
    mapped = stats.get("total_mapped")
    nodups = stats.get("total_nodups")
    cis = stats.get("cis")
    cis_1kb = stats.get("cis_1kb+")

    values = {
        "Raw_Read_Pairs": fq["raw"],
        "Mean_Read_Length": fq["mean_length"],
        "Base_Quality": fq["base_quality"],
        "Sequence_Diversity": fq["diversity"],
        "Pct_Reads_Q30": fq["q30"],
        "Mapped_Pairs_Pct": pct(mapped, total),
        "Duplication_Pct": pct(stats.get("total_dups"), mapped),
        "Unique_Contacts": nodups,
        "Cis_Pct": pct(cis, nodups),
        "Cis_Under_1kb_Pct": pct(cis - cis_1kb, nodups) if cis is not None and cis_1kb is not None else None,
        "Cis_Over_20kb_Pct": pct(stats.get("cis_20kb+"), nodups),
        "Trans_Pct": pct(stats.get("trans"), nodups),
    }

    frip = read_table(args.frip) if args.frip else {}
    values["FRiP"] = to_number(frip.get("FRiP"))
    values["Total_Peaks"] = to_number(frip.get("Peaks"))

    loops = read_table(args.loops) if args.loops else {}
    values["Loops_Q0.01"] = to_number(loops.get("Loops_Q0.01"))
    values["Merged_Loops_Q0.01"] = to_number(loops.get("Merged_Loops_Q0.01"))

    rows = build_rows(METRICS, values, thresholds, ungraded=CONTROL_UNGRADED if args.control else ())
    write_reports(args.sample, args.prefix, "hichip", "HiChIP", "HiChIP", rows)


if __name__ == "__main__":
    main()
