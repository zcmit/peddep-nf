#!/usr/bin/env python3
"""PEDDEP graded QC report for one ATAC-seq sample.

Same grading idea as the PEDDEP/SEAseq evaluation-statistics.pl report: every
graded metric gets a score from -2 to 2 and the overall grade is the mean
score (POOR / BELOW_AVERAGE / AVERAGE / GOOD / EXCELLENT). Metrics come first
from sequencing quality (FastQC on the raw reads), then alignment and library
complexity, then ATAC enrichment. Thresholds are read from a JSON file.

Writes <prefix>.atac_qc_report.txt, .html and <prefix>.atac_qc_metrics.json.
"""
import argparse
import json

from qc_common import (fastqc_summary, read_trim_reports, read_bowtie2_logs, read_idxstats,
                       read_picard_metrics, read_table, read_flagstat_read1, read_ataqv,
                       to_number, build_rows, write_reports)

# Display order and descriptions; a metric is graded when it has a thresholds entry
METRICS = [
    ("Raw_Read_Pairs", "Read pairs sequenced (all lanes)", "FastQC (raw)"),
    ("Mean_Read_Length", "Mean raw read length (bp)", "FastQC (raw)"),
    ("Base_Quality", "Per-base sequence quality, worst of all raw FASTQs", "FastQC (raw)"),
    ("Sequence_Diversity", "Overrepresented sequences, worst of all raw FASTQs", "FastQC (raw)"),
    ("Pct_Reads_Q30", "% reads with mean base quality >= 30 (worst mate)", "FastQC (raw)"),
    ("Adapter_Trimmed_Pct", "% reads with adapter trimmed", "Trim Galore"),
    ("Aligned_Percent", "% reads aligned (concordant pairs)", "bowtie2"),
    ("Mito_Percent", "% of mapped reads on the mitochondrial genome", "samtools idxstats"),
    ("Duplication_Pct", "% duplicate reads (marked, then removed)", "Picard MarkDuplicates"),
    ("NRF", "Non-redundant fraction (distinct / total fragments)", "library complexity"),
    ("PBC1", "PCR bottleneck coefficient 1 (once / distinct)", "library complexity"),
    ("PBC2", "PCR bottleneck coefficient 2 (once / twice)", "library complexity"),
    ("Final_Read_Pairs", "Read pairs kept after all filters", "samtools flagstat"),
    ("Median_Fragment", "Median fragment length (bp)", "fragment QC"),
    ("NFR_Fraction", "Fraction of fragments < 100 bp", "fragment QC"),
    ("Nucleosome_Pattern", "Nucleosome-free + mono-nucleosome peaks visible", "fragment QC"),
    ("TSS_Enrichment", "TSS enrichment score", "ataqv"),
    ("FRiP", "Fraction of reads in narrow peaks", "MACS3 + samtools"),
    ("Total_Peaks", "Narrow peaks called (unique intervals)", "MACS3"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--mito-name", default="chrM")
    parser.add_argument("--fastqc", nargs="+", required=True, help="raw-read FastQC zips (R1 and R2, all lanes)")
    parser.add_argument("--trim-logs", nargs="*", default=[])
    parser.add_argument("--bowtie2-logs", nargs="*", default=[])
    parser.add_argument("--idxstats")
    parser.add_argument("--markdup-metrics")
    parser.add_argument("--complexity")
    parser.add_argument("--filtered-flagstat")
    parser.add_argument("--frip")
    parser.add_argument("--ataqv")
    parser.add_argument("--fragment-qc")
    args = parser.parse_args()

    with open(args.thresholds) as handle:
        thresholds = json.load(handle)

    fq = fastqc_summary(args.fastqc)
    values = {
        "Raw_Read_Pairs": fq["raw"],
        "Mean_Read_Length": fq["mean_length"],
        "Base_Quality": fq["base_quality"],
        "Sequence_Diversity": fq["diversity"],
        "Pct_Reads_Q30": fq["q30"],
        "Adapter_Trimmed_Pct": read_trim_reports(args.trim_logs) if args.trim_logs else None,
        "Aligned_Percent": read_bowtie2_logs(args.bowtie2_logs) if args.bowtie2_logs else None,
        "Mito_Percent": read_idxstats(args.idxstats, set(args.mito_name.split(","))) if args.idxstats else None,
        "Duplication_Pct": read_picard_metrics(args.markdup_metrics) if args.markdup_metrics else None,
        "Final_Read_Pairs": read_flagstat_read1(args.filtered_flagstat) if args.filtered_flagstat else None,
        "TSS_Enrichment": read_ataqv(args.ataqv) if args.ataqv else None,
    }

    complexity = read_table(args.complexity) if args.complexity else {}
    for metric in ("NRF", "PBC1", "PBC2"):
        values[metric] = to_number(complexity.get(metric))

    fragments = read_table(args.fragment_qc) if args.fragment_qc else {}
    values["Median_Fragment"] = to_number(fragments.get("Median_Fragment"))
    values["NFR_Fraction"] = to_number(fragments.get("NFR_Fraction"))
    values["Nucleosome_Pattern"] = fragments.get("Nucleosome_Pattern") or None

    frip = read_table(args.frip) if args.frip else {}
    values["FRiP"] = to_number(frip.get("FRiP"))
    values["Total_Peaks"] = to_number(frip.get("Peaks"))

    rows = build_rows(METRICS, values, thresholds)
    write_reports(args.sample, args.prefix, "atac", "ATAC-seq", "ATAC-seq", rows)


if __name__ == "__main__":
    main()
