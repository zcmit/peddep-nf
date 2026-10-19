#!/usr/bin/env python3
"""PEDDEP graded QC report for one ChIP-seq or CUT&RUN sample (case or control).

The metrics are those of the PEDDEP/SEAseq report (raw reads, base quality,
sequence diversity, alignment, NRF/PBC, NSC/RSC, FRiP, peaks, ROSE regions),
graded -2..+2 with cutoffs from a JSON file; the overall grade is the mean.
For control libraries (Input / IgG) enrichment metrics are reported but not
graded, since a good control should show no enrichment.

Writes <prefix>.chip_qc_report.txt, .html and <prefix>.chip_qc_metrics.json.
"""
import argparse
import json

from qc_common import (fastqc_summary, read_trim_reports, read_bowtie2_logs, read_picard_metrics,
                       read_table, read_flagstat_reads, read_spp, count_table_rows, to_number,
                       build_rows, write_reports)

METRICS = [
    ("Raw_Reads", "Reads (single-end) or read pairs (paired-end) sequenced, all lanes", "FastQC (raw)"),
    ("Mean_Read_Length", "Mean raw read length (bp)", "FastQC (raw)"),
    ("Base_Quality", "Per-base sequence quality, worst of all raw FASTQs", "FastQC (raw)"),
    ("Sequence_Diversity", "Overrepresented sequences, worst of all raw FASTQs", "FastQC (raw)"),
    ("Pct_Reads_Q30", "% reads with mean base quality >= 30 (worst FASTQ)", "FastQC (raw)"),
    ("Adapter_Trimmed_Pct", "% reads with adapter trimmed", "Trim Galore"),
    ("Aligned_Percent", "% reads aligned", "bowtie2"),
    ("Duplication_Pct", "% duplicate reads (marked; removed for ChIP and controls)", "Picard MarkDuplicates"),
    ("NRF", "Non-redundant fraction", "library complexity"),
    ("PBC1", "PCR bottleneck coefficient 1 (once / distinct)", "library complexity"),
    ("PBC2", "PCR bottleneck coefficient 2 (once / twice)", "library complexity"),
    ("Final_Reads", "Reads (SE) or pairs (PE) kept after filtering", "samtools flagstat"),
    ("Fragment_Length", "Fragment length (bp): median insert size (PE) or cross-correlation estimate (SE)", "phantompeakqualtools"),
    ("NSC", "Normalized strand cross-correlation", "phantompeakqualtools"),
    ("RSC", "Relative strand cross-correlation", "phantompeakqualtools"),
    ("Peak_Type", "Peak set used for FRiP / Total_Peaks", "MACS3"),
    ("FRiP", "Fraction of reads in peaks", "MACS3 + samtools"),
    ("Total_Peaks", "Peaks called vs control (unique intervals)", "MACS3"),
    ("Stitched_Regions", "ROSE stitched regions (ROSE marks only)", "ROSE"),
    ("SE_Like_Regions", "ROSE super-enhancer-like regions (ROSE marks only)", "ROSE"),
]
CONTROL_UNGRADED = ("NSC", "RSC", "FRiP", "Total_Peaks", "Stitched_Regions", "SE_Like_Regions")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--assay-title", default="ChIP-seq")
    parser.add_argument("--single-end", action="store_true")
    parser.add_argument("--control", action="store_true", help="this library is an Input/IgG control")
    parser.add_argument("--peak-type", default="narrow")
    parser.add_argument("--fastqc", nargs="+", required=True)
    parser.add_argument("--trim-logs", nargs="*", default=[])
    parser.add_argument("--bowtie2-logs", nargs="*", default=[])
    parser.add_argument("--markdup-metrics")
    parser.add_argument("--complexity")
    parser.add_argument("--filtered-flagstat")
    parser.add_argument("--spp")
    parser.add_argument("--frip")
    parser.add_argument("--rose-stitched")
    parser.add_argument("--rose-super")
    args = parser.parse_args()

    with open(args.thresholds) as handle:
        thresholds = json.load(handle)

    fq = fastqc_summary(args.fastqc, single_end=args.single_end)
    values = {
        "Raw_Reads": fq["raw"],
        "Mean_Read_Length": fq["mean_length"],
        "Base_Quality": fq["base_quality"],
        "Sequence_Diversity": fq["diversity"],
        "Pct_Reads_Q30": fq["q30"],
        "Adapter_Trimmed_Pct": read_trim_reports(args.trim_logs) if args.trim_logs else None,
        "Aligned_Percent": read_bowtie2_logs(args.bowtie2_logs) if args.bowtie2_logs else None,
        "Duplication_Pct": read_picard_metrics(args.markdup_metrics) if args.markdup_metrics else None,
        "Final_Reads": read_flagstat_reads(args.filtered_flagstat) if args.filtered_flagstat else None,
        "Peak_Type": None if args.control else args.peak_type,
        "Stitched_Regions": count_table_rows(args.rose_stitched),
        "SE_Like_Regions": count_table_rows(args.rose_super),
    }

    complexity = read_table(args.complexity) if args.complexity else {}
    for metric in ("NRF", "PBC1", "PBC2"):
        values[metric] = to_number(complexity.get(metric))

    if args.spp:
        values["Fragment_Length"], values["NSC"], values["RSC"] = read_spp(args.spp)

    frip = read_table(args.frip) if args.frip else {}
    values["FRiP"] = to_number(frip.get("FRiP"))
    values["Total_Peaks"] = to_number(frip.get("Peaks"))

    rows = build_rows(METRICS, values, thresholds, ungraded=CONTROL_UNGRADED if args.control else ())
    write_reports(args.sample, args.prefix, "chip", args.assay_title, args.assay_title, rows)


if __name__ == "__main__":
    main()
