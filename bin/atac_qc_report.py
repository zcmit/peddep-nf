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
import html
import json
import os
import re
import sys
import zipfile

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

# Colour-blind palette used by the PEDDEP/SEAseq reports
COLORS = {-2: "#C53440", -1: "#E28B51", 0: "#CBBC6B", 1: "#97BEBF", 2: "#A77393"}
STATUS_SCORES = {"pass": 2, "warn": 0, "fail": -2}
STATUS_RANK = {"fail": 0, "warn": 1, "pass": 2}


# ----------------------------------------------------------------------------- parsers

def read_fastqc(zip_path):
    """Return dict with total, statuses, q30 fraction and mean length."""
    with zipfile.ZipFile(zip_path) as zf:
        name = next(n for n in zf.namelist() if n.endswith("fastqc_data.txt"))
        text = zf.read(name).decode()
    data = {"file": os.path.basename(zip_path), "statuses": {}}
    module = None
    qual_counts = {}
    len_counts = {}
    for line in text.splitlines():
        if line.startswith(">>END_MODULE"):
            module = None
            continue
        if line.startswith(">>"):
            parts = line[2:].split("\t")
            module = parts[0]
            if len(parts) > 1:
                data["statuses"][module] = parts[1]
            continue
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if module == "Basic Statistics" and fields[0] == "Total Sequences":
            data["total"] = int(fields[1])
        elif module == "Per sequence quality scores" and len(fields) == 2:
            qual_counts[float(fields[0])] = float(fields[1])
        elif module == "Sequence Length Distribution" and len(fields) == 2:
            bounds = [float(x) for x in fields[0].split("-")]
            len_counts[sum(bounds) / len(bounds)] = float(fields[1])
    total_q = sum(qual_counts.values())
    data["q30"] = 100.0 * sum(c for q, c in qual_counts.items() if q >= 30) / total_q if total_q else None
    total_l = sum(len_counts.values())
    data["mean_length"] = sum(l * c for l, c in len_counts.items()) / total_l if total_l else None
    return data


def read_trim_reports(paths):
    processed = with_adapter = 0
    for path in paths:
        with open(path) as handle:
            for line in handle:
                m = re.match(r"Total reads processed:\s+([\d,]+)", line)
                if m:
                    processed += int(m.group(1).replace(",", ""))
                m = re.match(r"Reads with adapters:\s+([\d,]+)", line)
                if m:
                    with_adapter += int(m.group(1).replace(",", ""))
    return 100.0 * with_adapter / processed if processed else None


def read_bowtie2_logs(paths):
    pairs = 0
    aligned = 0.0
    for path in paths:
        n = rate = None
        with open(path) as handle:
            for line in handle:
                m = re.match(r"\s*(\d+) reads; of these:", line)
                if m:
                    n = int(m.group(1))
                m = re.match(r"\s*([\d.]+)% overall alignment rate", line)
                if m:
                    rate = float(m.group(1))
        if n is not None and rate is not None:
            pairs += n
            aligned += n * rate
    return aligned / pairs if pairs else None


def read_idxstats(path, mito_names):
    total = mito = 0
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] == "*":
                continue
            mapped = int(fields[2])
            total += mapped
            if fields[0] in mito_names:
                mito += mapped
    return 100.0 * mito / total if total else None


def read_picard_metrics(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            header = lines[i + 1].split("\t")
            values = lines[i + 2].split("\t")
            row = dict(zip(header, values))
            value = row.get("PERCENT_DUPLICATION")
            return 100.0 * float(value) if value not in (None, "", "?") else None
    return None


def read_table(path):
    """First data row of a 2-line TSV as a dict."""
    with open(path) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if len(rows) < 2:
        return {}
    return dict(zip(rows[0], rows[1]))


def read_flagstat_read1(path):
    with open(path) as handle:
        for line in handle:
            m = re.match(r"(\d+) \+ \d+ read1", line)
            if m:
                return int(m.group(1))
    return None


def read_ataqv(path):
    with open(path) as handle:
        data = json.load(handle)
    entry = data[0] if isinstance(data, list) else data
    metrics = entry.get("metrics", entry)
    value = metrics.get("tss_enrichment")
    return float(value) if value not in (None, "") else None


def to_number(value):
    if value in (None, "", "NA", "nan", "inf"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


# ----------------------------------------------------------------------------- grading

def score_metric(value, rule):
    if value is None or rule is None:
        return None
    kind = rule.get("type", "higher")
    if kind == "status":
        return STATUS_SCORES.get(str(value).lower(), -2)
    cuts = rule["cutoffs"]
    if kind == "higher":
        score = -2
        for s, cut in zip((-1, 0, 1, 2), cuts):
            if value >= cut:
                score = s
        return score
    if kind == "lower":
        score = -2
        for s, cut in zip((-1, 0, 1, 2), cuts):
            if value <= cut:
                score = s
        return score
    raise ValueError("unknown threshold type %r" % kind)


def overall_grade(scores):
    if not scores:
        return "NA", None, None
    mean = sum(scores) / float(len(scores))
    rounded = round(mean, 1)
    label, level = "POOR", -2
    for threshold, name, lvl in ((-1, "BELOW_AVERAGE", -1), (0, "AVERAGE", 0), (1, "GOOD", 1), (2, "EXCELLENT", 2)):
        if rounded >= threshold:
            label, level = name, lvl
    return label, level, mean


def format_value(name, value):
    if value is None:
        return "NA"
    if isinstance(value, str):
        return value
    if name in ("Raw_Read_Pairs", "Final_Read_Pairs", "Total_Peaks", "Median_Fragment"):
        return "%d" % value
    if name in ("NRF", "PBC1", "PBC2", "FRiP", "NFR_Fraction"):
        return "%.4f" % value
    return "%.2f" % value


def describe_rule(rule):
    if not rule:
        return "informational"
    if rule.get("type") == "status":
        return "pass = 2, warn = 0, fail = -2"
    sign = ">=" if rule.get("type", "higher") == "higher" else "<="
    return ", ".join("%s %s for %+d" % (sign, c, s) for s, c in zip((-1, 0, 1, 2), rule["cutoffs"]))


# ----------------------------------------------------------------------------- output

def html_page(sample, label, level, rows):
    head = "".join("<th>%s</th>" % html.escape(r["name"]) for r in rows)
    cells = "".join(
        '<td style="background:%s">%s</td>' % (COLORS[r["score"]], html.escape(r["display"]))
        if r["score"] is not None else "<td>%s</td>" % html.escape(r["display"])
        for r in rows
    )
    overall_style = ' style="background:%s"' % COLORS[level] if level is not None else ""
    details = "".join(
        "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
            html.escape(r["name"]), html.escape(r["display"]),
            "" if r["score"] is None else "%+d" % r["score"],
            html.escape(r["rule"]), html.escape(r["source"]))
        for r in rows
    )
    legend = "".join(
        '<span class="key" style="background:%s">%s</span>' % (COLORS[s], n)
        for s, n in ((2, "EXCELLENT / +2"), (1, "GOOD / +1"), (0, "AVERAGE / 0"),
                     (-1, "BELOW_AVERAGE / -1"), (-2, "POOR / -2"))
    )
    return """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>{sample} ATAC-seq QC</title>
<style>
body {{ font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif; margin: 24px; color: #1a202c; }}
h1 {{ font-size: 20px; margin: 0 0 4px; }}
p.sub {{ color: #4a5568; margin: 0 0 16px; font-size: 13px; }}
table {{ border-collapse: collapse; font-size: 12px; margin-bottom: 20px; }}
th, td {{ border: 1px solid #cbd5e0; padding: 4px 8px; text-align: center; }}
th {{ background: #edf2f7; }}
.scroll {{ overflow-x: auto; }}
.key {{ display: inline-block; padding: 2px 8px; margin-right: 6px; font-size: 12px; border-radius: 3px; }}
table.details td {{ text-align: left; }}
</style></head><body>
<h1>{sample} &mdash; ATAC-seq QC report</h1>
<p class="sub">PEDDEP grading: each graded metric scores -2 to +2; the overall grade is the mean score.</p>
<div class="scroll"><table><tr><th>Sample_Name</th><th>Overall_Quality</th>{head}</tr>
<tr><td>{sample}</td><td{overall_style}><b>{label}</b></td>{cells}</tr></table></div>
<p>{legend}</p>
<h2 style="font-size:16px">Metric details</h2>
<table class="details"><tr><th>Metric</th><th>Value</th><th>Score</th><th>Thresholds</th><th>Source</th></tr>{details}</table>
</body></html>
""".format(sample=html.escape(sample), head=head, cells=cells, overall_style=overall_style,
           label=html.escape(label), legend=legend, details=details)


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

    values = {}

    fastqc = [read_fastqc(p) for p in args.fastqc]
    r1 = [d for d in fastqc if re.search(r"_R1(_\d+)?_fastqc\.zip$", d["file"])]
    if not r1:
        r1 = fastqc[: max(1, len(fastqc) // 2)]
    values["Raw_Read_Pairs"] = float(sum(d.get("total", 0) for d in r1))
    lengths = [d["mean_length"] for d in fastqc if d.get("mean_length") is not None]
    values["Mean_Read_Length"] = sum(lengths) / len(lengths) if lengths else None
    for metric, module in (("Base_Quality", "Per base sequence quality"),
                           ("Sequence_Diversity", "Overrepresented sequences")):
        statuses = [d["statuses"].get(module) for d in fastqc if d["statuses"].get(module)]
        values[metric] = min(statuses, key=lambda s: STATUS_RANK.get(s, 0)) if statuses else None
    q30 = [d["q30"] for d in fastqc if d.get("q30") is not None]
    values["Pct_Reads_Q30"] = min(q30) if q30 else None

    values["Adapter_Trimmed_Pct"] = read_trim_reports(args.trim_logs) if args.trim_logs else None
    values["Aligned_Percent"] = read_bowtie2_logs(args.bowtie2_logs) if args.bowtie2_logs else None
    values["Mito_Percent"] = read_idxstats(args.idxstats, set(args.mito_name.split(","))) if args.idxstats else None
    values["Duplication_Pct"] = read_picard_metrics(args.markdup_metrics) if args.markdup_metrics else None

    complexity = read_table(args.complexity) if args.complexity else {}
    for metric in ("NRF", "PBC1", "PBC2"):
        values[metric] = to_number(complexity.get(metric))

    values["Final_Read_Pairs"] = read_flagstat_read1(args.filtered_flagstat) if args.filtered_flagstat else None

    fragments = read_table(args.fragment_qc) if args.fragment_qc else {}
    values["Median_Fragment"] = to_number(fragments.get("Median_Fragment"))
    values["NFR_Fraction"] = to_number(fragments.get("NFR_Fraction"))
    values["Nucleosome_Pattern"] = fragments.get("Nucleosome_Pattern") or None

    values["TSS_Enrichment"] = read_ataqv(args.ataqv) if args.ataqv else None

    frip = read_table(args.frip) if args.frip else {}
    values["FRiP"] = to_number(frip.get("FRiP"))
    values["Total_Peaks"] = to_number(frip.get("Peaks"))

    rows = []
    for name, description, source in METRICS:
        rule = thresholds.get(name)
        value = values.get(name)
        score = score_metric(value, rule)
        rows.append({
            "name": name, "description": description, "source": source,
            "value": value, "display": format_value(name, value),
            "score": score, "graded": rule is not None, "rule": describe_rule(rule),
        })

    label, level, mean = overall_grade([r["score"] for r in rows if r["score"] is not None])

    with open(args.prefix + ".atac_qc_report.txt", "w") as out:
        out.write("# PEDDEP ATAC-seq QC Report\n")
        out.write("# Overall_Quality = mean score of graded metrics (-2..2): "
                  ">=2 EXCELLENT, >=1 GOOD, >=0 AVERAGE, >=-1 BELOW_AVERAGE, else POOR\n")
        out.write("\t".join(["Sample_Name", "Overall_Quality"] + [r["name"] for r in rows]) + "\n")
        out.write("\t".join([args.sample, label] + [r["display"] for r in rows]) + "\n")

    with open(args.prefix + ".atac_qc_report.html", "w") as out:
        out.write(html_page(args.sample, label, level, rows))

    with open(args.prefix + ".atac_qc_metrics.json", "w") as out:
        json.dump({
            "sample": args.sample,
            "overall": {"label": label, "score": level, "mean": mean},
            "metrics": [{k: r[k] for k in ("name", "value", "display", "score", "graded", "description", "source", "rule")}
                        for r in rows],
        }, out, indent=2)

    print("%s: %s (mean score %s)" % (args.sample, label, "NA" if mean is None else "%.2f" % mean), file=sys.stderr)


if __name__ == "__main__":
    main()
