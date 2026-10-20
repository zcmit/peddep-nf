"""Shared helpers for the PEDDEP graded QC reports (ATAC-seq, ChIP-seq, CUT&RUN).

Parsers for FastQC / Trim Galore / bowtie2 / samtools / Picard / ataqv / SPP
outputs, the -2..+2 scoring rules, and the text / HTML / JSON report writers.
Imported by atac_qc_report.py and chip_qc_report.py (same directory).
"""
import html
import json
import os
import re
import sys
import zipfile

COLORS = {-2: "#C53440", -1: "#E28B51", 0: "#CBBC6B", 1: "#97BEBF", 2: "#A77393"}
STATUS_SCORES = {"pass": 2, "warn": 0, "fail": -2}
STATUS_RANK = {"fail": 0, "warn": 1, "pass": 2}
INTEGER_METRICS = {
    "Raw_Read_Pairs", "Raw_Reads", "Final_Read_Pairs", "Final_Reads", "Total_Peaks",
    "Median_Fragment", "Fragment_Length", "Stitched_Regions", "SE_Like_Regions",
}
FOUR_DECIMAL_METRICS = {"NRF", "PBC1", "PBC2", "FRiP", "NFR_Fraction", "NSC", "RSC"}



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
    if name in INTEGER_METRICS:
        return "%d" % value
    if name in FOUR_DECIMAL_METRICS:
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

def html_page(sample, label, level, rows, assay_label):
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
<html lang="en"><head><meta charset="utf-8"><title>{sample} {assay} QC</title>
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
<h1>{sample} &mdash; {assay} QC report</h1>
<p class="sub">PEDDEP grading: each graded metric scores -2 to +2; the overall grade is the mean score.</p>
<div class="scroll"><table><tr><th>Sample_Name</th><th>Overall_Quality</th>{head}</tr>
<tr><td>{sample}</td><td{overall_style}><b>{label}</b></td>{cells}</tr></table></div>
<p>{legend}</p>
<h2 style="font-size:16px">Metric details</h2>
<table class="details"><tr><th>Metric</th><th>Value</th><th>Score</th><th>Thresholds</th><th>Source</th></tr>{details}</table>
</body></html>
""".format(sample=html.escape(sample), head=head, cells=cells, overall_style=overall_style,
           label=html.escape(label), legend=legend, details=details,
           assay=html.escape(assay_label))


def fastqc_summary(paths, single_end=False):
    """Raw read (pair) count, mean length, worst module statuses and worst %Q30 over FASTQs."""
    fastqc = [read_fastqc(p) for p in paths]
    if single_end:
        first = fastqc
    else:
        first = [d for d in fastqc if re.search(r"_R1(_\d+)?_fastqc\.zip$", d["file"])]
        if not first:
            first = fastqc[: max(1, len(fastqc) // 2)]
    out = {"raw": float(sum(d.get("total", 0) for d in first))}
    lengths = [d["mean_length"] for d in fastqc if d.get("mean_length") is not None]
    out["mean_length"] = sum(lengths) / len(lengths) if lengths else None
    for key, module in (("base_quality", "Per base sequence quality"),
                        ("diversity", "Overrepresented sequences")):
        statuses = [d["statuses"].get(module) for d in fastqc if d["statuses"].get(module)]
        out[key] = min(statuses, key=lambda s: STATUS_RANK.get(s, 0)) if statuses else None
    q30 = [d["q30"] for d in fastqc if d.get("q30") is not None]
    out["q30"] = min(q30) if q30 else None
    return out


def read_flagstat_reads(path):
    """Read pairs (read1 count) for paired-end data, else primary mapped reads."""
    pairs = read_flagstat_read1(path)
    if pairs:
        return pairs
    with open(path) as handle:
        for line in handle:
            m = re.match(r"(\d+) \+ \d+ (primary )?mapped", line)
            if m:
                return int(m.group(1))
    return None


def read_spp(path):
    """phantompeakqualtools run_spp.R output: (fragment length, NSC, RSC)."""
    with open(path) as handle:
        fields = handle.readline().rstrip("\n").split("\t")
    if len(fields) < 10:
        return None, None, None
    frag = to_number(fields[2].split(",")[0])
    return frag, to_number(fields[8]), to_number(fields[9])


def count_table_rows(path, header_lines=6):
    """ROSE tables start with a fixed header block (the PEDDEP report subtracts 6 lines)."""
    if not path or not os.path.exists(path):
        return None
    with open(path) as handle:
        n = sum(1 for _ in handle)
    return float(max(0, n - header_lines))


def build_rows(metrics, values, thresholds, ungraded=()):
    rows = []
    for name, description, source in metrics:
        rule = None if name in ungraded else thresholds.get(name)
        value = values.get(name)
        score = score_metric(value, rule)
        rows.append({
            "name": name, "description": description, "source": source,
            "value": value, "display": format_value(name, value),
            "score": score, "graded": rule is not None, "rule": describe_rule(rule),
        })
    return rows


def write_reports(sample, prefix, kind, assay_title, assay_label, rows):
    """<prefix>.<kind>_qc_report.txt / .html and <prefix>.<kind>_qc_metrics.json"""
    label, level, mean = overall_grade([r["score"] for r in rows if r["score"] is not None])
    with open("%s.%s_qc_report.txt" % (prefix, kind), "w") as out:
        out.write("# PEDDEP %s QC Report\n" % assay_title)
        out.write("# Overall_Quality = mean score of graded metrics (-2..2): "
                  ">=2 EXCELLENT, >=1 GOOD, >=0 AVERAGE, >=-1 BELOW_AVERAGE, else POOR\n")
        out.write("\t".join(["Sample_Name", "Overall_Quality"] + [r["name"] for r in rows]) + "\n")
        out.write("\t".join([sample, label] + [r["display"] for r in rows]) + "\n")
    with open("%s.%s_qc_report.html" % (prefix, kind), "w") as out:
        out.write(html_page(sample, label, level, rows, assay_label))
    with open("%s.%s_qc_metrics.json" % (prefix, kind), "w") as out:
        json.dump({
            "sample": sample,
            "overall": {"label": label, "score": level, "mean": mean},
            "metrics": [{k: r[k] for k in ("name", "value", "display", "score", "graded", "description", "source", "rule")}
                        for r in rows],
        }, out, indent=2)
    print("%s: %s (mean score %s)" % (sample, label, "NA" if mean is None else "%.2f" % mean), file=sys.stderr)
