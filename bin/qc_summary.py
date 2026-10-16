#!/usr/bin/env python3
"""Combine per-sample *_qc_metrics.json files (any assay) into one run-level table.

Writes <out-prefix>.tsv, <out-prefix>.html (coloured by score) and an HTML
section for MultiQC (custom content, file name must end in _mqc.html).
"""
import argparse
import html
import json

COLORS = {-2: "#C53440", -1: "#E28B51", 0: "#CBBC6B", 1: "#97BEBF", 2: "#A77393"}


def load(paths):
    samples = []
    for path in paths:
        with open(path) as handle:
            samples.append(json.load(handle))
    return sorted(samples, key=lambda s: s["sample"])


def table_html(samples):
    names = [m["name"] for m in samples[0]["metrics"]]
    head = "<tr><th>Sample_Name</th><th>Overall_Quality</th>%s</tr>" % "".join(
        "<th>%s</th>" % html.escape(n) for n in names)
    body = []
    for s in samples:
        level = s["overall"]["score"]
        overall = '<td style="background:%s"><b>%s</b></td>' % (COLORS[level], html.escape(s["overall"]["label"])) \
            if level is not None else "<td>%s</td>" % html.escape(s["overall"]["label"])
        cells = "".join(
            '<td style="background:%s">%s</td>' % (COLORS[m["score"]], html.escape(m["display"]))
            if m["score"] is not None else "<td>%s</td>" % html.escape(m["display"])
            for m in s["metrics"]
        )
        body.append("<tr><td>%s</td>%s%s</tr>" % (html.escape(s["sample"]), overall, cells))
    return '<table class="peddep-grades">%s%s</table>' % (head, "".join(body))


STYLE = """<style>
table.peddep-grades { border-collapse: collapse; font-size: 12px; }
table.peddep-grades th, table.peddep-grades td { border: 1px solid #cbd5e0; padding: 3px 7px; text-align: center; white-space: nowrap; }
table.peddep-grades th { background: #edf2f7; }
</style>"""

LEGEND = "<p>" + "".join(
    '<span style="display:inline-block;padding:2px 8px;margin-right:6px;background:%s">%s</span>' % (COLORS[s], n)
    for s, n in ((2, "EXCELLENT / +2"), (1, "GOOD / +1"), (0, "AVERAGE / 0"), (-1, "BELOW_AVERAGE / -1"), (-2, "POOR / -2"))
) + "</p>"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--multiqc", required=True)
    parser.add_argument("--title", default="ATAC-seq", help="assay name shown in headings")
    parser.add_argument("--section-id", default="peddep_atac_grades", help="MultiQC section id")
    parser.add_argument("--thresholds-name", default="atac_qc_thresholds.json",
                        help="thresholds file named in the MultiQC description")
    parser.add_argument("metrics", nargs="+")
    args = parser.parse_args()

    samples = load(args.metrics)
    names = [m["name"] for m in samples[0]["metrics"]]

    with open(args.out_prefix + ".tsv", "w") as out:
        out.write("\t".join(["Sample_Name", "Overall_Quality", "Overall_Mean_Score"] + names) + "\n")
        for s in samples:
            mean = s["overall"]["mean"]
            out.write("\t".join([s["sample"], s["overall"]["label"], "NA" if mean is None else "%.2f" % mean]
                                + [m["display"] for m in s["metrics"]]) + "\n")

    table = table_html(samples)
    with open(args.out_prefix + ".html", "w") as out:
        out.write("<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>%s QC summary</title>" % args.title)
        out.write("%s</head><body style=\"font-family:Helvetica,Arial,sans-serif;margin:24px\">" % STYLE)
        out.write("<h1 style=\"font-size:20px\">PEDDEP %s QC summary (%d samples)</h1>" % (args.title, len(samples)))
        out.write("<p>Each graded metric scores -2 to +2; Overall_Quality is the mean score.</p>")
        out.write("<div style=\"overflow-x:auto\">%s</div>%s</body></html>\n" % (table, LEGEND))

    with open(args.multiqc, "w") as out:
        out.write("<!--\nid: '%s'\nsection_name: 'PEDDEP QC grades'\n" % args.section_id)
        out.write("description: 'Each graded metric scores -2 to +2 (thresholds in assets/%s); " % args.thresholds_name)
        out.write("Overall_Quality is the mean score. Uncoloured columns are informational.'\n-->\n")
        out.write(STYLE + '<div style="overflow-x:auto">' + table + "</div>" + LEGEND + "\n")

if __name__ == "__main__":
    main()
