#!/usr/bin/env python3
"""Fragment-length QC for paired-end ATAC-seq.

Counts every proper pair once (read 1, |TLEN|) and reports:
  Median_Fragment     median fragment length
  NFR_Fraction        fraction of fragments < 100 bp (nucleosome-free)
  Mono_Fraction       fraction of fragments 180-247 bp (mono-nucleosome)
  Nucleosome_Pattern  pass: both a nucleosome-free peak and a mono-nucleosome
                      peak are visible; warn: only one; fail: neither
Writes <prefix>.fragment_qc.txt, <prefix>.fragment_lengths.txt and a PNG.
"""
import argparse
import sys

import pysam

MAX_LEN = 1000


def smooth(values, width=11):
    half = width // 2
    out = []
    for i in range(len(values)):
        lo, hi = max(0, i - half), min(len(values), i + half + 1)
        out.append(sum(values[lo:hi]) / float(hi - lo))
    return out


def argmax(values, lo, hi):
    best = lo
    for i in range(lo, hi):
        if values[i] > values[best]:
            best = i
    return best


def nucleosome_pattern(hist, total):
    if total < 1000:
        return "fail"
    s = smooth(hist)
    nfr_i = argmax(s, 30, 121)
    mono_i = argmax(s, 150, 301)
    has_nfr = 35 < nfr_i < 120
    has_mono = False
    if 155 < mono_i < 300:
        valley = min(s[min(nfr_i, mono_i):mono_i])
        has_mono = s[mono_i] > 1.1 * valley
    if has_nfr and has_mono:
        return "pass"
    if has_nfr or has_mono:
        return "warn"
    return "fail"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    hist = [0] * (MAX_LEN + 1)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if (read.is_read1 and read.is_proper_pair and not read.is_secondary
                    and not read.is_supplementary and read.template_length != 0):
                hist[min(abs(read.template_length), MAX_LEN)] += 1

    total = sum(hist)
    median = 0
    if total:
        running = 0
        for length, n in enumerate(hist):
            running += n
            if running * 2 >= total:
                median = length
                break
    nfr = sum(hist[:100]) / total if total else 0.0
    mono = sum(hist[180:248]) / total if total else 0.0
    pattern = nucleosome_pattern(hist, total)

    with open(args.prefix + ".fragment_qc.txt", "w") as out:
        out.write("Sample\tFragments\tMedian_Fragment\tNFR_Fraction\tMono_Fraction\tNucleosome_Pattern\n")
        out.write("%s\t%d\t%d\t%.4f\t%.4f\t%s\n" % (args.sample, total, median, nfr, mono, pattern))
    with open(args.prefix + ".fragment_lengths.txt", "w") as out:
        out.write("Size\tOccurrences\n")
        for length, n in enumerate(hist):
            if n:
                out.write("%d\t%d\n" % (length, n))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(range(len(hist)), [100.0 * n / total if total else 0 for n in hist], color="#2b6cb0", lw=1.2)
        ax.axvspan(0, 100, color="#a0aec0", alpha=0.15, lw=0)
        ax.axvspan(180, 247, color="#a0aec0", alpha=0.3, lw=0)
        ax.set_xlim(0, 800)
        ax.set_xlabel("Fragment length (bp)")
        ax.set_ylabel("% of fragments")
        ax.set_title("%s  (median %d bp, NFR %.1f%%, pattern: %s)" % (args.sample, median, 100 * nfr, pattern), fontsize=9)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        fig.tight_layout()
        fig.savefig(args.prefix + ".fragment_sizes.png", dpi=150)
    except ImportError:
        print("matplotlib not available; writing an empty PNG placeholder", file=sys.stderr)
        open(args.prefix + ".fragment_sizes.png", "wb").close()


if __name__ == "__main__":
    main()
