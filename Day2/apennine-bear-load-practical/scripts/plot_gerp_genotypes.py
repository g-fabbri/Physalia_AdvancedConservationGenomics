#!/usr/bin/env python3
"""Plot GERP distributions and derived-site counts above a GERP threshold."""

import argparse
import csv
import gzip
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(
        path, encoding="utf-8"
    )


def read_samples(path):
    with open(path, encoding="utf-8") as stream:
        return [line.strip() for line in stream if line.strip()]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("gerp", help="Apennine-coordinate BED or BED.GZ with GERP in column 4")
parser.add_argument("gt", help="GenoLoader POP_OUT .gt file")
parser.add_argument("abb", help="ABB sample list")
parser.add_argument("sbb", help="SBB sample list")
parser.add_argument("output_tsv", help="Output per-sample TSV")
parser.add_argument("output_prefix", help="Prefix for the two PDF figures")
args = parser.parse_args()

groups = {"ABB": read_samples(args.abb), "SBB": read_samples(args.sbb)}
samples = {sample: group for group, members in groups.items() for sample in members}
if len(samples) != sum(len(members) for members in groups.values()):
    parser.error("ABB and SBB lists overlap or contain duplicate IDs")

# BED end equals the original one-based VCF position for a one-base interval.
scores = {}
with open_text(args.gerp) as stream:
    for line_number, line in enumerate(stream, 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            parser.error(f"GERP row {line_number} has fewer than four columns")
        try:
            start = int(fields[1])
            end = int(fields[2])
            score = float(fields[3])
        except ValueError:
            parser.error(f"Non-numeric coordinate or score in GERP row {line_number}")
        if end - start != 1:
            parser.error(f"GERP row {line_number} is not a one-base BED interval")
        key = (fields[0], end)
        if key in scores:
            parser.error(f"Duplicate GERP score for {fields[0]}:{end}")
        scores[key] = score

gerp_threshold = 2.0
site_scores = []
# called scored sites, derived copies, heterozygous sites > threshold,
# homozygous-derived sites > threshold
summary = defaultdict(lambda: [0, 0, 0, 0])
gt_rows = 0
trusted_rows = 0
scored_rows = 0
seen_scored_sites = set()

with open(args.gt, encoding="utf-8", newline="") as stream:
    reader = csv.DictReader(stream, delimiter="\t")
    needed = set(samples) | {"scaffold", "position", "flag"}
    missing = needed - set(reader.fieldnames or [])
    if missing:
        parser.error("Missing .gt columns: " + ", ".join(sorted(missing)))

    for row in reader:
        gt_rows += 1
        flag = row["flag"]
        if not (flag.startswith("unfolded") or flag == "allFix"):
            continue
        trusted_rows += 1
        try:
            key = (row["scaffold"], int(row["position"]))
        except ValueError:
            continue
        if key not in scores:
            continue
        if key in seen_scored_sites:
            parser.error(
                f"Duplicate trusted GenoLoader row at {key[0]}:{key[1]}; "
                "resolve duplicate positions before summing genotypes"
            )
        score = scores[key]
        scored_rows += 1
        site_scores.append(score)
        seen_scored_sites.add(key)

        for sample, group in samples.items():
            dosage = row[sample]
            if dosage not in {"0", "1", "2"}:
                continue
            values = summary[(group, sample)]
            values[0] += 1
            values[1] += int(dosage)
            if score > gerp_threshold and dosage == "1":
                values[2] += 1
            elif score > gerp_threshold and dosage == "2":
                values[3] += 1

with open(args.output_tsv, "w", encoding="utf-8", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(
        [
            "population",
            "sample",
            "called_scored_sites",
            "derived_copies",
            "heterozygous_derived_sites_GERP_gt2",
            "homozygous_derived_sites_GERP_gt2",
            "total_derived_sites_GERP_gt2",
            "total_GERP_gt2_sites_per_called_scored_site",
        ]
    )
    for group in ("ABB", "SBB"):
        for sample in groups[group]:
            called, copies, hetero, homo = summary[(group, sample)]
            total = hetero + homo
            writer.writerow(
                [
                    group,
                    sample,
                    called,
                    copies,
                    hetero,
                    homo,
                    total,
                    f"{total / called:.6g}" if called else "NA",
                ]
            )

colors = {"ABB": "#B22222", "SBB": "#4682B4"}

# Figure 1: raw score distribution, including negative values and without
# separating sites by SnpEff impact.
distribution_pdf = f"{args.output_prefix}_score_distribution.pdf"
figure, axis = plt.subplots(figsize=(7, 5))
axis.hist(site_scores, bins=50, color="#6A7D89", edgecolor="white", linewidth=0.3)
axis.axvline(0, color="black", linestyle="--", linewidth=1)
axis.axvline(
    gerp_threshold, color="#B22222", linestyle=":", linewidth=1.5, label="GERP = 2"
)
axis.set_xlabel("GERP score")
axis.set_ylabel("Number of scored SNPs")
axis.legend(frameon=False)
axis.grid(axis="y", color="0.92", linewidth=0.7)
axis.spines["top"].set_visible(False)
axis.spines["right"].set_visible(False)
figure.suptitle("GERP scores at uniquely mapped Scaffold_25 SNPs")
figure.tight_layout(rect=(0, 0, 1, 0.95))
figure.savefig(distribution_pdf)
plt.close(figure)

# Figure 2: derived-site counts across all annotations, with GERP > 2.
derived_sites_pdf = f"{args.output_prefix}_derived_sites_GERP_gt2.pdf"
figure, axis = plt.subplots(figsize=(8, 5.5))
offsets = {"ABB": -0.16, "SBB": 0.16}
categories = (
    ("Heterozygous", 2),
    ("Homozygous derived", 3),
    ("Total", None),
)
for category_index, (_label, metric_index) in enumerate(categories):
    for group in ("ABB", "SBB"):
        if metric_index is None:
            values = [
                summary[(group, sample)][2] + summary[(group, sample)][3]
                for sample in groups[group]
            ]
        else:
            values = [summary[(group, sample)][metric_index] for sample in groups[group]]
        center = category_index + offsets[group]
        x_values = [
            center + (index - (len(values) - 1) / 2) * 0.025
            for index in range(len(values))
        ]
        axis.scatter(
            x_values,
            values,
            s=30,
            alpha=0.8,
            color=colors[group],
            edgecolor="white",
            linewidth=0.4,
            label=group if category_index == 0 else None,
            zorder=2,
        )
        mean = sum(values) / len(values)
        axis.plot(
            [center - 0.1, center + 0.1],
            [mean, mean],
            color=colors[group],
            linewidth=2.3,
            zorder=3,
        )
axis.set_xticks(range(len(categories)), [label for label, _ in categories])
axis.set_ylabel("Number of derived sites with GERP > 2")
axis.grid(axis="y", color="0.92", linewidth=0.7)
axis.spines["top"].set_visible(False)
axis.spines["right"].set_visible(False)
axis.legend(frameon=False)
figure.suptitle("Derived sites at constrained positions on Scaffold_25", y=0.96)
figure.text(
    0.5,
    0.01,
    "Total sites = heterozygous sites + homozygous-derived sites; points are individuals.",
    ha="center",
    fontsize=9,
)
figure.tight_layout(rect=(0, 0.05, 1, 0.93))
figure.savefig(derived_sites_pdf)
plt.close(figure)

print(f"Read {len(scores)} unique GERP positions")
print(f"Read {gt_rows} GenoLoader rows; retained {trusted_rows} trusted rows")
print(f"Matched {len(seen_scored_sites)} unique trusted sites to GERP scores")
print(f"Wrote {args.output_tsv}")
print(f"Wrote {distribution_pdf}")
print(f"Wrote {derived_sites_pdf}")
