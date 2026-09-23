#!/usr/bin/env python3
"""Plot GERP distributions and constraint-weighted derived genotypes."""

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

impact_classes = ("MODIFIER", "LOW", "MODERATE", "HIGH")
site_scores = defaultdict(list)
# called sites, derived copies, heterozygous score, homozygous-copy score
summary = defaultdict(lambda: [0, 0, 0.0, 0.0])
gt_rows = 0
trusted_rows = 0
scored_rows = 0
seen_scored_sites = set()

with open(args.gt, encoding="utf-8", newline="") as stream:
    reader = csv.DictReader(stream, delimiter="\t")
    needed = set(samples) | {"scaffold", "position", "effect", "flag"}
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
        effect = row["effect"].upper()
        if effect not in impact_classes:
            continue
        if key in seen_scored_sites:
            parser.error(
                f"Duplicate trusted GenoLoader row at {key[0]}:{key[1]}; "
                "resolve duplicate positions before summing genotypes"
            )
        score = scores[key]
        scored_rows += 1
        site_scores[effect].append(score)
        seen_scored_sites.add(key)

        # Positive GERP is treated as constraint evidence. Negative values are
        # retained in the distribution plot but set to zero for the derived
        # constraint-weighted genotype summary.
        positive_score = max(score, 0.0)
        for sample, group in samples.items():
            dosage = row[sample]
            if dosage not in {"0", "1", "2"}:
                continue
            values = summary[(group, sample, effect)]
            values[0] += 1
            values[1] += int(dosage)
            if dosage == "1":
                values[2] += positive_score
            elif dosage == "2":
                values[3] += 2 * positive_score

with open(args.output_tsv, "w", encoding="utf-8", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(
        [
            "population",
            "sample",
            "impact",
            "called_scored_sites",
            "derived_copies",
            "heterozygous_positive_GERP",
            "homozygous_derived_positive_GERP",
            "total_derived_positive_GERP",
            "total_per_called_scored_site",
        ]
    )
    for group in ("ABB", "SBB"):
        for sample in groups[group]:
            for effect in impact_classes:
                called, copies, hetero, homo = summary[(group, sample, effect)]
                total = hetero + homo
                writer.writerow(
                    [
                        group,
                        sample,
                        effect,
                        called,
                        copies,
                        f"{hetero:.6g}",
                        f"{homo:.6g}",
                        f"{total:.6g}",
                        f"{total / called:.6g}" if called else "NA",
                    ]
                )

colors = {"ABB": "#B22222", "SBB": "#4682B4"}

# Figure 1: raw score distribution, including negative values.
distribution_pdf = f"{args.output_prefix}_score_distribution.pdf"
all_scores = [score for effect in impact_classes for score in site_scores[effect]]
figure, axes = plt.subplots(1, 2, figsize=(10, 4.8))
axes[0].hist(all_scores, bins=50, color="#6A7D89", edgecolor="white", linewidth=0.3)
axes[0].axvline(0, color="black", linestyle="--", linewidth=1)
axes[0].set_xlabel("GERP score")
axes[0].set_ylabel("Number of scored SNPs")
axes[0].set_title("All retained sites")

box_values = [site_scores[effect] for effect in impact_classes]
axes[1].boxplot(box_values, labels=impact_classes, showfliers=False)
axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1].set_ylabel("GERP score")
axes[1].set_title("Scores by SnpEff impact")
axes[1].tick_params(axis="x", rotation=25)
for axis in axes:
    axis.grid(axis="y", color="0.92", linewidth=0.7)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
figure.suptitle("GERP scores at uniquely mapped Scaffold_25 SNPs")
figure.tight_layout(rect=(0, 0, 1, 0.94))
figure.savefig(distribution_pdf)
plt.close(figure)


def plot_points(axis, metric_index, ylabel):
    offsets = {"ABB": -0.16, "SBB": 0.16}
    for effect_index, effect in enumerate(impact_classes):
        for group in ("ABB", "SBB"):
            values = [
                summary[(group, sample, effect)][metric_index]
                for sample in groups[group]
            ]
            center = effect_index + offsets[group]
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
                label=group if effect_index == 0 else None,
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
    axis.set_xticks(range(len(impact_classes)), impact_classes, rotation=25)
    axis.set_ylabel(ylabel)
    axis.grid(axis="y", color="0.92", linewidth=0.7)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


# Figure 2: positive-GERP weighted derived-allele contributions.
weighted_pdf = f"{args.output_prefix}_derived_genotype_scores.pdf"
figure, axes = plt.subplots(1, 3, figsize=(13, 4.8))
plot_points(axes[0], 2, "Heterozygous contribution")
plot_points(axes[1], 3, "Homozygous-derived contribution")

# Total is derived from the two stored components rather than a list index.
offsets = {"ABB": -0.16, "SBB": 0.16}
for effect_index, effect in enumerate(impact_classes):
    for group in ("ABB", "SBB"):
        values = [
            summary[(group, sample, effect)][2] + summary[(group, sample, effect)][3]
            for sample in groups[group]
        ]
        center = effect_index + offsets[group]
        x_values = [
            center + (index - (len(values) - 1) / 2) * 0.025
            for index in range(len(values))
        ]
        axes[2].scatter(
            x_values,
            values,
            s=30,
            alpha=0.8,
            color=colors[group],
            edgecolor="white",
            linewidth=0.4,
            zorder=2,
        )
        mean = sum(values) / len(values)
        axes[2].plot(
            [center - 0.1, center + 0.1],
            [mean, mean],
            color=colors[group],
            linewidth=2.3,
            zorder=3,
        )
axes[2].set_xticks(range(len(impact_classes)), impact_classes, rotation=25)
axes[2].set_ylabel("Total derived positive-GERP score")
axes[2].grid(axis="y", color="0.92", linewidth=0.7)
axes[2].spines["top"].set_visible(False)
axes[2].spines["right"].set_visible(False)

handles, labels = axes[0].get_legend_handles_labels()
figure.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
figure.suptitle("Constraint-weighted derived genotypes on Scaffold_25", y=0.95)
figure.text(
    0.5,
    0.01,
    "Positive GERP only; homozygous-derived genotypes contribute two allele copies.",
    ha="center",
    fontsize=9,
)
figure.tight_layout(rect=(0, 0.05, 1, 0.89))
figure.savefig(weighted_pdf)
plt.close(figure)

print(f"Read {len(scores)} unique GERP positions")
print(f"Read {gt_rows} GenoLoader rows; retained {trusted_rows} trusted rows")
print(f"Matched {len(seen_scored_sites)} unique trusted sites to GERP scores")
print(f"Wrote {args.output_tsv}")
print(f"Wrote {distribution_pdf}")
print(f"Wrote {weighted_pdf}"
