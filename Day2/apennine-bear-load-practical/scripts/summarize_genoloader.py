#!/usr/bin/env python3
"""Summarize and plot per-sample derived-allele proxies from GenoLoader."""

import argparse
import csv
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_samples(path):
    with open(path, encoding="utf-8") as stream:
        return [line.strip() for line in stream if line.strip()]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("gt", help="GenoLoader POP_OUT .gt file")
parser.add_argument("abb", help="ABB sample list")
parser.add_argument("sbb", help="SBB sample list")
parser.add_argument("output_tsv", help="Output per-sample TSV")
parser.add_argument("output_pdf", help="Output comparison PDF")
args = parser.parse_args()

groups = {"ABB": read_samples(args.abb), "SBB": read_samples(args.sbb)}
samples = {sample: group for group, members in groups.items() for sample in members}
if len(samples) != sum(map(len, groups.values())):
    parser.error("ABB and SBB lists overlap or contain duplicate IDs")

counts = defaultdict(lambda: [0, 0, 0])  # called, derived copies, hom-derived
rows = 0
trusted = 0
with open(args.gt, encoding="utf-8", newline="") as stream:
    reader = csv.DictReader(stream, delimiter="\t")
    needed = set(samples) | {"flag", "effect", "vartype"}
    missing = needed - set(reader.fieldnames or [])
    if missing:
        parser.error("Missing .gt columns: " + ", ".join(sorted(missing)))

    for row in reader:
        rows += 1
        flag = row["flag"]
        # Outgroup monomorphic ('unfolded...') or every group fixed identically.
        # Exclude ingroup-major-allele fallbacks and polymorphic-outgroup flags.
        if not (flag.startswith("unfolded") or flag == "allFix"):
            continue
        trusted += 1

        effect = row["effect"].upper()
        vartype = row["vartype"].lower()
        if effect == "HIGH":
            category = "HIGH"
        elif "missense" in vartype:
            category = "missense"
        elif "synonymous" in vartype:
            category = "synonymous"
        else:
            continue

        for sample in samples:
            dosage = row[sample]
            if dosage not in {"0", "1", "2"}:
                continue
            values = counts[(samples[sample], sample, category)]
            values[0] += 1
            values[1] += int(dosage)
            values[2] += dosage == "2"

summary = []
with open(args.output_tsv, "w", encoding="utf-8", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(["population", "sample", "category", "called_sites",
                     "derived_copies", "homozygous_derived",
                     "derived_copies_per_called_site",
                     "homozygous_derived_per_called_site"])
    for group in ("ABB", "SBB"):
        for sample in groups[group]:
            for category in ("HIGH", "missense", "synonymous"):
                called, copies, homo = counts[(group, sample, category)]
                copy_ratio = copies / called if called else None
                homo_ratio = homo / called if called else None
                writer.writerow([
                    group, sample, category, called, copies, homo,
                    f"{copy_ratio:.6g}" if copy_ratio is not None else "NA",
                    f"{homo_ratio:.6g}" if homo_ratio is not None else "NA",
                ])
                summary.append({
                    "population": group,
                    "sample": sample,
                    "category": category,
                    "copy_ratio": copy_ratio,
                    "homo_ratio": homo_ratio,
                })


def plot_metric(axis, key, ylabel):
    categories = ("HIGH", "missense", "synonymous")
    colors = {"ABB": "#B22222", "SBB": "#4682B4"}
    offsets = {"ABB": -0.16, "SBB": 0.16}
    seen_labels = set()

    for category_index, category in enumerate(categories):
        for group in ("ABB", "SBB"):
            values = [
                row[key] for row in summary
                if row["category"] == category
                and row["population"] == group
                and row[key] is not None
            ]
            if not values:
                continue
            center = category_index + offsets[group]
            if len(values) == 1:
                x_values = [center]
            else:
                x_values = [
                    center + (index - (len(values) - 1) / 2) * 0.025
                    for index in range(len(values))
                ]
            label = group if group not in seen_labels else None
            axis.scatter(
                x_values, values, s=34, alpha=0.8, color=colors[group],
                edgecolor="white", linewidth=0.4, label=label, zorder=2,
            )
            seen_labels.add(group)
            mean = sum(values) / len(values)
            axis.plot(
                [center - 0.10, center + 0.10], [mean, mean],
                color=colors[group], linewidth=2.4, zorder=3,
            )

    axis.set_xticks(range(len(categories)), categories)
    axis.set_ylabel(ylabel)
    axis.grid(axis="y", color="0.9", linewidth=0.8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


figure, axes = plt.subplots(1, 2, figsize=(10, 4.8))
plot_metric(axes[0], "copy_ratio", "Derived copies per called site")
plot_metric(axes[1], "homo_ratio", "Homozygous-derived sites per called site")
handles, labels = axes[0].get_legend_handles_labels()
figure.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
figure.suptitle("Scaffold_25 derived-allele burden proxies", y=0.94)
figure.text(
    0.5, 0.015,
    "Points are individuals; horizontal bars are population means.",
    ha="center", fontsize=9,
)
figure.tight_layout(rect=(0, 0.05, 1, 0.88))
figure.savefig(args.output_pdf)
plt.close(figure)

print(f"Read {rows} loci; retained {trusted} with conservative outgroup flags")
print(f"Wrote {args.output_tsv}")
print(f"Wrote {args.output_pdf}")
