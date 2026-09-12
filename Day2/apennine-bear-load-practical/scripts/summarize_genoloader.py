#!/usr/bin/env python3
"""Summarize per-sample derived-allele proxies from a GenoLoader .gt table."""

import argparse
import csv
from collections import defaultdict


def read_samples(path):
    with open(path, encoding="utf-8") as stream:
        return [line.strip() for line in stream if line.strip()]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("gt", help="GenoLoader POP_OUT .gt file")
parser.add_argument("abb", help="ABB sample list")
parser.add_argument("sbb", help="SBB sample list")
parser.add_argument("output", help="Output TSV")
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
        elif vartype == "missense":
            category = "missense"
        elif vartype == "synonymous":
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

with open(args.output, "w", encoding="utf-8", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(["population", "sample", "category", "called_sites",
                     "derived_copies", "homozygous_derived",
                     "derived_copies_per_called_site"])
    for group in ("ABB", "SBB"):
        for sample in groups[group]:
            for category in ("HIGH", "missense", "synonymous"):
                called, copies, homo = counts[(group, sample, category)]
                ratio = f"{copies / called:.6g}" if called else "NA"
                writer.writerow([group, sample, category, called, copies,
                                 homo, ratio])

print(f"Read {rows} loci; retained {trusted} with conservative outgroup flags")
print(f"Wrote {args.output}")
