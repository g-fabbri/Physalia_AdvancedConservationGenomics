# Optional Part 3B — currentNe2: a contemporary LD estimate

This extension is **outside the timed practical**. GONE2 uses LD at different recombination distances to fit a recent trajectory. currentNe2 also uses LD, but reports a **contemporary** effective-size estimate. Because the methods share a broad signal, agreement is not independent confirmation. We reuse the single-scaffold GONE2 input here only to learn the workflow; a substantive currentNe2 analysis should use many chromosomes, more independent markers, and assess sample size and relatedness.

Install and test currentNe2 using the [software setup](../software/README.md). Run these commands from the Day 1 course directory after completing [GONE2](02-gone2.md).

## Step 1 — Reuse the ABB and SBB genotypes

**Purpose:** put each population's existing PED/MAP pair in its own currentNe2 input location without adding a new SNP filter.

```bash
CHROM=Scaffold_25
GONE2_INPUT=results/gone2/input
OUTDIR=results/currentne2
mkdir -p "$OUTDIR"

for POPULATION in ABB SBB; do
  cp "$GONE2_INPUT/${POPULATION}_${CHROM}.ped" "$OUTDIR/${POPULATION}_${CHROM}.ped"
  cp "$GONE2_INPUT/${POPULATION}_${CHROM}.map" "$OUTDIR/${POPULATION}_${CHROM}.map"
  printf '%s individuals: ' "$POPULATION"
  wc -l < "$OUTDIR/${POPULATION}_${CHROM}.ped"
  printf '%s markers: ' "$POPULATION"
  wc -l < "$OUTDIR/${POPULATION}_${CHROM}.map"
done
```

**Expected:** 10 ABB and 8 SBB individuals, with a positive marker count for each. Both MAPs describe the same physical scaffold, numbered `1`. Do **not** concatenate ABB and SBB: this step estimates each population separately.

## Step 2 — Estimate contemporary size

**Purpose:** apply the same physical-to-genetic-map assumption to each population. `-r 1` means 1 cM/Mb; `-t 2` uses two threads. We do not use `-x`, which would introduce a metapopulation model, or supply a full-genome length for this one-scaffold demonstration.

```bash
CURRENTNE2_BIN="$PWD/software/currentNe2/currentne2"

for POPULATION in ABB SBB; do
  "$CURRENTNE2_BIN" -r 1 -t 2 \
    "$OUTDIR/${POPULATION}_${CHROM}.ped"
done

ls -lh "$OUTDIR"/*_currentNe2_OUTPUT.txt
```

**Expected:** one output text file per population. currentNe2's documented default is to write beside the input with the `_currentNe2_OUTPUT.txt` suffix. If the run fails, inspect the complete terminal message before comparing estimates.

## Step 3 — Read and discuss

```bash
cat "$OUTDIR/ABB_${CHROM}_currentNe2_OUTPUT.txt"
cat "$OUTDIR/SBB_${CHROM}_currentNe2_OUTPUT.txt"
```

Record the reported contemporary estimate and any uncertainty or diagnostic information. Ask whether ABB and SBB differ, whether the most recent part of GONE2 points in the same direction, and whether sampling noise or shared LD assumptions could explain a mismatch. **Do not present either one-scaffold value as the current effective size of the whole population.** currentNe2 can also use information from pairs of *different* chromosomes when a suitable genome-wide input is available; this classroom input cannot supply that comparison.

Return to the [synthesis](04-synthesis.md).
