# Optional Part 3B — currentNe2: a contemporary LD estimate

GONE2 uses LD at different recombination distances to fit a recent trajectory. currentNe2 also uses LD, but reports a **contemporary** effective-size estimate. Because the methods share a broad signal, agreement is not independent confirmation.

There is an important software limitation for this exercise: currentNe2 recognises chromosome-map information only when the input contains at least **two chromosomes**. A MAP containing only Scaffold_25 is read, but currentNe2 then reports `Number of chromosomes: Not given` and cannot use `-r` to estimate Ne from its physical positions. We can still demonstrate the program by supplying the approximate genetic span of Scaffold_25 explicitly, but currentNe2 will then assume that the markers are evenly distributed across that span. This is not a substitute for a multi-chromosome analysis.

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

## Step 2 — Calculate the teaching scaffold span

**Purpose:** approximate the genetic span covered by the Scaffold_25 markers under the same assumption used in the GONE2 lesson: 1 cM/Mb.

```bash
MAP="$OUTDIR/ABB_${CHROM}.map"

GENOME_MORGANS=$(awk '
  NR==1 {minimum=$4; maximum=$4}
  $4 < minimum {minimum=$4}
  $4 > maximum {maximum=$4}
  END {printf "%.8f", (maximum-minimum)/100000000}
' "$MAP")

printf 'Approximate analysed span: %s Morgans\n' "$GENOME_MORGANS"
```

The fourth MAP column is the physical base-pair position. At 1 cM/Mb, a span of 100 Mb corresponds to 100 cM, or 1 Morgan; therefore the base-pair span is divided by 100,000,000.

**Expected:** one positive decimal value. This is the span between the first and last marker, not the complete genetic length of the bear genome and not necessarily the full length of Scaffold_25.

## Step 3 — Estimate an illustrative contemporary size

**Purpose:** run currentNe2 using the supplied genetic span. Because currentNe2 does not accept a one-chromosome MAP as chromosome-map information, it will assume that the markers are evenly distributed over this total length.

```bash
CURRENTNE2_BIN="$PWD/software/currentNe2/currentne2"

for POPULATION in ABB SBB; do
  "$CURRENTNE2_BIN" -t 2 \
    "$OUTDIR/${POPULATION}_${CHROM}.ped" \
    "$GENOME_MORGANS"
done

ls -lh "$OUTDIR"/*_currentNe2_OUTPUT.txt
```

`-t 2` uses two threads. The final positional value, `GENOME_MORGANS`, supplies the analysed genetic span; it is not another option flag. We do not use `-r`, because physical marker locations are not used in this one-chromosome fallback. We also do not use `-x`, which requires chromosome assignments and fits a structured metapopulation model.

**Expected:** one output text file per population. currentNe2's documented default is to write beside the input with the `_currentNe2_OUTPUT.txt` suffix. In the output, `Genome size in Morgans` should now be the supplied positive value, and an Ne estimate should replace `Ne cannot be estimated because there is no map information`.

If the output still reports a genome size of `0.00`, check that `GENOME_MORGANS` is defined in the current terminal and appears after the PED filename in the command.

## Step 4 — Read and discuss

```bash
cat "$OUTDIR/ABB_${CHROM}_currentNe2_OUTPUT.txt"
cat "$OUTDIR/SBB_${CHROM}_currentNe2_OUTPUT.txt"
```

First check the preprocessing section:

- the total and effective numbers of individuals should be 10 for ABB and 8 for SBB;
- the input SNP count should match the PED, while the included count contains only polymorphic markers with less than 20% missing data;
- the proportion of missing data should be zero for these prepared files;
- `Genome size in Morgans` should be positive.

Then inspect observed d², expected and observed heterozygosity, the F statistic, inferred full-sibling pairs, and the reported contemporary Ne. A negative F, as observed in the first ABB attempt, means an excess of heterozygotes relative to Hardy–Weinberg expectations; it is a diagnostic that may reflect sampling, filtering, family composition, or biological processes and should not be silently ignored. The inferred sibling pairs are model-based candidates, not verified pedigrees.

Record the estimate and any uncertainty or diagnostic information. Ask whether ABB and SBB differ, whether the most recent part of GONE2 points in the same direction, and whether sampling noise or their shared LD assumptions could explain a mismatch.

**Do not present either value as the current effective size of the whole population.** This fallback ignores the actual spacing among markers and analyses only one scaffold. A substantive currentNe2 analysis should contain multiple autosomes so the program can use physical or genetic marker positions and compare within- and between-chromosome LD. It should also assess sample size, relatedness, filtering, and population structure.

