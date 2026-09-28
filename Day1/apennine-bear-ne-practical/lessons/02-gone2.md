# Part 2 — GONE2: comparing recent history in two bear populations


## Start your terminal

From the Day 1 course directory, activate the course environment and identify the local GONE2 executable:

~~~bash
conda activate bear-ne-practical
COURSE_DIR=$(pwd)
GONE2_BIN="$COURSE_DIR/software/GONE2/gone2"
~~~

`conda activate` makes BCFtools, PLINK, and the other shared dependencies available. GONE2 is compiled separately inside the course `software/` directory, so `GONE2_BIN` records its complete path. Run this block whenever you open a new terminal.

**Check:**

~~~bash
command -v bcftools
command -v plink
test -x "$GONE2_BIN" && echo "GONE2 ready: $GONE2_BIN"
~~~

The first two commands should print paths inside the `bear-ne-practical` environment. The last command should print the complete path to the executable. If it prints nothing, confirm that you are in the course root and that GONE2 was compiled before class.

> GONE2 estimates recent effective population size from linkage disequilibrium (LD) at different recombination distances. We analyse Apennine brown bears (ABB) and Slovak brown bears (SBB) separately and compare their trajectories.

This is a one-scaffold teaching analysis. The small samples and limited genomic coverage mean that diagnostics are as important as the estimated curves.

~~~text
18-individual VCF
        ↓ separate ABB and SBB
population VCFs
        ↓ filter and convert with PLINK
PED + MAP
        ↓ GONE2
recent Ne trajectories and diagnostics
~~~

Software is prepared before class. Installation instructions are in the [software README](../software/README.md).

## Step 1 — Prepare ABB and SBB VCFs

**Purpose:** separate the two populations using the sample lists created during shared QC.

**Input:** the checked 18-individual VCF and the ABB/SBB sample lists.

~~~bash
CHROM=Scaffold_25
ALL_VCF=data/UrArMa_18i_s25.vcf.gz
ABB_VCF=data/ABB_s25.vcf.gz
SBB_VCF=data/SBB_s25.vcf.gz
~~~

~~~bash
bcftools view -S data/apennine.samples "$ALL_VCF" \
  -Oz -o "$ABB_VCF"
bcftools view -S data/slovak.samples "$ALL_VCF" \
  -Oz -o "$SBB_VCF"
~~~

~~~bash
bcftools index -t "$ABB_VCF"
bcftools index -t "$SBB_VCF"
~~~

**Expected:** **ABB_s25.vcf.gz** contains the 10 ABB individuals and **SBB_s25.vcf.gz** contains the 8 SBB individuals.


## Step 2 — Create the PED/MAP files

**Purpose:** convert each population VCF into the PED/MAP format required by GONE2, then adjust the chromosome labels.

PLINK's text format uses two files with the same prefix:

- The **PED file** contains one row per individual. The first six columns describe the sample and family, followed by two allele entries for every marker.
- The **MAP file** contains one row per marker: chromosome, marker ID, genetic-map position, and physical base-pair position.

The marker order in the MAP file must exactly match the genotype order in every PED row.

### 2.1 — Define the input directory

```bash
INPUTDIR=results/gone2/input
mkdir -p "$INPUTDIR"
```

### 2.2 — Convert the VCF files with PLINK

**Purpose:** retain canonical, biallelic SNPs and convert each population VCF into PED/MAP files.

```bash
for POPULATION in ABB SBB; do
  VCF="data/${POPULATION}_s25.vcf.gz"
  PREFIX="$INPUTDIR/${POPULATION}_${CHROM}"

  plink \
    --vcf "$VCF" \
    --double-id \
    --allow-extra-chr \
    --snps-only just-acgt \
    --biallelic-only strict \
    --recode \
    --out "${PREFIX}_original_label"
done
```

The loop applies exactly the same conversion to ABB and SBB.

- `--double-id` uses the VCF sample ID as both the family and individual ID.
- `--allow-extra-chr` permits scaffold names such as `Scaffold_25`.
- `--snps-only just-acgt` retains SNPs containing only A, C, G or T.
- `--biallelic-only strict` retains sites with exactly two alleles.
- `--recode` writes text-format PED and MAP files.

The temporary output files retain the original scaffold label:

```text
ABB_Scaffold_25_original_label.ped
ABB_Scaffold_25_original_label.map
SBB_Scaffold_25_original_label.ped
SBB_Scaffold_25_original_label.map
```

### 2.3 — Prepare the files for GONE2

**Purpose:** preserve the PED genotypes while replacing the scaffold name in the MAP file with chromosome code `1`, which GONE2 accepts.

```bash
for POPULATION in ABB SBB; do
  PREFIX="$INPUTDIR/${POPULATION}_${CHROM}"

  cp "${PREFIX}_original_label.ped" "${PREFIX}.ped"

  awk 'BEGIN {OFS="\t"} {
    $1=1
    print $1,$2,$3,$4
  }' "${PREFIX}_original_label.map" > "${PREFIX}.map"
done
```

The `cp` command creates the final PED file without modifying its genotypes.

The `awk` command processes the MAP file:

- `$1=1` changes the first column from `Scaffold_25` to chromosome code `1`.
- `$2`, `$3`, and `$4` preserve the marker ID, genetic-map position and physical position.
- `OFS="\t"` writes a tab-separated MAP file.

Only the chromosome label is changed. Marker IDs, marker order, genetic positions and physical positions remain unchanged.

### 2.4 — Check the final files

```bash
wc -l \
  "$INPUTDIR/ABB_${CHROM}.ped" \
  "$INPUTDIR/SBB_${CHROM}.ped"
```

```bash
cut -f1 "$INPUTDIR/ABB_${CHROM}.map" | sort -u
cut -f1 "$INPUTDIR/SBB_${CHROM}.map" | sort -u
```

The PED files should contain 10 ABB individuals and 8 SBB individuals. Both MAP checks should print only:

```text
1
```


## Step 3 — Run GONE2

**Purpose:** estimate a recent Ne trajectory for each population under identical settings.

~~~bash
COURSE_DIR=$(pwd)
GONE2_BIN="$COURSE_DIR/software/GONE2/gone2"
OUTDIR="$COURSE_DIR/results/gone2"
~~~

~~~bash
"$GONE2_BIN" \
  -g 0 \
  -r 1 \
  -u 0.02 \
  -t 2 \
  -S 1 \
  -E \
  -o "$OUTDIR/ABB" \
  "$COURSE_DIR/$INPUTDIR/ABB_${CHROM}.ped"
~~~

~~~bash
"$GONE2_BIN" \
  -g 0 \
  -r 1 \
  -u 0.02 \
  -t 2 \
  -S 1 \
  -E \
  -o "$OUTDIR/SBB" \
  "$COURSE_DIR/$INPUTDIR/SBB_${CHROM}.ped"
~~~

| Option | Meaning |
|---|---|
| **-g 0** | Treat genotypes as unphased diploids |
| **-r 1** | Assume a constant recombination rate of 1 cM/Mb |
| **-u 0.02** | Use LD bins only up to a recombination fraction of 0.02 |
| **-t 2** | Use two threads |
| **-S 1** | Fix the random seed for reproducibility |
| **-E** | Request variation among genetic-algorithm rounds |
| **-o** | Set the output prefix |

### Why use `-u 0.02`?

GONE2 relates LD between marker pairs to the time in the past that generated that LD. As a useful approximation, LD at recombination fraction `c` is most informative about approximately `1/(2c)` generations ago. The default upper bound is `0.05`, corresponding roughly to information from about 10 generations ago. Here we set the upper bound to `0.02`, corresponding roughly to **25 generations ago**.

This choice excludes the most weakly linked marker pairs and focuses the fit on more tightly linked pairs. For this one-scaffold exercise, it reduces sensitivity to noisy LD at larger distances and to error from assuming a uniform recombination rate of 1 cM/Mb instead of using a detailed genetic map. The trade-off is important: `-u 0.02` provides **less information about the most recent generations** than the default `0.05`; it is not universally a better value. It should be reported explicitly, and a formal analysis should compare plausible `-u` settings and, where possible, use an empirical recombination map.

Do not use **-x** in the main run. It fits a structured metapopulation model, while this exercise initially treats ABB and SBB as separate populations.

**Expected:** progress messages for reading the data, measuring d², estimating Ne, and writing output files.

## Step 4 — Inspect the results

**Purpose:** understand the structure of the output before drawing or interpreting the demographic trajectories.

GONE2 writes several files for each population. The two most important for this practical are **GONE2_STATS**, which documents the run and its diagnostics, and **GONE2_Ne**, which contains the estimated trajectory.

### Step 4.1 — Read the run summary and diagnostics

**Input:** the two `GONE2_STATS` files.

~~~bash
cat "$OUTDIR/ABB_GONE2_STATS"
~~~

~~~bash
cat "$OUTDIR/SBB_GONE2_STATS"
~~~

The exact layout can vary slightly among GONE2 versions, but the file records information such as:

| Item | Meaning and interpretation |
|---|---|
| **Individuals and loci** | Number of samples and markers actually read. Unexpected values may indicate a problem in the PED/MAP preparation. |
| **Genome or chromosome length** | Genetic length inferred from the MAP positions and the recombination-rate assumption. An unrealistic value changes how LD distance is translated into time. |
| **Fis** | Departure from Hardy–Weinberg genotype proportions. A high positive value can reflect inbreeding, population structure, related individuals, or genotype/missing-data problems. |
| **LD-fit information** | Describes how well the fitted model reproduces the observed decay of LD. A poor fit suggests that the inferred trajectory does not adequately explain the data. |
| **Selected model/combination** | The optimisation solution retained by GONE2. Different retained solutions across sensitivity runs can indicate instability. |
| **Warnings** | Messages about structure, excessive Fis, chromosome length, marker limits, or model failure. These must be investigated before biological interpretation. |

**Expected:** both files should report a completed analysis rather than an error. Confirm that ABB has 10 individuals, SBB has 8, and that each run contains a plausible nonzero number of loci.

**Check:** if GONE2 produces only a `GONE2_STATS` file and no `GONE2_Ne` file, inspect the end of the summary for the reason:

~~~bash
ls -lhrt "$OUTDIR"
~~~


### Step 4.2 — Inspect the Ne trajectories

**Input:** the two `GONE2_Ne` tables.

~~~bash
head "$OUTDIR/ABB_GONE2_Ne"
~~~

~~~bash
head "$OUTDIR/SBB_GONE2_Ne"
~~~

Each row represents one point in the reconstructed demographic trajectory:

| Column | Meaning |
|---|---|
| **Generation** | Time before the sampled generation. Small values are more recent; larger values are further in the past. |
| **Ne** or a column beginning with **Ne** | Estimated effective population size at that generation. Ne is the size of an idealised population experiencing the observed genetic drift, not a direct census count. |

Depending on the GONE2 version and options, the table may contain additional Ne-related columns. Display the column names with:

~~~bash
head -n 1 "$OUTDIR/ABB_GONE2_Ne"
head -n 1 "$OUTDIR/SBB_GONE2_Ne"
~~~

**Expected:** generation values should increase into the past, while Ne should be positive. Large jumps between adjacent generations—especially near the youngest or oldest boundary—should be treated cautiously because resolution is not uniform through time.

**Check:** confirm that both trajectory files exist, are nonempty, and contain more than a header:

~~~bash
wc -l "$OUTDIR/ABB_GONE2_Ne" "$OUTDIR/SBB_GONE2_Ne"
~~~

### Step 4.3 — Recognise the supporting LD file

Principal outputs:

| File ending | Content |
|---|---|
| **GONE2_Ne** | Generation and estimated Ne |
| **GONE2_d2** | Observed and predicted LD by recombination bin |
| **GONE2_STATS** | Inputs, parameters, Hardy–Weinberg diagnostics, runtime, and warnings |

The `GONE2_d2` file contains the observed LD statistic and the values predicted by the fitted demographic model across recombination-distance bins. Students do not need to plot it during this short exercise, but it is useful for checking whether a visually attractive Ne trajectory is actually supported by a reasonable fit to the LD data.

Before interpreting either population, examine the number of individuals and loci, Fis, inferred genetic length, fit information, and every warning in `GONE2_STATS`. The `GONE2_Ne` file is the model result; the `GONE2_STATS` and `GONE2_d2` files help determine whether that result is trustworthy.

## Step 5 — Plot ABB and SBB

**Purpose:** compare ABB and SBB over a chosen number of generations and save the figure directly as a PDF.

The provided [R script](../scripts/plot_gone2.R) reads both `GONE2_Ne` tables, finds the Ne column, removes non-positive or non-finite values, restricts the trajectories to the requested generations, and draws them with a logarithmic Ne axis. The first command-line number is the youngest generation to display and the second is the oldest.

For example, plot generations 1–100 from the course root:

~~~bash
Rscript scripts/plot_gone2.R 1 100
~~~

To choose another interval, replace `1 100`; for example, `10 75` plots generations 10–75. With no numbers, the script defaults to generations 1–100.

**Expected:** the script prints the selected interval and creates **results/gone2/GONE2_ABB_SBB_generations_1_100.pdf**. The present is toward the left and older generations are toward the right.

**Check:**

~~~bash
ls -lh results/gone2/GONE2_ABB_SBB_generations_1_100.pdf
~~~

The PDF should have a nonzero size and contain both coloured trajectories. Do not request generations outside those available in the two result tables.


