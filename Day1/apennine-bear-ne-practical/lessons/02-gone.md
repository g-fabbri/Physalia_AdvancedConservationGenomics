# GONE — recent demographic history from LD

Estimated practical time: 30 minutes.

GONE uses linkage disequilibrium at different recombination distances to infer effective population size over recent generations. Drift creates LD; recombination removes it. The method therefore requires genotypes from **multiple individuals**, not the single MSMC2 genome.

## What students need to understand

The original Apennine bear driver script performs many bookkeeping operations internally. We will not reproduce all of them. Our practical focuses on decisions that can change the biological result:

```text
population sample and variant filtering
                 ↓
chromosome and recombination information
                 ↓
LD-distance and replicate parameters
                 ↓
GONE estimate and diagnostic checks
```

The driver will handle temporary control files, splitting, LD-bin aggregation, and output renaming. Those are useful for software development, but not central learning objectives here.

## Predict before running

Which factors can generate LD besides small population size? Consider relatives, migrants, pooled populations, selection, physical linkage, and genotype errors.

## Step 1 — Define the dataset

Replace `chrN` with the chromosome selected by the instructor:

```bash
CHROM=chrN
VCF=data/teaching/population.${CHROM}.vcf.gz
OUTDIR=results/gone
PREFIX="$OUTDIR/population_${CHROM}"
mkdir -p "$OUTDIR"
```

Check the input:

```bash
ls -lh "$VCF"
bcftools query -l "$VCF" | wc -l
bcftools view -m2 -M2 -v snps "$VCF" -Ou | bcftools view -H | wc -l
```

**Expected:** a non-empty VCF, the number of sampled individuals, and the number of biallelic SNPs. Record both counts.

### Question

Is the number of individuals sufficient to estimate LD precisely? More SNPs cannot completely compensate for very few independently sampled individuals.

## Step 2 — Filter variants and create PED/MAP input

```bash
plink \
  --vcf "$VCF" \
  --double-id \
  --allow-extra-chr \
  --snps-only just-acgt \
  --biallelic-only strict \
  --geno 0.10 \
  --mac 2 \
  --recode \
  --out "${PREFIX}_original_label"
```

| Flag | Effect | Decision question |
|---|---|---|
| `--double-id` | Uses sample ID as family and individual ID | Are relatives known and removed? |
| `--allow-extra-chr` | Accepts non-human chromosome labels | Are only autosomal markers included? |
| `--snps-only just-acgt` | Keeps canonical SNP alleles | Why exclude indels? |
| `--biallelic-only strict` | Keeps exactly two alleles | Is this appropriate for GONE? |
| `--geno 0.10` | Removes loci with >10% missing genotypes | Is 10% reasonable for this sample size? |
| `--mac 2` | Requires at least two minor-allele copies | Would three copies be more defensible? |
| `--recode` | Writes PED/MAP files | Does the GONE version expect these? |

Inspect the output:

```bash
tail -20 "${PREFIX}_original_label.log"
wc -l "${PREFIX}_original_label.ped" "${PREFIX}_original_label.map"
head "${PREFIX}_original_label.map"
```

**Expected:** PED rows equal sampled individuals and MAP rows equal retained loci. PLINK should finish without `Error:`.

### Why use MAC instead of MAF here?

For 10 diploid individuals, `--maf 0.02` represents less than one chromosome copy and removes almost nothing. `--mac 2` expresses the desired minimum allele count directly. Calculate what `--mac 2` and `--mac 3` mean as sample frequencies in your dataset.

## Step 3 — Create a one-chromosome GONE MAP

The supplied Apennine script assumes chromosomes are numbered consecutively from 1 and obtains the chromosome count from the final MAP row. Therefore, a selected chromosome labelled `12` would incorrectly make it expect 12 chromosomes.

Copy the PED and recode only column 1 of the derived MAP:

```bash
cp "${PREFIX}_original_label.ped" "${PREFIX}.ped"
awk 'BEGIN {OFS="\t"} {$1=1; print $1,$2,$3,$4}' \
  "${PREFIX}_original_label.map" > "${PREFIX}.map"
```

Validate:

```bash
cut -f1 "${PREFIX}.map" | sort -u
head "${PREFIX}.map"
wc -l "${PREFIX}.ped" "${PREFIX}.map"
```

**Expected:** the chromosome check prints only `1`. Physical positions and locus count remain unchanged. Record the original biological chromosome name separately.

## Step 4 — Choose the parameters

The supplied Apennine bear analysis used:

| Parameter | Original value | Biological or computational meaning |
|---|---:|---|
| `PHASE` | `2` | Genotypes have unknown phase |
| `cMMb` | `1` | Assumed cM/Mb if MAP genetic distances are zero |
| `DIST` | `1` | Haldane mapping correction |
| `NGEN` | `2000` | Generations represented by LD bins |
| `NBIN` | `400` | 400 bins: five generations per bin |
| `MAF` | `0.0` | No additional internal frequency filter |
| `ZERO` | `1` | Allow missing genotype codes |
| `maxNCHROM` | `-99` | Analyze all detected chromosomes—one here |
| `maxNSNP` | `50000` | Approximate SNP cap per chromosome |
| `hc` | `0.01` | Maximum recombination fraction used |
| `REPS` | `40` | Replicate estimates |
| `threads` | `10` | Parallel workers |

Discuss before editing:

1. If genetic positions in the MAP are zero, what evidence supports `cMMb=1`?
2. Why might `hc=0.01` have been chosen instead of the general `0.05` recommendation?
3. Does `NGEN=2000` mean the data truly resolve 2,000 generations?
4. Which parameters affect the inference, and which primarily affect runtime?

For the live run, keep the biological settings and change only:

- `REPS=5` to shorten computation;
- `threads` to the number of cores allocated by the instructor.

Compare the live output with a precomputed run using the original `REPS=40`.

## Step 5 — Create a private run directory

The original driver uses generic temporary filenames and removes previous results. Each group therefore needs a private copy of GONE.

```bash
COURSE_DIR=$(pwd)
DATA_DIR="$COURSE_DIR/results/gone"
GONE_SOURCE="$COURSE_DIR/software/GONE"
GROUP=group01
RUN_DIR="$COURSE_DIR/work/gone_${CHROM}_${GROUP}"
FILE=population_chrN
mkdir -p "$COURSE_DIR/work"
cp -R "$GONE_SOURCE" "$RUN_DIR"
```

Replace `group01` with your group and `population_chrN` with the basename used in Step 3. `RUN_DIR` must be new; do not reuse another group's directory.

Enter it and inspect the configuration:

```bash
cd "$RUN_DIR"
less INPUT_PARAMETERS_FILE
```

Press `q`, then edit `REPS` and `threads`:

```bash
nano INPUT_PARAMETERS_FILE
```

Save using `Ctrl-O`, Enter, then exit with `Ctrl-X`. Confirm:

```bash
grep -E '^(PHASE|cMMb|DIST|NGEN|NBIN|MAF|ZERO|maxNCHROM|maxNSNP|hc|REPS|threads)=' \
  INPUT_PARAMETERS_FILE
```

## Step 6 — Run GONE

```bash
bash script_GONE.sh "$FILE" "$DATA_DIR"
```

The first argument is the PED/MAP basename without an extension. The second is the directory containing the inputs and receiving results.

Expected progress:

```text
DIVIDE .ped AND .map FILES IN CHROMOSOMES
RUNNING ANALYSIS OF CHROMOSOMES ...
CHROMOSOME ANALYSES took ... seconds
Running GONE
GONE run took ... seconds
END OF ANALYSES
```

Internally, the driver runs four important components:

| Component | Role |
|---|---|
| `MANAGE_CHROMOSOMES2` | Prepares and optionally subsamples chromosome data |
| `LD_SNP_REAL3` | Calculates LD in recombination-distance bins |
| `SUMM_REP_CHROM3` | Standardizes and combines LD summaries |
| `GONEparallel.sh` | Performs replicate demographic inference |

Students do not need to manage the temporary files produced between these programs.

## Step 7 — Validate outputs

Return to the repository:

```bash
cd "$COURSE_DIR"
ls -lh "$DATA_DIR"
```

Expected principal outputs:

```text
OUTPUT_population_chrN
Output_Ne_population_chrN
Output_d2_population_chrN
outfileHWD
timefile
seedfile
TEMPORARY_FILES/
```

Inspect them:

```bash
head "$DATA_DIR/Output_Ne_${FILE}"
cat "$DATA_DIR/timefile"
cat "$DATA_DIR/outfileHWD"
```

**Check:** the `Ne` file is non-empty and contains multiple generations; `timefile` shows completed stages; `outfileHWD` contains the Hardy–Weinberg diagnostic.

## Step 8 — Plot the trajectory

Start R:

```bash
R
```

Read the result first, without assuming its columns:

```r
gone <- read.table("results/gone/Output_Ne_population_chrN", header=TRUE)
names(gone)
head(gone)
```

If the first row of the file contains numbers rather than column names, read it with `header=FALSE` and assign names after inspecting it:

```r
gone <- read.table("results/gone/Output_Ne_population_chrN", header=FALSE)
names(gone)[1:2] <- c("Generation", "Ne")
```

Identify the columns representing generation and `Ne`, then plot them. Replace the placeholder column names below with the names printed by `names(gone)`:

```r
plot(gone$Generation, gone$Ne,
     type="l", log="y", lwd=2,
     xlab="Generations before present",
     ylab="Effective population size")
```

If the most recent generation should appear on the right, reverse the x-axis:

```r
plot(gone$Generation, gone$Ne,
     type="l", log="y", lwd=2,
     xlim=rev(range(gone$Generation)),
     xlab="Generations before present",
     ylab="Effective population size")
```

Do not copy the example column names blindly. Inspect the real output first.

## Interpret before moving on

1. Could structure, relatives, or migration produce the recent pattern?
2. Does the Hardy–Weinberg diagnostic support panmixia?
3. How variable is the five-replicate classroom result compared with 40 replicates?
4. What information is lost by using one chromosome?
5. Which portion of the trajectory would you avoid interpreting?

Continue to [NeEstimator](03-neestimator.md).


