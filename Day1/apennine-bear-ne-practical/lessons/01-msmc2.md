# MSMC2 — historical effective population size

Estimated practical time: 30 minutes.

MSMC2 uses the spacing of heterozygous sites and the amount of callable sequence to infer coalescence rates through time.

## Predict

1. Is a region without heterozygous sites informative if it was not callable?
2. Which ends of a one-chromosome trajectory will be least reliable?
3. Should a chromosome result be noisier than a genome-wide result?

## Prepare the input

### Step 1 — Create an output directory

**Purpose:** keep derived files separate from immutable teaching inputs.

```bash
CHROM=chrN
VCF=data/teaching/single_bear.${CHROM}.vcf.gz
MASK=data/teaching/single_bear.${CHROM}.callable.bed.gz
OUTDIR=results/msmc2
mkdir -p "$OUTDIR"
```

**Expected:** no output. `mkdir -p` creates the directory if needed and does not complain if it already exists.

### Step 2 — Convert VCF plus mask to multihetsep

**Input:** one diploid VCF and its callable-region mask.

```bash
generate_multihetsep.py \
  --chr "$CHROM" \
  --mask "$MASK" \
  "$VCF" \
  > "$OUTDIR/single_bear.${CHROM}.multihetsep.txt"
```

**Expected:** the command writes a new text file and normally prints little or nothing to the terminal. A typical record has four fields:

```text
chrN    68306    44    TC
chrN    87563    259    AG
```

These values illustrate the format only. The fields represent chromosome, position, callable distance from the preceding segregating site, and observed alleles/haplotypes.

| Argument | Meaning | Decision to justify |
|---|---|---|
| `--chr` | Selects/labels the chromosome | Must match VCF labels exactly |
| `--mask` | Includes callable intervals | How was callability defined? |
| VCF path | Single diploid genotype input | Are depth and quality adequate? |
| `>` | Redirects output to a file | Is overwriting acceptable? |

Inspect it:

```bash
head "$OUTDIR/single_bear.${CHROM}.multihetsep.txt"
wc -l "$OUTDIR/single_bear.${CHROM}.multihetsep.txt"
```

**Check:** the file must be non-empty, positions should increase, chromosome labels should be consistent, and the allele field should not be missing. `wc -l` approximates the number of segregating records.

If the file is empty, check sample name, chromosome label, mask overlap, genotype filtering, and whether the VCF contains heterozygous variants.

### Question: one mask or two?

The official tool can intersect sample callability with a mappability mask:

```bash
generate_multihetsep.py \
  --chr "$CHROM" \
  --mask data/teaching/sample_callable.bed.gz \
  --mask data/teaching/mappable_regions.bed.gz \
  "$VCF" > "$OUTDIR/two_masks.multihetsep.txt"
```

Why is this preferable to treating every reference position as observable?

## Run MSMC2

### Step 3 — Examine available parameters

```bash
msmc2 --help | less
```

Press `q` to exit `less`. Locate `-t`, `-p`, and `-o` in the help before continuing.

### Step 4 — Fit a teaching model

**Input:** the multihetsep file created above.

```bash
msmc2 \
  -t 2 \
  -p '1*2+15*1+1*2' \
  -o "$OUTDIR/bear_chr" \
  "$OUTDIR/single_bear.${CHROM}.multihetsep.txt"
```

**Expected:** MSMC2 reports optimization progress and creates files beginning `bear_chr`, including a final table. Confirm:

```bash
ls -lh "$OUTDIR"/bear_chr*
head "$OUTDIR/bear_chr.final.txt"
```

Expected table shape:

```text
time_index  left_time_boundary  right_time_boundary  lambda_00
0           ...                 ...                  ...
```

**Check:** time boundaries and rate values should be numeric and the final file should contain multiple rows. A completed run is not automatically a trustworthy run.

| Flag | Meaning | Question |
|---|---|---|
| `-t 2` | Two compute threads | How many cores are allocated? |
| `-p` | Groups time segments into free parameters | Is the dataset rich enough for more parameters? |
| `-o` | Output prefix | Can runs with different settings be distinguished? |

The teaching pattern combines two intervals at each end and estimates 15 middle intervals separately. More parameters do not automatically improve inference; sparse data can make a flexible curve unstable.

If necessary, stop with `Ctrl-C` and use the precomputed checkpoint:

```bash
cp results/msmc2/precomputed/msmc2.final.txt "$OUTDIR/bear_chr.final.txt"
```

## Scale and plot in R

### Step 5 — Convert scaled units to biological units

**Input:** the MSMC2 final table, mutation rate per site per generation, and generation time in years.

```bash
R
```

```r
x <- read.table("results/msmc2/bear_chr.final.txt", header=TRUE)
names(x)
mu <- 4.5e-9
generation_time <- 10
midpoint <- sqrt(x$left_time_boundary * x$right_time_boundary)
years <- midpoint / mu * generation_time
lambda <- if ("lambda_00" %in% names(x)) x$lambda_00 else x$lambda
Ne <- 1 / (2 * mu * lambda)
plot(years, Ne, type="s", log="xy",
     xlab="Years before present", ylab="Effective population size")
```

**Expected:** `names(x)` lists the time-boundary and lambda columns. The plot has logarithmic time and population-size axes and a step-like trajectory.

**Check:** `years` and `Ne` should be positive and finite:

```r
summary(years)
summary(Ne)
```

`NA`, `Inf`, or non-positive values indicate a parsing, column-selection, or numerical problem that must be resolved before interpretation.

The values of `mu` and `generation_time` are placeholders for teaching, not recommended bear parameters.

### Parameter challenge

Repeat with alternative values and explain:

- Is `mu` measured per site per generation?
- Is it estimated for bears or borrowed from another mammal?
- Does generation time represent this population?
- Which axes change with `mu`?
- Which axis changes with generation time?

Exit R using `q()`. Continue to [GONE](02-gone.md).

