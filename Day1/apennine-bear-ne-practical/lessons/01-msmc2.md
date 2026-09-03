# MSMC2 — historical effective population size

Estimated practical time: 30 minutes.

MSMC2 uses the spacing of heterozygous sites and the amount of callable sequence to infer coalescence rates through time.

## Predict

1. Is a region without heterozygous sites informative if it was not callable?
2. Which ends of a one-chromosome trajectory will be least reliable?
3. Should a chromosome result be noisier than a genome-wide result?

## Prepare the input

```bash
CHROM=chrN
VCF=data/teaching/single_bear.${CHROM}.vcf.gz
MASK=data/teaching/single_bear.${CHROM}.callable.bed.gz
OUTDIR=results/msmc2
mkdir -p "$OUTDIR"
```

```bash
generate_multihetsep.py \
  --chr "$CHROM" \
  --mask "$MASK" \
  "$VCF" \
  > "$OUTDIR/single_bear.${CHROM}.multihetsep.txt"
```

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

```bash
msmc2 --help | less
```

```bash
msmc2 \
  -t 2 \
  -p '1*2+15*1+1*2' \
  -o "$OUTDIR/bear_chr" \
  "$OUTDIR/single_bear.${CHROM}.multihetsep.txt"
```

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

The values of `mu` and `generation_time` are placeholders for teaching, not recommended bear parameters.

### Parameter challenge

Repeat with alternative values and explain:

- Is `mu` measured per site per generation?
- Is it estimated for bears or borrowed from another mammal?
- Does generation time represent this population?
- Which axes change with `mu`?
- Which axis changes with generation time?

Exit R using `q()`. Continue to [GONE](02-gone.md).

