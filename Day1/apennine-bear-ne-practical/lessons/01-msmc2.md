# MSMC2 — comparing historical population-size trajectories

Estimated practical time: 30 minutes.

We begin with one representative from each population:

| Population | Individual |
|---|---|
| Apennine | 4573 |
| Slovak | U1916 |

Our question is:

> Do these two genomes record different histories of coalescence and effective population size?

MSMC2 does not read BAM or VCF files directly. It reads **multihetsep**, a format combining segregating-site information with the amount of callable sequence between sites.

~~~text
BAM + reference
      ↓ variant calling and callability assessment
VCF + callable mask
      ↓ generate_multihetsep.py
multihetsep
      ↓ MSMC2
historical trajectory
~~~

## A note about the common mask

Both bears were aligned to the Apennine reference. We use **UrArMa_callable.bed.gz** as the common teaching mask. This is defensible only if it represents reference mappability or regions callable in both samples. Alignment to the same reference alone does not make a sample-specific depth mask transferable.

We will additionally identify and exclude missing genotypes separately for each bear.

## Step 1 — Define files and individuals

**Purpose:** name the inputs and representatives.

**Input:** two prepared single-sample VCFs, the common callable mask, and Scaffold_34.

~~~bash
CHROM=Scaffold_34
APN_RAW=data/UrArMa_4573_s34.vcf.gz
SVK_RAW=data/UrArMa_U1916_s34.vcf.gz
MASK=data/UrArMa_callable.bed.gz
APN_ID=4573
SVK_ID=U1916
OUTDIR=results/msmc2
mkdir -p "$OUTDIR"
~~~

**Expected:** no terminal output.

**Check:**

~~~bash
ls -lh "$APN_RAW" "$SVK_RAW" "$MASK"
bcftools query -l "$APN_RAW" | grep -x "$APN_ID"
bcftools query -l "$SVK_RAW" | grep -x "$SVK_ID"
~~~

Both identifiers should be printed.

## Step 2 — Verify the prepared representative genomes

**Purpose:** confirm that each prepared VCF contains only its intended representative.

**Input:** `UrArMa_4573_s34.vcf.gz` and `UrArMa_U1916_s34.vcf.gz`.

~~~bash
bcftools query -l "$APN_RAW"
bcftools query -l "$SVK_RAW"
bcftools index -n "$APN_RAW"
bcftools index -n "$SVK_RAW"
~~~

**Expected:** sample IDs `4573` and `U1916`, followed by one record count per VCF.

**Check:**

~~~bash
bcftools query -l "$APN_RAW" | wc -l
bcftools query -l "$SVK_RAW" | wc -l
~~~

Each command should print `1`.


## Step 3 — Create the multihetsep files

**Purpose:** combine allele and callability information in MSMC2 format.

**Input:** one complete VCF, the common positive mask, and the individual's missing-site negative mask.

~~~bash
generate_multihetsep.py \
  --chr "$CHROM" \
  --mask "$MASK" \
  --negative_mask "$OUTDIR/${APN_ID}.missing.bed.gz" \
  "$OUTDIR/${APN_ID}.complete.vcf.gz" \
  > "$OUTDIR/${APN_ID}.${CHROM}.multihetsep.txt"

generate_multihetsep.py \
  --chr "$CHROM" \
  --mask "$MASK" \
  --negative_mask "$OUTDIR/${SVK_ID}.missing.bed.gz" \
  "$OUTDIR/${SVK_ID}.complete.vcf.gz" \
  > "$OUTDIR/${SVK_ID}.${CHROM}.multihetsep.txt"
~~~

**Expected:** the program reports **2 haplotypes** for each diploid genome and creates two non-empty files.

Typical format:

~~~text
Scaffold_34    68306    44     TC
Scaffold_34    87563    259    AG
~~~

Columns represent chromosome, segregating-site position, callable sites since the previous segregating site, and observed haplotype alleles.

**Check:**

~~~bash
head "$OUTDIR/${APN_ID}.${CHROM}.multihetsep.txt"
head "$OUTDIR/${SVK_ID}.${CHROM}.multihetsep.txt"
wc -l "$OUTDIR/"*.multihetsep.txt
~~~

Positions should increase, chromosome labels should match, and neither file should be empty.

## Step 4 — Run MSMC2 separately

**Purpose:** estimate a historical trajectory for each diploid genome.

~~~bash
msmc2 -t 2 -p '1*2+15*1+1*2' \
  -o "$OUTDIR/APN_4573" \
  "$OUTDIR/${APN_ID}.${CHROM}.multihetsep.txt"

msmc2 -t 2 -p '1*2+15*1+1*2' \
  -o "$OUTDIR/SVK_U1916" \
  "$OUTDIR/${SVK_ID}.${CHROM}.multihetsep.txt"
~~~

| Flag | Meaning |
|---|---|
| **-t 2** | Use two compute threads |
| **-p** | Group atomic time intervals into estimated parameters |
| **-o** | Set the population-specific output prefix |

**Expected:** two final tables:

~~~text
results/msmc2/APN_4573.final.txt
results/msmc2/SVK_U1916.final.txt
~~~

**Check:**

~~~bash
head "$OUTDIR/APN_4573.final.txt"
head "$OUTDIR/SVK_U1916.final.txt"
~~~

Both should contain time boundaries and a lambda column.

## Step 5 — Scale and compare the trajectories

**Purpose:** apply the same mutation rate and generation time so the curves are directly comparable.

Start R:

~~~bash
R
~~~

~~~r
apn <- read.table("results/msmc2/APN_4573.final.txt", header=TRUE)
svk <- read.table("results/msmc2/SVK_U1916.final.txt", header=TRUE)

mu <- 4.5e-9
generation_time <- 10

scale_msmc <- function(x) {
  lambda <- if ("lambda_00" %in% names(x)) x$lambda_00 else x$lambda
  midpoint <- sqrt(x$left_time_boundary * x$right_time_boundary)
  data.frame(
    years = midpoint / mu * generation_time,
    Ne = 1 / (2 * mu * lambda)
  )
}

apn_scaled <- scale_msmc(apn)
svk_scaled <- scale_msmc(svk)

plot(apn_scaled$years, apn_scaled$Ne,
     type="s", log="xy", lwd=2, col="firebrick",
     xlab="Years before present", ylab="Effective population size")
lines(svk_scaled$years, svk_scaled$Ne,
      type="s", lwd=2, col="steelblue")
legend("topleft",
       legend=c("Apennine 4573", "Slovak U1916"),
       col=c("firebrick", "steelblue"), lwd=2)
~~~

**Expected:** two trajectories on identical logarithmic axes.

**Check:**

~~~r
summary(apn_scaled)
summary(svk_scaled)
~~~

All plotted values should be positive and finite. The values of **mu** and generation time are teaching placeholders and must be justified for biological use.

Exit R with **q()**.

## Questions for discussion

1. Where do the two curves differ, and where do they overlap?
2. Which parts of each curve are least reliable?
3. Could unequal coverage or callability imitate a population difference?
4. What does the third multihetsep column contribute that a variant-only VCF does not?
5. How would alternative mutation rates and generation times change the figure?
6. Why should a one-scaffold trajectory not be treated as a final demographic reconstruction?

The first piece of evidence is now complete. Continue to [recent population history with GONE](02-gone.md).
