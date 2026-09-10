# MSMC2 — comparing historical population-size trajectories

Estimated practical time: 30 minutes.

We compare one representative genome from each population:

| Population | Meaning | Individual |
|---|---|---|
| ABB | Apennine brown bear | 4573 |
| SBB | Slovak brown bear | U1916 |

Our biological question is:

> Do the ABB and SBB genomes record different histories of coalescence and effective population size?

MSMC2 does not read BAM or VCF files directly. It reads **multihetsep**, a format that combines segregating sites with the amount of callable sequence between them.

~~~text
already-filtered VCF + callable mask
               ↓ generate_multihetsep.py
           multihetsep
               ↓ MSMC2
     historical Ne trajectory
~~~

The preliminary QC has already confirmed the files, sample IDs, chromosome labels, indexes, and genotype filtering. We therefore begin by preparing the MSMC2 input rather than repeating those checks.

## A note about the common mask

Both bears were aligned to the Apennine reference, and this exercise uses **UrArMa_callable.bed.gz** as a common teaching mask. This is appropriate if the file describes reference mappability or regions callable in both genomes. In a full analysis, sample-specific callability masks—or their intersection for a comparison—are preferable when callability was estimated from read depth and genotype quality.

## Step 1 — Define the analysis variables

**Purpose:** assign short, meaningful names to the prepared inputs and create an output directory.

**Input:** the two already-filtered single-sample VCFs, the common callable mask, and Scaffold_34.

~~~bash
CHROM=Scaffold_34
ABB_VCF=data/UrArMa_4573_s34.vcf.gz
SBB_VCF=data/UrArMa_U1916_s34.vcf.gz
ABB_ID=4573
SBB_ID=U1916
MASK=data/UrArMa_callable.bed.gz
OUTDIR=results/msmc2

mkdir -p "$OUTDIR"
~~~

**Expected:** nothing is printed. The commands assign shell variables and create **results/msmc2** if necessary.

**Check:**

~~~bash
printf 'ABB: %s (%s)\nSBB: %s (%s)\nMask: %s\n' \
  "$ABB_ID" "$ABB_VCF" "$SBB_ID" "$SBB_VCF" "$MASK"
~~~

Expected shape:

~~~text
ABB: 4573 (data/UrArMa_4573_s34.vcf.gz)
SBB: U1916 (data/UrArMa_U1916_s34.vcf.gz)
Mask: data/UrArMa_callable.bed.gz
~~~

## Step 2 — Create the multihetsep files

**Purpose:** combine genotype information and callable sequence in the format required by MSMC2.

**Input:** one already-filtered single-sample VCF and the common positive mask for each run.

~~~bash
generate_multihetsep.py \
  --chr "$CHROM" \
  --mask "$MASK" \
  "$ABB_VCF" \
  > "$OUTDIR/ABB_${ABB_ID}.${CHROM}.multihetsep.txt"

generate_multihetsep.py \
  --chr "$CHROM" \
  --mask "$MASK" \
  "$SBB_VCF" \
  > "$OUTDIR/SBB_${SBB_ID}.${CHROM}.multihetsep.txt"
~~~

| Argument | Meaning |
|---|---|
| **--chr "$CHROM"** | Process only Scaffold_34 |
| **--mask "$MASK"** | Count only positions included in the callable mask |
| **"$ABB_VCF" / "$SBB_VCF"** | Read the prepared diploid genome for that population |
| **> output file** | Save the generated multihetsep text |

**Expected:** each run reports that it is generating input for **2 haplotypes**, because each VCF contains one diploid individual. Two non-empty multihetsep files are created.

Typical format:

~~~text
Scaffold_34    68306    44     TC
Scaffold_34    87563    259    AG
~~~

The four columns contain:

1. chromosome or scaffold;
2. position of the segregating site;
3. number of callable sites since the previous segregating site;
4. alleles observed on the two haplotypes.

**Check:**

~~~bash
head "$OUTDIR/ABB_${ABB_ID}.${CHROM}.multihetsep.txt"
head "$OUTDIR/SBB_${SBB_ID}.${CHROM}.multihetsep.txt"

wc -l \
  "$OUTDIR/ABB_${ABB_ID}.${CHROM}.multihetsep.txt" \
  "$OUTDIR/SBB_${SBB_ID}.${CHROM}.multihetsep.txt"
~~~

Positions should increase, chromosome labels should equal **Scaffold_34**, and both files should contain records.

If the program reports `invalid literal for int() with base 10: '.'`, a missing genotype remains in the VCF. Return to the input-QC step rather than replacing missing genotypes with **0/0**.

## Step 3 — Run MSMC2 separately for ABB and SBB

**Purpose:** estimate one historical coalescence-rate trajectory from each diploid genome.

**Input:** the two multihetsep files.

~~~bash
msmc2 -t 2 -p '1*2+15*1+1*2' \
  -o "$OUTDIR/ABB_${ABB_ID}" \
  "$OUTDIR/ABB_${ABB_ID}.${CHROM}.multihetsep.txt"

msmc2 -t 2 -p '1*2+15*1+1*2' \
  -o "$OUTDIR/SBB_${SBB_ID}" \
  "$OUTDIR/SBB_${SBB_ID}.${CHROM}.multihetsep.txt"
~~~

| Flag | Meaning |
|---|---|
| **-t 2** | Use two compute threads |
| **-p '1*2+15*1+1*2'** | Group adjacent atomic time intervals that share an estimated rate |
| **-o** | Set the output prefix for the population |

The time-pattern string estimates 17 free rate parameters: the first and last parameters each cover two atomic intervals, while the 15 middle parameters each cover one. With only one scaffold, a simpler pattern can be more stable than a highly parameterized model.

**Expected:** MSMC2 writes several files for each prefix. The principal result tables are:

~~~text
results/msmc2/ABB_4573.final.txt
results/msmc2/SBB_U1916.final.txt
~~~

**Check:**

~~~bash
head "$OUTDIR/ABB_${ABB_ID}.final.txt"
head "$OUTDIR/SBB_${SBB_ID}.final.txt"
~~~

Both tables should contain time boundaries and a coalescence-rate column named **lambda** or **lambda_00**.

## Step 4 — Scale and compare the trajectories

**Purpose:** apply the same mutation rate and generation time to both results so their scales are directly comparable.

Start R:

~~~bash
R
~~~

Then paste:

~~~r
abb <- read.table("results/msmc2/ABB_4573.final.txt", header=TRUE)
sbb <- read.table("results/msmc2/SBB_U1916.final.txt", header=TRUE)

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

abb_scaled <- scale_msmc(abb)
sbb_scaled <- scale_msmc(sbb)

plot(abb_scaled$years, abb_scaled$Ne,
     type="s", log="xy", lwd=2, col="firebrick",
     xlab="Years before present", ylab="Effective population size")
lines(sbb_scaled$years, sbb_scaled$Ne,
      type="s", lwd=2, col="steelblue")
legend("topleft",
       legend=c("ABB: 4573", "SBB: U1916"),
       col=c("firebrick", "steelblue"), lwd=2)
~~~

**Expected:** two stepwise trajectories on identical logarithmic axes.

**Check:** both curves should appear and both axes should contain positive values. If R warns about non-positive values, inspect the corresponding **final.txt** file before interpreting the graph.

The numerical dates and population sizes depend directly on **mu** and **generation_time**. These teaching values must be replaced or justified for a formal analysis.

Exit R without saving the workspace:

~~~r
q(save="no")
~~~

## Interpretation questions

1. Where do the ABB and SBB trajectories begin to differ, and where do they overlap?
2. What aspects of isolation, connectivity, or bottleneck history might explain the contrast?
3. Why should the recent ends of single-genome MSMC2 trajectories be interpreted cautiously?
4. How might using only Scaffold_34 affect the smoothness and uncertainty of the curves?
5. Which assumptions are shared by both curves, and which sources of bias could differ between ABB and SBB?
6. What happens to the time axis if generation time increases? What happens to both axes if the mutation rate changes?

Record the main comparison and its limitations in the [answer sheet](../answers/student_answers.md). Continue to [GONE](02-gone.md).
