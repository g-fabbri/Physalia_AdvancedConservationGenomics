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
  filtered VCF + callable mask
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

mu <- 1.82e-8
generation_time <- 11

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

## Optional extension — block bootstrap

**Purpose:** explore how strongly the inferred trajectories depend on which genomic regions were sampled.

The official [multihetsep bootstrap utility](https://github.com/stschiff/msmc-tools/blob/master/multihetsep_bootstrap.py) resamples genomic blocks with replacement and constructs pseudo-scaffolds. MSMC2 is then run independently on each bootstrap replicate. This extension is computationally longer than the main practical; ten replicates are useful for demonstration, whereas a formal analysis should use substantially more.

### Step A — Choose the block layout

**Input:** the callable mask and the multihetsep files produced above.

Use 5 Mb blocks, the default used by **multihetsep_bootstrap.py**, and calculate how many blocks are needed to approximate the length of Scaffold_34:

~~~bash
CHUNK_SIZE=5000000
N_BOOT=10
BOOTDIR="$OUTDIR/bootstrap"

SCAFFOLD_END=$(zcat "$MASK" |
  awk -v chrom="$CHROM" '$1 == chrom && $3 > end {end=$3} END {print end}')

N_CHUNKS=$(( (SCAFFOLD_END + CHUNK_SIZE - 1) / CHUNK_SIZE ))

mkdir -p "$BOOTDIR"
printf 'Scaffold length: %s bp\nBootstrap blocks per replicate: %s\n' \
  "$SCAFFOLD_END" "$N_CHUNKS"
~~~

**Expected:** a positive scaffold length and a positive number of blocks. The final block may be shorter than 5 Mb, so the reconstructed length is approximate.

### Step B — Generate bootstrap multihetsep files

~~~bash
multihetsep_bootstrap.py \
  -n "$N_BOOT" \
  -s "$CHUNK_SIZE" \
  --chunks_per_chromosome "$N_CHUNKS" \
  --nr_chromosomes 1 \
  --seed 12345 \
  "$BOOTDIR/ABB" \
  "$OUTDIR/ABB_${ABB_ID}.${CHROM}.multihetsep.txt"

multihetsep_bootstrap.py \
  -n "$N_BOOT" \
  -s "$CHUNK_SIZE" \
  --chunks_per_chromosome "$N_CHUNKS" \
  --nr_chromosomes 1 \
  --seed 12345 \
  "$BOOTDIR/SBB" \
  "$OUTDIR/SBB_${SBB_ID}.${CHROM}.multihetsep.txt"
~~~

| Option | Meaning |
|---|---|
| **-n 10** | Generate ten bootstrap replicates |
| **-s 5000000** | Resample blocks of 5 Mb |
| **--chunks_per_chromosome** | Preserve approximately the original scaffold length |
| **--nr_chromosomes 1** | Create one pseudo-scaffold per replicate |
| **--seed 12345** | Make the classroom resampling reproducible |

**Expected:** directories named **ABB_1** to **ABB_10** and **SBB_1** to **SBB_10**, each containing one bootstrap multihetsep file.

**Check:**

~~~bash
find "$BOOTDIR" -name 'bootstrap_multihetsep.chr1.txt' | sort
find "$BOOTDIR" -name 'bootstrap_multihetsep.chr1.txt' | wc -l
~~~

The count should be **20**: ten ABB replicates and ten SBB replicates.

### Step C — Run MSMC2 on every replicate

~~~bash
for REP in $(seq 1 "$N_BOOT"); do
  msmc2 -t 2 -p '1*2+15*1+1*2' \
    -o "$BOOTDIR/ABB_${REP}/ABB_${REP}" \
    "$BOOTDIR/ABB_${REP}/bootstrap_multihetsep.chr1.txt"
done

for REP in $(seq 1 "$N_BOOT"); do
  msmc2 -t 2 -p '1*2+15*1+1*2' \
    -o "$BOOTDIR/SBB_${REP}/SBB_${REP}" \
    "$BOOTDIR/SBB_${REP}/bootstrap_multihetsep.chr1.txt"
done
~~~

**Expected:** each bootstrap directory receives an MSMC2 result, including one **.final.txt** file.

**Check:**

~~~bash
find "$BOOTDIR" -name '*.final.txt' | wc -l
~~~

The count should be **20**. If classroom time is limited, run two or three replicates together and leave the remaining runs as an exercise.

### Step D — Display bootstrap variation

Start R again. The following code writes the figure explicitly to **results/msmc2/MSMC2_ABB_SBB_bootstrap.pdf**:

~~~r
mu <- 4.5e-9
generation_time <- 10

scale_msmc <- function(x) {
  lambda <- if ("lambda_00" %in% names(x)) x$lambda_00 else x$lambda
  midpoint <- sqrt(x$left_time_boundary * x$right_time_boundary)
  y <- data.frame(
    years = midpoint / mu * generation_time,
    Ne = 1 / (2 * mu * lambda)
  )
  y[is.finite(y$years) & is.finite(y$Ne) &
      y$years > 0 & y$Ne > 0, ]
}

read_scaled <- function(filename) {
  scale_msmc(read.table(filename, header=TRUE))
}

abb <- read_scaled("results/msmc2/ABB_4573.final.txt")
sbb <- read_scaled("results/msmc2/SBB_U1916.final.txt")

abb_files <- Sys.glob(
  "results/msmc2/bootstrap/ABB_*/*.final.txt")
sbb_files <- Sys.glob(
  "results/msmc2/bootstrap/SBB_*/*.final.txt")

if (length(abb_files) == 0 || length(sbb_files) == 0) {
  stop("No bootstrap final.txt files found; check BOOTDIR and Step C")
}

abb_boot <- lapply(abb_files, read_scaled)
sbb_boot <- lapply(sbb_files, read_scaled)
all_curves <- c(list(abb, sbb), abb_boot, sbb_boot)

x_limits <- range(unlist(lapply(all_curves, function(x) x$years)))
y_limits <- range(unlist(lapply(all_curves, function(x) x$Ne)))

pdf("results/msmc2/MSMC2_ABB_SBB_bootstrap.pdf",
    width=7, height=5)
plot(abb$years, abb$Ne, type="n", log="xy",
     xlab="Years before present", ylab="Effective population size",
     xlim=x_limits, ylim=y_limits)

for (x in abb_boot) {
  lines(x$years, x$Ne, type="s",
        col=adjustcolor("firebrick", alpha.f=0.20))
}

for (x in sbb_boot) {
  lines(x$years, x$Ne, type="s",
        col=adjustcolor("steelblue", alpha.f=0.20))
}

lines(abb$years, abb$Ne, type="s", lwd=3, col="firebrick")
lines(sbb$years, sbb$Ne, type="s", lwd=3, col="steelblue")
legend("topleft", legend=c("ABB: 4573", "SBB: U1916"),
       col=c("firebrick", "steelblue"), lwd=3)

dev.off()
file.info("results/msmc2/MSMC2_ABB_SBB_bootstrap.pdf")$size
~~~

**Expected:** **dev.off()** prints the name of the closed graphics device, followed by a positive PDF file size. The PDF contains the original ABB and SBB estimates as thick lines surrounded by faint bootstrap trajectories.

Return to the terminal and check the result:

~~~bash
ls -lh results/msmc2/MSMC2_ABB_SBB_bootstrap.pdf
~~~

If R stops with **No bootstrap final.txt files found**, run this in the terminal and compare the paths with Step C:

~~~bash
find results/msmc2/bootstrap -name '*.final.txt'
~~~

Do not call **pdf()** after the plotting commands: that would create a new, empty PDF device.

The spread of ten replicates is a teaching visualization, not a precise confidence interval. Bootstrapping one scaffold measures sensitivity to blocks within that scaffold; it cannot compensate for limited genome coverage, systematic callability bias, or uncertainty in mutation rate and generation time.

## Interpretation questions

1. Where do the ABB and SBB trajectories begin to differ, and where do they overlap?
2. What aspects of isolation, connectivity, or bottleneck history might explain the contrast?
3. Why should the recent ends of single-genome MSMC2 trajectories be interpreted cautiously?
4. How might using only Scaffold_34 affect the smoothness and uncertainty of the curves?
5. Which assumptions are shared by both curves, and which sources of bias could differ between ABB and SBB?
6. What happens to the time axis if generation time increases? What happens to both axes if the mutation rate changes?
7. In which periods are the bootstrap trajectories most variable? What does that imply about confidence in the ABB–SBB contrast?

Record the main comparison and its limitations in the [answer sheet](../answers/student_answers.md). Continue to [GONE](02-gone.md).

