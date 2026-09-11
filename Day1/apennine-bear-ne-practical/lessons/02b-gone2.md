# Optional Part 2B — GONE2: testing the updated LD method

Estimated practical time: 20–30 minutes, excluding installation.

GONE2 is the newer implementation of LD-based demographic inference. We keep the original GONE lesson and run GONE2 separately so that the two programs can be tested on the same ABB and SBB data.

This lesson uses the PED/MAP files already created in [Part 2 — GONE](02-gone.md):

~~~text
results/gone/ABB_Scaffold_34.ped
results/gone/ABB_Scaffold_34.map
results/gone/SBB_Scaffold_34.ped
results/gone/SBB_Scaffold_34.map
~~~

GONE2 can read PED/MAP directly; it does not use the original GONE driver or **INPUT_PARAMETERS_FILE**.

## Instructor setup — install GONE2

Installation is completed once before class, not by every student:

~~~bash
git clone --depth 1 https://github.com/esrud/GONE2.git software/GONE2
cd software/GONE2
make gone
cd ../..
~~~

**Check:**

~~~bash
test -x software/GONE2/gone2 && echo "GONE2 is ready"
~~~

**Expected:** **GONE2 is ready**. On a cluster, the instructor may need to load a suitable C++ compiler before running **make gone**.

## Step 1 — Define the inputs

**Purpose:** use the same population files and mapping assumptions used for the original GONE test.

~~~bash
COURSE_DIR=$(pwd)
GONE2_BIN="$COURSE_DIR/software/GONE2/gone2"
INPUTDIR="$COURSE_DIR/results/gone"
OUTDIR="$COURSE_DIR/results/gone2"

ABB_PED="$INPUTDIR/ABB_Scaffold_34.ped"
SBB_PED="$INPUTDIR/SBB_Scaffold_34.ped"

mkdir -p "$OUTDIR"
~~~

**Expected:** no terminal output.

## Step 2 — Run ABB and SBB

**Purpose:** infer the recent Ne trajectory independently for each population.

~~~bash
"$GONE2_BIN" \
  -g 0 \
  -r 1 \
  -u 0.01 \
  -t 2 \
  -S 1 \
  -E \
  -o "$OUTDIR/ABB" \
  "$ABB_PED"

"$GONE2_BIN" \
  -g 0 \
  -r 1 \
  -u 0.01 \
  -t 2 \
  -S 1 \
  -E \
  -o "$OUTDIR/SBB" \
  "$SBB_PED"
~~~

| Option | Meaning |
|---|---|
| **-g 0** | Treat the genotypes as unphased diploids |
| **-r 1** | Assume a constant recombination rate of 1 cM/Mb |
| **-u 0.01** | Use 0.01 as the upper recombination fraction, matching **hc=0.01** in the GONE test |
| **-t 2** | Use two threads |
| **-S 1** | Fix the random seed for a reproducible classroom run |
| **-E** | Include the range among genetic-algorithm rounds in the Ne output |
| **-o** | Set the output prefix |

We do not add a second MAF filter because the PED/MAP inputs were already filtered with PLINK. Do not use **-x** here: that option fits a structured metapopulation model, whereas ABB and SBB are being analysed separately as populations.

**Expected:** progress messages for input reading, d² calculation, Ne estimation, and output writing.

## Step 3 — Inspect the results

~~~bash
ls -lh "$OUTDIR"
head "$OUTDIR/ABB_GONE2_Ne"
head "$OUTDIR/SBB_GONE2_Ne"
cat "$OUTDIR/ABB_GONE2_STATS"
cat "$OUTDIR/SBB_GONE2_STATS"
~~~

The current GONE2 release normally creates:

~~~text
ABB_GONE2_Ne
ABB_GONE2_d2
ABB_GONE2_STATS
SBB_GONE2_Ne
SBB_GONE2_d2
SBB_GONE2_STATS
~~~

- **GONE2_Ne** contains Generation and Ne; with **-E**, it also contains the minimum and maximum across genetic-algorithm rounds.
- **GONE2_d2** contains observed and predicted LD by recombination bin.
- **GONE2_STATS** records the command, runtime, input summary, Hardy–Weinberg information, and warnings.

If only a **GONE2_STATS** file is produced, read it carefully: GONE2 records why the analysis could not estimate a trajectory.

## Step 4 — Plot the GONE2 trajectories

Start R from the course root:

~~~bash
R
~~~

~~~r
abb <- read.table("results/gone2/ABB_GONE2_Ne", header=TRUE)
sbb <- read.table("results/gone2/SBB_GONE2_Ne", header=TRUE)

names(abb)
head(abb)

pdf("results/gone2/GONE2_ABB_SBB.pdf", width=7, height=5)
plot(abb$Generation, abb$Ne,
     type="l", log="y", lwd=2, col="firebrick",
     xlim=rev(range(c(abb$Generation, sbb$Generation))),
     ylim=range(c(abb$Ne, sbb$Ne), finite=TRUE),
     xlab="Generations before present",
     ylab="Effective population size")
lines(sbb$Generation, sbb$Ne,
      lwd=2, col="steelblue")
legend("topright", legend=c("ABB", "SBB"),
       col=c("firebrick", "steelblue"), lwd=2)
dev.off()
~~~

**Expected:** **results/gone2/GONE2_ABB_SBB.pdf**, with the present on the left and older generations on the right.

If the Ne column has a different capitalization in the tested release, use the names printed by **names(abb)** in the plotting commands.

## Questions for comparison

1. Do GONE and GONE2 infer the same direction and approximate timing of recent changes?
2. Are differences caused by the data, or by differences in model fitting and regularisation?
3. Are the GONE2 minimum and maximum estimates narrow enough to support the main pattern?
4. What warnings appear in **GONE2_STATS** for ABB and SBB?
5. How do the small samples—10 ABB and 8 SBB individuals—and the use of one scaffold limit both methods?
6. Why is matching the recombination assumptions important when comparing GONE with GONE2?

Return to [Part 2 — GONE](02-gone.md) or continue to [Part 3 — NeEstimator](03-neestimator.md).
