# Optional Part 2B — original GONE comparison

Estimated practical time: 20–30 minutes.

This optional lesson runs the original GONE workflow on the same ABB and SBB PED/MAP files prepared in [Part 2 — GONE2](02-gone2.md). It is retained while the two implementations are being tested.

GONE estimates recent effective population size from LD bins but uses a shell driver, several compiled programs, and an external parameter file. Software installation is described only in the [software README](../software/README.md).

## Step 1 — Reuse the GONE2 inputs

**Purpose:** copy the prepared files into a separate output directory so original GONE does not mix its outputs with GONE2.

~~~bash
COURSE_DIR=$(pwd)
INPUTDIR="$COURSE_DIR/results/gone2/input"
ABB_DIR="$COURSE_DIR/results/gone/ABB"
SBB_DIR="$COURSE_DIR/results/gone/SBB"
GONE_DIR="$COURSE_DIR/software/GONE"

mkdir -p "$ABB_DIR" "$SBB_DIR"
cp "$INPUTDIR"/ABB_Scaffold_34.ped "$ABB_DIR"/
cp "$INPUTDIR"/ABB_Scaffold_34.map "$ABB_DIR"/
cp "$INPUTDIR"/SBB_Scaffold_34.ped "$SBB_DIR"/
cp "$INPUTDIR"/SBB_Scaffold_34.map "$SBB_DIR"/
~~~

**Expected:** separate ABB and SBB directories under **results/gone**. This prevents the driver's generic diagnostic files from overwriting one another.

## Step 2 — Inspect the parameters and run GONE

Move into the GONE directory:

~~~bash
cd "$GONE_DIR"
cat INPUT_PARAMETERS_FILE
~~~

The classroom file contains:

| Parameter | Value | Meaning |
|---|---:|---|
| **PHASE** | **2** | Phase is unknown |
| **cMMb** | **1** | Assume 1 cM/Mb |
| **DIST** | **1** | Apply the Haldane mapping correction |
| **NGEN** | **2000** | Generations represented in the LD bins |
| **NBIN** | **400** | Number of LD bins |
| **MAF** | **0.0** | No additional MAF filter inside GONE |
| **maxNSNP** | **50000** | Approximate maximum SNPs sampled per chromosome |
| **hc** | **0.01** | Maximum recombination fraction |
| **REPS** | **5** | Replicate estimates for the live exercise |
| **threads** | **2** | Parallel workers |

Run ABB:

~~~bash
FILE=ABB_Scaffold_34
bash script_GONE.sh "$FILE" "$ABB_DIR"
~~~

Then run SBB:

~~~bash
FILE=SBB_Scaffold_34
bash script_GONE.sh "$FILE" "$SBB_DIR"
~~~

**Expected:** each run reports chromosome preparation, LD calculation, GONE fitting, and **END OF ANALYSES**.

## Step 3 — Inspect the outputs

~~~bash
cd "$COURSE_DIR"
head results/gone/ABB/Output_Ne_ABB_Scaffold_34
head results/gone/SBB/Output_Ne_SBB_Scaffold_34
cat results/gone/ABB/OUTPUT_ABB_Scaffold_34
cat results/gone/SBB/OUTPUT_SBB_Scaffold_34
~~~

Principal outputs:

| File | Meaning |
|---|---|
| **Output_Ne_...** | Ne estimates across generations |
| **Output_d2_...** | Observed LD information |
| **OUTPUT_...** | Run summary and diagnostic information |

Do not compare GONE and GONE2 as if every setting were identical. In particular, GONE uses **hc=0.01**, whereas GONE2 defaults to an upper recombination fraction of 0.05.

## Step 4 — Plot the original GONE trajectories

Inspect the first rows to determine whether the files contain headers:

~~~bash
head results/gone/ABB/Output_Ne_ABB_Scaffold_34
head results/gone/SBB/Output_Ne_SBB_Scaffold_34
~~~

Start R and adjust **header** if required:

~~~r
abb <- read.table("results/gone/ABB/Output_Ne_ABB_Scaffold_34",
                  header=TRUE)
sbb <- read.table("results/gone/SBB/Output_Ne_SBB_Scaffold_34",
                  header=TRUE)

names(abb)[1:2] <- c("Generation", "Ne")
names(sbb)[1:2] <- c("Generation", "Ne")

pdf("results/gone/GONE_ABB_SBB.pdf", width=7, height=5)
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

## Questions for comparison

1. Do GONE and GONE2 infer the same direction of population change?
2. Which time periods are stable across the two implementations?
3. How might different recombination limits and fitting procedures affect the curves?
4. Do the diagnostics support treating ABB and SBB as panmictic populations?
5. Which conclusions disappear when only one scaffold and these small samples are used?

Return to [Part 2 — GONE2](02-gone2.md) or continue to [Part 3 — NeEstimator](03-neestimator.md).
