# Part 2 — GONE: comparing recent history in two bear populations

Estimated practical time: 30 minutes.

GONE estimates recent effective population size from linkage disequilibrium (LD) measured at different recombination distances. We will estimate ABB and SBB separately and then compare the trajectories.

In this exercise, we use one chromosome to reduce runtime. The result is a teaching demonstration rather than a complete demographic reconstruction.

The shared QC section has already verified the input VCF, identified 10 Apennine brown bears (ABB) and 8 Slovak brown bears (SBB), and created the two population sample lists. We therefore begin with population subsetting rather than checking the files again.

## Workflow

```text
population VCF
      ↓
separate ABB and SBB samples
      ↓
filter variants and create PED/MAP files
      ↓
adapt the MAP for a one-chromosome run
      ↓
inspect GONE parameters
      ↓
run GONE
      ↓
check and plot the results
```

## Step 1 — Prepare the population VCFs

**Purpose**

Create separate ABB and SBB VCFs from the validated sample lists, then select one population for the first run.

**Input**

- The previously checked 18-individual VCF.
- **data/apennine.samples** and **data/slovak.samples**, created during shared QC.

**Command**

~~~bash
CHROM=Scaffold_34
ALL_VCF=data/UrArMa_18i_s34.vcf.gz
OUTDIR=results/gone
mkdir -p "$OUTDIR"

bcftools view -S data/apennine.samples "$ALL_VCF" \
  -Oz -o data/ABB_s34.vcf.gz
bcftools view -S data/slovak.samples "$ALL_VCF" \
  -Oz -o data/SBB_s34.vcf.gz

bcftools index -t data/ABB_s34.vcf.gz
bcftools index -t data/SBB_s34.vcf.gz
~~~

**Expected output**

- **data/ABB_s34.vcf.gz**, containing the 10 ABB individuals;
- **data/SBB_s34.vcf.gz**, containing the 8 SBB individuals;
- one tabix index for each VCF.

Begin with ABB. After completing the workflow, repeat it with **POPULATION=SBB**:

~~~bash
POPULATION=ABB
VCF=data/${POPULATION}_s34.vcf.gz
PREFIX="$OUTDIR/${POPULATION}_${CHROM}"
~~~

**Check**

~~~bash
printf 'Population: %s\nPopulation VCF: %s\nOutput prefix: %s\n' \
  "$POPULATION" "$VCF" "$PREFIX"
~~~

The population, VCF, and output prefix should all refer to the same assigned population.

## Step 2 — Create the GONE PED/MAP input

**Purpose**

Retain suitable biallelic SNPs and convert the VCF into the PED/MAP format required by this GONE workflow.

**Input**

The population VCF.

**Command**

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

**Flags**

| Flag | Meaning |
|---|---|
| `--vcf` | Read genotypes from the VCF |
| `--double-id` | Use the sample ID as both family and individual ID |
| `--allow-extra-chr` | Accept non-human chromosome labels |
| `--snps-only just-acgt` | Retain canonical A/C/G/T SNPs |
| `--biallelic-only strict` | Retain loci with exactly two alleles |
| `--geno 0.10` | Remove loci missing in more than 10% of individuals |
| `--mac 2` | Require at least two copies of the minor allele |
| `--recode` | Write PED/MAP text files |
| `--out` | Set the output basename |

**Expected output**

```text
results/gone/ABB_Scaffold_34_original_label.ped
results/gone/ABB_Scaffold_34_original_label.map
results/gone/ABB_Scaffold_34_original_label.log
```

PLINK should finish with a message indicating that the PED and MAP files were written.

The example names above are for the ABB run; the SBB run produces the same names beginning with **SBB**.

**Check**

```bash
tail -20 "${PREFIX}_original_label.log"
wc -l "${PREFIX}_original_label.ped"
wc -l "${PREFIX}_original_label.map"
```

- PED rows should equal the number of individuals.
- MAP rows should equal the number of retained SNPs reported by PLINK.
- The log should not contain `Error:`.

### Relabel the selected chromosome

**Purpose**

Make the selected chromosome compatible with the supplied Apennine bear driver.

The driver assumes chromosomes are numbered consecutively beginning with `1`. It reads the last chromosome code as the total number of chromosomes. Therefore, a selected chromosome labelled `12` would incorrectly make the driver expect 12 chromosomes.

Only the chromosome identifier in this derived teaching MAP is changed. Marker IDs and positions remain unchanged.

**Input**

The PED/MAP files produced by PLINK.

**Command**

```bash
cp "${PREFIX}_original_label.ped" "${PREFIX}.ped"

awk 'BEGIN {OFS="\t"} {$1=1; print $1,$2,$3,$4}' \
  "${PREFIX}_original_label.map" > "${PREFIX}.map"
```

**Expected output**

```text
results/gone/ABB_Scaffold_34.ped
results/gone/ABB_Scaffold_34.map
```

Every row in the new MAP should begin with `1`.

**Check**

```bash
cut -f1 "${PREFIX}.map" | sort -u
head "${PREFIX}.map"
wc -l "${PREFIX}.ped" "${PREFIX}.map"
```

- The first command should print only `1`.
- The new and original files should contain the same numbers of individuals and loci.
- MAP column 4 should still contain the original physical positions.

## Step 3 — Set the parameters and run GONE

**Purpose**

Review the main parameters, then calculate LD and infer recent effective population size.

**Input**

The supplied `software/GONE/INPUT_PARAMETERS_FILE`.

The classroom parameter file contains:

| Parameter | Value | Meaning |
|---|---:|---|
| `PHASE` | `2` | Genotype phase is unknown |
| `cMMb` | `1` | Assume 1 cM/Mb when MAP genetic distances are unavailable |
| `DIST` | `1` | Apply the Haldane mapping correction |
| `NGEN` | `2000` | Number of generations represented in LD bins |
| `NBIN` | `400` | Number of LD bins |
| `MAF` | `0.0` | Apply no additional MAF filter inside GONE |
| `ZERO` | `1` | Allow missing genotype codes |
| `maxNCHROM` | `-99` | Analyze all chromosomes detected—one in this exercise |
| `maxNSNP` | `50000` | Approximate maximum SNPs sampled per chromosome |
| `hc` | `0.01` | Maximum recombination fraction analyzed |
| `REPS` | `5` | Number of replicate estimates in the live exercise |
| `threads` | `2` | Number of parallel workers in the live exercise |

The repository already contains these classroom settings. The original Apennine analysis used `REPS=40` and `threads=10`; the smaller classroom values reduce runtime. A precomputed 40-replicate run can be used for the final comparison.

**Command**

Define the course, GONE, and data directories:

```bash
COURSE_DIR=$(pwd)
GONE_DIR="$COURSE_DIR/software/GONE"
DATA_DIR="$COURSE_DIR/results/gone"
FILE="${POPULATION}_${CHROM}"
```

Move into the GONE directory and display the supplied classroom parameter file:

```bash
cd "$GONE_DIR"
cat INPUT_PARAMETERS_FILE
```

**Expected output**

```text
PHASE=2
cMMb=1
DIST=1
NGEN=2000
NBIN=400
MAF=0.0
ZERO=1
maxNCHROM=-99
maxNSNP=50000
hc=0.01
REPS=5
threads=2
```

Students do not need to copy, recreate, or edit this file. The command only makes its settings visible before the run.

### Run GONE

From inside **software/GONE**, calculate LD by recombination-distance bin and infer recent effective population size:

The instructor must have installed the platform-specific programs before class. Confirm that the essential files are executable:

~~~bash
test -x PROGRAMMES/MANAGE_CHROMOSOMES2 &&
test -x PROGRAMMES/LD_SNP_REAL3 &&
test -x PROGRAMMES/SUMM_REP_CHROM3 &&
test -x PROGRAMMES/GONEparallel.sh &&
echo "GONE programs are ready"
~~~

**Expected:** **GONE programs are ready**. If nothing is printed, stop and ask the instructor; the GONE installation is incomplete.

```bash
bash script_GONE.sh "$FILE" "$DATA_DIR"
```

The first argument is the PED/MAP basename without its extension. The second is the directory containing the input and receiving the output.

**Expected output**

The terminal should display stages similar to:

```text
DIVIDE .ped AND .map FILES IN CHROMOSOMES
RUNNING ANALYSIS OF CHROMOSOMES ...
CHROMOSOME ANALYSES took ... seconds
Running GONE
GONE run took ... seconds
END OF ANALYSES
```

The driver internally runs:

| Component | Function |
|---|---|
| `MANAGE_CHROMOSOMES2` | Prepare and optionally subsample chromosome data |
| `LD_SNP_REAL3` | Calculate LD in recombination-distance bins |
| `SUMM_REP_CHROM3` | Standardize and combine LD summaries |
| `GONEparallel.sh` | Perform replicate demographic inference |

**Check**

If the workflow stops, read the final terminal lines. Common causes are missing executables, incompatible compiled binaries, incorrect input basenames, and non-numeric chromosome codes.

## Step 4 — Inspect the outputs

**Purpose**

Confirm successful completion before interpreting the trajectory.

**Input**

Files returned to `results/gone/` by the driver.

**Command**

```bash
cd "$COURSE_DIR"
ls -lh "$DATA_DIR"
```

Principal expected files:

```text
OUTPUT_ABB_Scaffold_34
Output_Ne_ABB_Scaffold_34
Output_d2_ABB_Scaffold_34
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

**Expected output**

- `Output_Ne_${FILE}`: estimates across generations.
- `Output_d2_${FILE}`: observed LD information.
- `outfileHWD`: Hardy–Weinberg deviation diagnostic.
- `timefile`: completed stages and elapsed times.
- `seedfile`: random seed used for SNP sampling.

**Check**

- All principal files should be non-empty.
- The `Ne` output should contain multiple generations.
- `timefile` should finish with `END OF ANALYSES`.

## Step 5 — Plot the GONE trajectory

**Purpose**

Visualize how inferred effective population size changes through recent generations.

**Input**

The ABB or SBB `Output_Ne` file.

**Command**

Inspect the file before deciding whether it contains a header:

```bash
head "$DATA_DIR/Output_Ne_${FILE}"
export GONE_RESULT="$DATA_DIR/Output_Ne_${FILE}"
```

Start R:

```bash
R
```

If the first row contains column names:

```r
gone_file <- Sys.getenv("GONE_RESULT")
gone <- read.table(gone_file, header=TRUE)
names(gone)
head(gone)
```

If the first row is numeric:

```r
gone_file <- Sys.getenv("GONE_RESULT")
gone <- read.table(gone_file, header=FALSE)
names(gone)[1:2] <- c("Generation", "Ne")
```

Plot after replacing the example column names if necessary:

```r
plot(gone$Generation, gone$Ne,
     type="l", log="y", lwd=2,
     xlim=rev(range(gone$Generation)),
     xlab="Generations before present",
     ylab="Effective population size")
```

**Expected output**

A line plot with generations before present on the horizontal axis and effective population size on a logarithmic vertical axis.

**Check**

- Both axes should contain finite positive values.
- The number of plotted rows should match the number of rows in the GONE result.
- Do not interpret abrupt changes until the diagnostics and replicate stability have been examined.

Exit R with `q()`.

After both population runs are available, place or link their final files under `results/gone/comparison/` and overlay them using the same axes. Do not compare curves plotted with different parameter settings.

## Questions for discussion

1. Why does GONE require several individuals while MSMC2 can use one diploid genome?
2. Is `--mac 2` appropriate for this sample size? What changes with `--mac 3`?
3. If MAP column 3 is zero, how strongly is the result dependent on `cMMb=1`?
4. Why might `hc=0.01` have been selected instead of `0.05`?
5. Does requesting 2,000 generations mean that the data resolve the full interval?
6. Could relatives, population structure, or recent migrants imitate a population decline?
7. Does the Hardy–Weinberg diagnostic support the assumption of a panmictic sample?
8. How different are the live five-replicate and precomputed 40-replicate estimates?
9. What information and precision are lost by analyzing only one chromosome?
10. Which portion of the trajectory is sufficiently stable to interpret?

Continue to [NeEstimator](03-neestimator.md).
