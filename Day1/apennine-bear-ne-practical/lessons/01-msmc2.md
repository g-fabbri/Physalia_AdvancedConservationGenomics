# Part 1 — MSMC2: comparing historical population-size trajectories

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


## Step 1 — Define the analysis variables

**Purpose:** assign short, meaningful names to the prepared inputs and create an output directory.

**Input:** the two already-filtered single-sample VCFs, the common callable mask, and Scaffold_25.

~~~bash
CHROM=Scaffold_25
ABB_VCF=data/UrArMa_4573_s25.vcf.gz
SBB_VCF=data/UrArMa_U1916_s25.vcf.gz
ABB_ID=4573
SBB_ID=U1916
MASK=data/UrArMa_callable_s25.bed.gz
OUTDIR=results/msmc2

mkdir -p "$OUTDIR"
~~~

**Expected:** nothing is printed. The commands assign shell variables and create **results/msmc2** if necessary.

Run these commands from the course root and keep using the same terminal: the variables are needed in Steps 2 and 3. The VCFs already contain one selected bear each; **OUTDIR** keeps intermediate files and final estimates together without changing the input VCFs.

**Check:**

~~~bash
printf 'ABB: %s (%s)\nSBB: %s (%s)\nMask: %s\n' \
  "$ABB_ID" "$ABB_VCF" "$SBB_ID" "$SBB_VCF" "$MASK"
~~~

Expected shape:

~~~text
ABB: 4573 (data/UrArMa_4573_s25.vcf.gz)
SBB: U1916 (data/UrArMa_U1916_s25.vcf.gz)
Mask: data/UrArMa_callable_s25.bed.gz
~~~

## Step 2 — Create the multihetsep files

**Purpose:** combine genotype information and callable sequence in the format required by MSMC2.

**Input:** one already-filtered single-sample VCF and the common positive mask for each run.

MSMC2 needs to know both where heterozygous sites occur and how much sequence could have been observed between them. The VCF supplies genotypes at variant positions; the BED mask supplies the callable intervals. We run the conversion separately because ABB and SBB have different genotypes, even though they use the same teaching mask.

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
| **--chr "$CHROM"** | Process only Scaffold_25 |
| **--mask "$MASK"** | Count only positions included in the callable mask |
| **"$ABB_VCF" / "$SBB_VCF"** | Read the prepared diploid genome for that population |
| **> output file** | Save the generated multihetsep text |

**Expected:** each run reports that it is generating input for **2 haplotypes**, because each VCF contains one diploid individual. Two non-empty multihetsep files are created.

Typical format:

~~~text
Scaffold_25    68306    44     TC
Scaffold_25    87563    259    AG
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

Positions should increase, chromosome labels should equal **Scaffold_25**, and both files should contain records.

The number of lines is the number of multihetsep records, not the number of callable bases. The third column carries information about callable sequence between records, which would be lost if we simply handed MSMC2 a list of SNP positions.

If the program reports `invalid literal for int() with base 10: '.'`, a missing genotype remains in the VCF. Return to the input-QC step rather than replacing missing genotypes with **0/0**.

## Step 3 — Run MSMC2 separately for ABB and SBB

**Purpose:** estimate one historical coalescence-rate trajectory from each diploid genome.

**Input:** the two multihetsep files.

Each run treats the two haplotypes of one bear as the genetic sample. The `-p` pattern constrains adjacent time intervals to share rates, reducing the number of independently fitted values for this one-scaffold demonstration. Use the **same** pattern and thread count for ABB and SBB so the comparison does not also change the model settings.

~~~bash
msmc2_Linux -t 2 -p '1*2+15*1+1*2' \
  -o "$OUTDIR/ABB_${ABB_ID}" \
  "$OUTDIR/ABB_${ABB_ID}.${CHROM}.multihetsep.txt"

msmc2_Linux -t 2 -p '1*2+15*1+1*2' \
  -o "$OUTDIR/SBB_${SBB_ID}" \
  "$OUTDIR/SBB_${SBB_ID}.${CHROM}.multihetsep.txt"
~~~

| Flag | Meaning |
|---|---|
| **-t 2** | Use two compute threads |
| **-p '1x2+15x1+1x2'** | Group adjacent atomic time intervals that share an estimated rate |
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

Both tables should contain time boundaries and a coalescence-rate column named **lambda**.

The time boundaries and rates are in MSMC2's scaled units. They are not yet calendar years or directly readable values of effective population size.

## Step 4 — Scale and compare the trajectories

**Purpose:** convert the fitted time boundaries to years and the coalescence rate to effective population size, then save an ABB–SBB comparison figure.

**Input:** the two **.final.txt** files from Step 3 and the provided [plotting script](../scripts/plot_msmc2.R).

The script uses a mutation rate of **1.82 × 10⁻⁸ per site per generation** and a generation time of **11 years** for both bears. It uses the geometric midpoint of each fitted time interval, calculates `Ne = 1 / (2 × mu × lambda)`, and plots only positive, finite values on logarithmic axes. These are explicit scaling assumptions, not values estimated by MSMC2; justify or revise them for a formal analysis.

From the course root, run:

~~~bash
Rscript scripts/plot_msmc2.R
~~~

**Expected:** the script prints the path **results/msmc2/MSMC2_ABB_4573_SBB_U1916.pdf**. The PDF contains two stepwise trajectories on the same axes. No interactive R session or `q(save="no")` command is needed.

**Check:**

~~~bash
ls -lh results/msmc2/MSMC2_ABB_4573_SBB_U1916.pdf
~~~

The file should have a nonzero size. Inspect the figure: both curves should appear, and the time and Ne axes should contain positive values. The youngest and oldest intervals are often less reliable than the central part of a one-scaffold trajectory.

## Optional extension — block bootstrap

**Purpose:** explore how strongly the inferred trajectories depend on which genomic regions were sampled.

The official [multihetsep bootstrap utility](https://github.com/stschiff/msmc-tools/blob/master/multihetsep_bootstrap.py) resamples genomic blocks with replacement and constructs pseudo-scaffolds. MSMC2 is then run independently on each bootstrap replicate. This extension is computationally longer than the main practical; ten replicates are useful for demonstration, whereas a formal analysis should use substantially more.

The bootstrap begins from the **multihetsep files**, not from the original VCFs. Each resampled file represents a different selection of blocks from the same scaffold; comparing the fitted curves reveals sensitivity to that selection.

### Step A — Choose the block layout

**Input:** the callable mask and the multihetsep files produced above.

Use 5 Mb blocks, the default used by **multihetsep_bootstrap.py**, and calculate how many blocks are needed to approximate the length of Scaffold_25:

The BED intervals can have gaps, so summing callable bases would not give the scaffold coordinate span. Here we take the largest BED end coordinate on **Scaffold_25**, then round up to the number of 5 Mb blocks needed to cover that span.

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

**Check:** if either value is zero, confirm that the mask contains **Scaffold_25** and that **CHROM** and **MASK** still have the values from Step 1.

### Step B — Generate bootstrap multihetsep files

**Purpose:** make ten pseudo-scaffolds for each bear by drawing 5 Mb blocks with replacement. The fixed seed makes the block selections reproducible. ABB and SBB are resampled separately because their heterozygous sites differ.

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

**Purpose:** fit the same MSMC2 time model to each resampled input. A bootstrap replicate is not a new animal; it is another block sample from the original scaffold. The loops may take substantially longer than the two main runs.

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

**Purpose:** draw the original ABB and SBB estimates as thick lines and the bootstrap estimates as faint lines. The [bootstrap plotting script](../scripts/plot_msmc2_bootstrap.R) uses the same mutation rate and generation time as the main plot. It checks that bootstrap results exist **before** opening a PDF device, avoiding an empty PDF if Step C was skipped.

From the course root, run:

~~~bash
Rscript scripts/plot_msmc2_bootstrap.R
~~~

**Expected:** the script reports how many ABB and SBB replicates it read and prints the path **results/msmc2/MSMC2_ABB_SBB_bootstrap.pdf**.

**Check:**

~~~bash
ls -lh results/msmc2/MSMC2_ABB_SBB_bootstrap.pdf
~~~

If R stops with **Bootstrap results not found**, compare the actual paths with Step C:

~~~bash
find results/msmc2/bootstrap -name '*.final.txt'
~~~

The spread of ten replicates is a teaching visualization, not a precise confidence interval. Bootstrapping one scaffold measures sensitivity to blocks within that scaffold; it cannot compensate for limited genome coverage, systematic callability bias, or uncertainty in mutation rate and generation time.

## Interpretation questions

1. Where do the ABB and SBB trajectories begin to differ, and where do they overlap?
2. What aspects of isolation, connectivity, or bottleneck history might explain the contrast?
3. Why should the recent ends of single-genome MSMC2 trajectories be interpreted cautiously?
4. How might using only Scaffold_25 affect the smoothness and uncertainty of the curves?
5. Which assumptions are shared by both curves, and which sources of bias could differ between ABB and SBB?
6. What happens to the time axis if generation time increases? What happens to both axes if the mutation rate changes?
7. In which periods are the bootstrap trajectories most variable? What does that imply about confidence in the ABB–SBB contrast?

Record the main comparison and its limitations in the [answer sheet](../answers/student_answers.md). Continue to [GONE2](02-gone2.md).
