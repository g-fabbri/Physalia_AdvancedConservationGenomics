# Part 2 — GONE2: comparing recent history in two bear populations

Estimated practical time: 30 minutes.

## Start your terminal

From the Day 1 course directory, activate the course environment and identify the local GONE2 executable:

~~~bash
conda activate bear-ne-practical
COURSE_DIR=$(pwd)
GONE2_BIN="$COURSE_DIR/software/GONE2/gone2"
~~~

`conda activate` makes BCFtools, PLINK, and the other shared dependencies available. GONE2 is compiled separately inside the course `software/` directory, so `GONE2_BIN` records its complete path. Run this block whenever you open a new terminal.

**Check:**

~~~bash
command -v bcftools
command -v plink
test -x "$GONE2_BIN" && echo "GONE2 ready: $GONE2_BIN"
~~~

The first two commands should print paths inside the `bear-ne-practical` environment. The last command should print the complete path to the executable. If it prints nothing, confirm that you are in the course root and that GONE2 was compiled before class.

> GONE2 estimates recent effective population size from linkage disequilibrium (LD) at different recombination distances. We analyse Apennine brown bears (ABB) and Slovak brown bears (SBB) separately and compare their trajectories.

This is a one-scaffold teaching analysis. The small samples and limited genomic coverage mean that diagnostics are as important as the estimated curves.

~~~text
18-individual VCF
        ↓ separate ABB and SBB
population VCFs
        ↓ filter and convert with PLINK
PED + MAP
        ↓ GONE2
recent Ne trajectories and diagnostics
~~~

Software is prepared before class. Installation instructions are in the [software README](../software/README.md).

## Step 1 — Prepare ABB and SBB VCFs

**Purpose:** separate the two populations using the sample lists created during shared QC.

**Input:** the checked 18-individual VCF and the ABB/SBB sample lists.

~~~bash
CHROM=Scaffold_25
ALL_VCF=data/UrArMa_18i_s25.vcf.gz
ABB_VCF=data/ABB_s25.vcf.gz
SBB_VCF=data/SBB_s25.vcf.gz

bcftools view -S data/apennine.samples "$ALL_VCF" \
  -Oz -o "$ABB_VCF"
bcftools view -S data/slovak.samples "$ALL_VCF" \
  -Oz -o "$SBB_VCF"

bcftools index -t "$ABB_VCF"
bcftools index -t "$SBB_VCF"
~~~

**Expected:** **ABB_s25.vcf.gz** contains the 10 ABB individuals and **SBB_s25.vcf.gz** contains the 8 SBB individuals.

## Step 2 — Create the PED/MAP files

**Purpose:** apply the same variant filters to each population and create the PLINK files read by GONE2.

PLINK's text format uses two matching files with the same prefix:

- the **PED file** contains one row per individual. Its first six columns describe the sample and family, followed by two allele entries for every marker;
- the **MAP file** contains one row per marker, giving its chromosome, marker ID, genetic-map position, and physical base-pair position.

The marker order in the MAP file must exactly match the genotype order in every PED row. Together, the files tell GONE2 which alleles each bear carries and where those markers occur along the scaffold.

~~~bash
INPUTDIR=results/gone2/input
mkdir -p "$INPUTDIR"

for POPULATION in ABB SBB; do
  VCF="data/${POPULATION}_s25.vcf.gz"
  PREFIX="$INPUTDIR/${POPULATION}_${CHROM}"

  plink \
    --vcf "$VCF" \
    --double-id \
    --allow-extra-chr \
    --snps-only just-acgt \
    --biallelic-only strict \
    --recode \
    --out "${PREFIX}_original_label"

  cp "${PREFIX}_original_label.ped" "${PREFIX}.ped"
  awk 'BEGIN {OFS="\t"} {$1=1; print $1,$2,$3,$4}' \
    "${PREFIX}_original_label.map" > "${PREFIX}.map"
done
~~~

The loop performs exactly the same commands for ABB and SBB. The derived MAP uses chromosome code **1**, which GONE2 accepts. Marker IDs and physical positions are unchanged.

| PLINK option | Meaning |
|---|---|
| **--snps-only just-acgt** | Retain canonical A/C/G/T SNPs |
| **--biallelic-only strict** | Retain strictly biallelic sites |

**Expected:**

~~~text
results/gone2/input/ABB_Scaffold_25.ped
results/gone2/input/ABB_Scaffold_25.map
results/gone2/input/SBB_Scaffold_25.ped
results/gone2/input/SBB_Scaffold_25.map
~~~

**Check:**

~~~bash
wc -l "$INPUTDIR/ABB_${CHROM}.ped" "$INPUTDIR/SBB_${CHROM}.ped"
cut -f1 "$INPUTDIR/ABB_${CHROM}.map" | sort -u
cut -f1 "$INPUTDIR/SBB_${CHROM}.map" | sort -u
~~~

The PED files should contain 10 and 8 rows. Each MAP check should print only **1**.

## Step 3 — Run GONE2

**Purpose:** estimate a recent Ne trajectory for each population under identical settings.

~~~bash
COURSE_DIR=$(pwd)
GONE2_BIN="$COURSE_DIR/software/GONE2/gone2"
OUTDIR="$COURSE_DIR/results/gone2"

"$GONE2_BIN" \
  -g 0 \
  -r 1 \
  -u 0.02 \
  -t 2 \
  -S 1 \
  -E \
  -o "$OUTDIR/ABB" \
  "$COURSE_DIR/$INPUTDIR/ABB_${CHROM}.ped"

"$GONE2_BIN" \
  -g 0 \
  -r 1 \
  -u 0.02 \
  -t 2 \
  -S 1 \
  -E \
  -o "$OUTDIR/SBB" \
  "$COURSE_DIR/$INPUTDIR/SBB_${CHROM}.ped"
~~~

| Option | Meaning |
|---|---|
| **-g 0** | Treat genotypes as unphased diploids |
| **-r 1** | Assume a constant recombination rate of 1 cM/Mb |
| **-u 0.02** | Use LD bins only up to a recombination fraction of 0.02 |
| **-t 2** | Use two threads |
| **-S 1** | Fix the random seed for reproducibility |
| **-E** | Request variation among genetic-algorithm rounds |
| **-o** | Set the output prefix |

### Why use `-u 0.02`?

GONE2 relates LD between marker pairs to the time in the past that generated that LD. As a useful approximation, LD at recombination fraction `c` is most informative about approximately `1/(2c)` generations ago. The default upper bound is `0.05`, corresponding roughly to information from about 10 generations ago. Here we set the upper bound to `0.02`, corresponding roughly to **25 generations ago**.

This choice excludes the most weakly linked marker pairs and focuses the fit on more tightly linked pairs. For this one-scaffold exercise, it reduces sensitivity to noisy LD at larger distances and to error from assuming a uniform recombination rate of 1 cM/Mb instead of using a detailed genetic map. The trade-off is important: `-u 0.02` provides **less information about the most recent generations** than the default `0.05`; it is not universally a better value. It should be reported explicitly, and a formal analysis should compare plausible `-u` settings and, where possible, use an empirical recombination map.

Do not use **-x** in the main run. It fits a structured metapopulation model, while this exercise initially treats ABB and SBB as separate populations.

**Expected:** progress messages for reading the data, measuring d², estimating Ne, and writing output files.

## Step 4 — Inspect the results

~~~bash
head "$OUTDIR/ABB_GONE2_Ne"
head "$OUTDIR/SBB_GONE2_Ne"
cat "$OUTDIR/ABB_GONE2_STATS"
cat "$OUTDIR/SBB_GONE2_STATS"
~~~

Principal outputs:

| File ending | Content |
|---|---|
| **GONE2_Ne** | Generation and estimated Ne |
| **GONE2_d2** | Observed and predicted LD by recombination bin |
| **GONE2_STATS** | Inputs, parameters, Hardy–Weinberg diagnostics, runtime, and warnings |

If only a **GONE2_STATS** file is produced, read its failure explanation. Before interpreting a curve, examine the number of individuals and SNPs, estimated Fis, inferred genome length, and any population-structure warning.

## Step 5 — Plot ABB and SBB

**Purpose:** compare ABB and SBB over a chosen number of generations and save the figure directly as a PDF.

The provided [R script](../scripts/plot_gone2.R) reads both `GONE2_Ne` tables, finds the Ne column, removes non-positive or non-finite values, restricts the trajectories to the requested generations, and draws them with a logarithmic Ne axis. The first command-line number is the youngest generation to display and the second is the oldest.

For example, plot generations 1–100 from the course root:

~~~bash
Rscript scripts/plot_gone2.R 1 100
~~~

To choose another interval, replace `1 100`; for example, `10 75` plots generations 10–75. With no numbers, the script defaults to generations 1–100.

**Expected:** the script prints the selected interval and creates **results/gone2/GONE2_ABB_SBB_generations_1_100.pdf**. The present is toward the left and older generations are toward the right.

**Check:**

~~~bash
ls -lh results/gone2/GONE2_ABB_SBB_generations_1_100.pdf
~~~

The PDF should have a nonzero size and contain both coloured trajectories. Do not request generations outside those available in the two result tables.

## Questions for discussion

1. Which parts of the ABB and SBB trajectories differ?
2. Do the **GONE2_STATS** diagnostics support a panmictic population model?
3. Could relatedness, structure, or SNP ascertainment imitate a recent population change?
4. What information is lost by analysing only one scaffold?
5. Why are 10 ABB and 8 SBB individuals insufficient for strong biological conclusions?
6. Which additional chromosomes, individuals, or maps would most improve the analysis?

Continue to the [synthesis](04-synthesis.md). For further work, try the optional [SMC++](03-smcpp-optional.md) or [currentNe2](03b-currentne2-optional.md) extensions.
