# Part 2 — GONE2: comparing recent history in two bear populations

Estimated practical time: 30 minutes.

GONE2 estimates recent effective population size from linkage disequilibrium (LD) at different recombination distances. We analyse Apennine brown bears (ABB) and Slovak brown bears (SBB) separately and compare their trajectories.

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
    --geno 0.10 \
    --mac 2 \
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
| **--geno 0.10** | Remove sites missing in more than 10% of individuals |
| **--mac 2** | Require at least two copies of the minor allele |

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
  -t 2 \
  -S 1 \
  -E \
  -o "$OUTDIR/ABB" \
  "$COURSE_DIR/$INPUTDIR/ABB_${CHROM}.ped"

"$GONE2_BIN" \
  -g 0 \
  -r 1 \
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
| **-t 2** | Use two threads |
| **-S 1** | Fix the random seed for reproducibility |
| **-E** | Request variation among genetic-algorithm rounds |
| **-o** | Set the output prefix |

The upper recombination fraction is not specified, so GONE2 uses its default **-u 0.05**. Changing `-u` can be explored separately as an instructor sensitivity analysis, not part of the student lesson.

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

Start R from the course root:

~~~bash
R
~~~

~~~r
abb <- read.table("results/gone2/ABB_GONE2_Ne", header=TRUE)
sbb <- read.table("results/gone2/SBB_GONE2_Ne", header=TRUE)

abb_ne <- abb[[grep("^Ne", names(abb), value=TRUE)[1]]]
sbb_ne <- sbb[[grep("^Ne", names(sbb), value=TRUE)[1]]]

pdf("results/gone2/GONE2_ABB_SBB.pdf", width=7, height=5)
plot(abb$Generation, abb_ne,
     type="l", log="y", lwd=2, col="firebrick",
     xlim=rev(range(c(abb$Generation, sbb$Generation))),
     ylim=range(c(abb_ne, sbb_ne), finite=TRUE),
     xlab="Generations before present",
     ylab="Effective population size")
lines(sbb$Generation, sbb_ne,
      lwd=2, col="steelblue")
legend("topright", legend=c("ABB", "SBB"),
       col=c("firebrick", "steelblue"), lwd=2)
dev.off()
~~~

**Expected:** **results/gone2/GONE2_ABB_SBB.pdf**. The present is on the left and older generations are on the right.

## Questions for discussion

1. Which parts of the ABB and SBB trajectories differ?
2. Do the **GONE2_STATS** diagnostics support a panmictic population model?
3. Could relatedness, structure, or SNP ascertainment imitate a recent population change?
4. What information is lost by analysing only one scaffold?
5. Why are 10 ABB and 8 SBB individuals insufficient for strong biological conclusions?
6. Which additional chromosomes, individuals, or maps would most improve the analysis?

Continue to the [synthesis](04-synthesis.md). For further work, try the optional [SMC++](03-smcpp-optional.md) or [currentNe2](03b-currentne2-optional.md) extensions.
