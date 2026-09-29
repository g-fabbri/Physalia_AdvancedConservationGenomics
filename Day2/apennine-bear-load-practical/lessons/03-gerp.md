# Optional Part 3 — Examine GERP conservation scores

## What are we doing in this lesson?

In this lesson, we will use **GERP scores** to identify derived variants located at evolutionarily conserved genomic positions.

A positive GERP score indicates that a position has changed less often across species than expected under neutrality. Higher scores therefore suggest stronger evolutionary constraint and potentially greater functional importance.

GERP and SnpEff provide complementary information:

- **SnpEff** predicts how a variant may affect an annotated gene or transcript;
- **GERP** measures how conserved the genomic position is across species;

By combining these results, we can ask whether ABB and SBB differ in the number of derived variants found at strongly conserved positions.

### Choosing an available GERP track

GERP scores are calculated from multi-species genome alignments and are normally distributed as precomputed tracks for particular reference species and genome assemblies.

Before starting an analysis, we should check:

1. whether a GERP track already exists for our study species;
2. which genome assembly and version it uses;
3. whether its chromosome or scaffold names match our data.

If a score track exists for the same species and assembly used by the VCF, the positions can be compared directly.

In our case, a precomputed GERP track is not available for the Apennine-bear `mUrsArc1.1` assembly. We therefore use a track from a closely related species: the polar bear. Ensembl provides a 91-mammal GERP track using the polar-bear `UrsMar_1.0` assembly.

Because our SNPs and the GERP scores use different assemblies, we must align Apennine Scaffold 25 to the polar-bear genome and translate the SNP coordinates before extracting the scores.

Using the polar-bear track does not mean that the polar-bear allele is assumed to be ancestral. The polar-bear assembly is used only as the coordinate system for the available precomputed GERP scores.

A high GERP score does not prove that a variant is deleterious. It provides evidence that the position may be functionally important, but it does not measure selection coefficients, dominance or actual fitness effects.

## Overview of the analysis

The analysis has seven main steps:

1. **Align the assemblies:** align Apennine Scaffold 25 to the polar-bear genome.
2. **Create a chain file:** convert the alignment into a coordinate map that `liftOver` can use.
3. **Lift the SNP positions:** translate the Apennine SNP coordinates into polar-bear coordinates.
4. **Extract GERP scores:** retrieve conservation scores from the polar-bear bigWig.
5. **Return to Apennine coordinates:** attach each score to its original Scaffold 25 position.
6. **Join with GenoLoader:** combine the GERP scores with polarized genotypes.
7. **Summarize and plot:** compare constrained derived variants between ABB and SBB.

The data move through the workflow as follows:

```text
Apennine SNP position
        ↓ liftOver
Polar-bear position
        ↓ bigWig lookup
Polar-bear GERP score
        ↓ original site ID
Apennine SNP position + GERP score
        ↓ join with GenoLoader
GERP score + derived genotypes
```


## What is a bigWig?

A **bigWig** (`.bw`) is a compressed and indexed binary file containing numerical values along a genome. In this exercise, those values are GERP conservation scores.

Unlike a text BED or bedGraph file, a bigWig is not intended to be read with commands such as `head`. Its internal index allows `bigWigToBedGraph` to extract a selected genomic interval without converting the complete file.

The resulting bedGraph is a readable text table with four columns:

```text
chromosome    start    end    score
```

More information is available in the [UCSC bigWig guide](https://genome.ucsc.edu/goldenPath/help/bigWig).

## Why do we align to the polar-bear genome?

Our SNPs use coordinates from the **Apennine-bear mUrsArc1.1 assembly**, whereas the available GERP scores use coordinates from the **polar-bear UrsMar_1.0 assembly**.

For example, the Apennine position:

```text
Scaffold_25:100000
```

cannot be found directly in the polar-bear GERP file because the two assemblies use different sequence names and coordinates.

We use the Ensembl release 114 [91-mammal GERP bigWig for polar bear](https://ftp.ensembl.org/pub/release-114/compara/conservation_scores/91_mammals.gerp_conservation_score/):

```text
gerp_conservation_scores.ursus_maritimus.UrsMar_1.0.bw
```

We also use the polar-bear FASTA from the [same Ensembl release](https://ftp.ensembl.org/pub/release-114/fasta/ursus_maritimus/dna/):

```text
Ursus_maritimus.UrsMar_1.0.dna.toplevel.fa
```

The FASTA and bigWig both use GenBank-style sequence identifiers such as:

```text
AVOR01003489.1
```

This agreement is essential because `bigWigToBedGraph` matches sequence names exactly. The RefSeq version of the same assembly uses identifiers beginning with `NW_` and would return empty results unless those names were translated.

We therefore align Apennine Scaffold 25 to the polar-bear assembly, use `liftOver` to translate each SNP position, and retrieve the GERP score from the corresponding polar-bear coordinate. A stable site identifier then allows us to return the score to the original Apennine position.

The alignment is used only to translate coordinates. It does not calculate new GERP scores, and it does not assume that the polar-bear allele is ancestral.

## Start your terminal

From your Day 2 directory, run:

```bash
conda activate bear-load-practical

COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$PATH"
```

Define the input files:

```bash
CHROM=Scaffold_25

APP_FA="$COURSE_DIR/data/mUrsArc1.1.genome.s25.fasta"
POLAR_FA="$COURSE_DIR/gerp_input/Ursus_maritimus.UrsMar_1.0.dna.toplevel.fa"

VCF="$COURSE_DIR/data/Bears_4pops_s25.vcf.gz"
GT="$COURSE_DIR/results/snpeff/Bears_4pops_s25.ann.POP_OUT.gt"

GERP_BW="$COURSE_DIR/gerp_input/gerp_conservation_scores.ursus_maritimus.UrsMar_1.0.bw"
```

The variables identify:

- `CHROM`: the Apennine scaffold analysed in the lesson;
- `APP_FA`: the Apennine reference sequence for Scaffold 25;
- `POLAR_FA`: the polar-bear genome matching the GERP track;
- `VCF`: the Scaffold 25 variants from the four bear populations;
- `GT`: the polarized genotype table produced by GenoLoader;
- `GERP_BW`: the polar-bear bigWig containing the GERP scores.

Define the output directories:

```bash
OUTDIR="$COURSE_DIR/results/gerp"
PREP="$OUTDIR/preparation"

mkdir -p "$PREP"
```

- `OUTDIR` stores the final tables and figures;
- `PREP` stores intermediate alignment, liftOver and score-extraction files.

Check that the required programs are available:

```bash
command -v \
  minimap2 \
  transanno \
  liftOver \
  bigWigInfo \
  bigWigToBedGraph \
  bedtools
```

This should print the path of every program. If one is missing, check that the correct Conda environment is active.

## Step 1 — Map Apennine Scaffold 25 to the polar-bear genome

**Purpose:** create a coordinate map between Apennine Scaffold 25 and the polar-bear assembly used by the GERP track.

### 1.1 — Check the bigWig sequence names

```bash
bigWigInfo -chroms "$GERP_BW" | head -n 30
```

The bigWig should contain polar-bear sequence names beginning with identifiers such as `AVOR` and 'KK'.

Check the polar-bear FASTA:

```bash
grep '^>' "$POLAR_FA" | head
```

The FASTA and bigWig must use matching sequence names.

### 1.2 — Align the two assemblies

```bash
minimap2 -cx asm20 --cs -t 2 \
  "$POLAR_FA" \
  "$APP_FA" \
  > "$PREP/Apennine_to_polar.paf"
```

Here:

- the polar-bear genome is the target coordinate system;
- Apennine Scaffold 25 is the query;
- `asm20` uses parameters suitable for relatively divergent assemblies;
- `-t 2` uses two CPU threads.

Check the alignment:

```bash
head "$PREP/Apennine_to_polar.paf"
```

The PAF file describes matching sequence segments, but it cannot yet be used by `liftOver`.

### 1.3 — Convert the alignment to a chain file

```bash
transanno minimap2chain \
  "$PREP/Apennine_to_polar.paf" \
  --output "$PREP/Apennine_to_polar.chain"
```

A chain file records corresponding blocks between two genome assemblies, including their coordinates, orientations and gaps.

Check the file:

```bash
grep '^chain[[:space:]]' \
  "$PREP/Apennine_to_polar.chain" | head
```

The chain must translate **Apennine coordinates into polar-bear coordinates**.

> Steps 1–5 are instructor preparation. Students may begin at Step 6 using the supplied Apennine-coordinate GERP file.

## Step 2 — Convert the VCF positions to BED

**Purpose:** represent every VCF SNP as a one-base interval that can be processed by `liftOver`.

```bash
bcftools query -f '%CHROM\t%POS\n' "$VCF" |
  awk -v chrom="$CHROM" '
    BEGIN {
      OFS="\t"
    }
    $1==chrom {
      print $1,$2-1,$2,$1":"$2
    }
  ' |
  sort -k1,1 -k2,2n -u \
  > "$PREP/Apennine_sites.bed"
```

VCF positions are one-based, while BED starts are zero-based. Therefore, a VCF position of `100` becomes the BED interval `99–100`.

The fourth column stores the original Apennine position as a stable site identifier.

Check the file:

```bash
head "$PREP/Apennine_sites.bed"
wc -l "$PREP/Apennine_sites.bed"
```

## Step 3 — Lift the SNPs to polar-bear coordinates

### 3.1 — Run liftOver

```bash
liftOver -multiple \
  "$PREP/Apennine_sites.bed" \
  "$PREP/Apennine_to_polar.chain" \
  "$PREP/polar_sites_all.bed" \
  "$PREP/unmapped.bed"
```

This produces:

- `polar_sites_all.bed`: SNPs mapped to polar-bear coordinates;
- `unmapped.bed`: SNPs that could not be mapped.

The `-multiple` option reports alternative mappings. Multiple mappings may occur in repetitive or duplicated regions.

### 3.2 — Keep only unique one-base mappings

```bash
awk '
  BEGIN {
    OFS="\t"
  }
  {
    count[$4]++
    line[$4]=$0
  }
  END {
    for (id in count)
      if (count[id]==1)
        print line[id]
  }
' "$PREP/polar_sites_all.bed" |
  awk 'BEGIN{OFS="\t"} $3-$2==1 {print}' |
  sort -k1,1 -k2,2n \
  > "$PREP/polar_sites_unique.bed"
```

This removes:

- sites mapping to more than one polar-bear position;
- mappings that are no longer exactly one base long.

Check the result:

```bash
printf 'Input sites: '
wc -l < "$PREP/Apennine_sites.bed"

printf 'Unique mapped sites: '
wc -l < "$PREP/polar_sites_unique.bed"

printf 'Duplicate IDs remaining: '
cut -f4 "$PREP/polar_sites_unique.bed" |
  sort |
  uniq -d |
  wc -l
```

The final check should print `0`.

## Step 4 — Extract GERP scores

**Purpose:** retrieve scores only from the polar-bear regions containing mapped SNPs.

### 4.1 — Identify the required polar-bear regions

```bash
awk '
  BEGIN {
    OFS="\t"
  }
  {
    if (!($1 in minimum) || $2<minimum[$1])
      minimum[$1]=$2

    if (!($1 in maximum) || $3>maximum[$1])
      maximum[$1]=$3
  }
  END {
    for (chromosome in minimum)
      print chromosome,minimum[chromosome],maximum[chromosome]
  }
' "$PREP/polar_sites_unique.bed" \
  > "$PREP/polar_regions.tsv"
```

```bash
cat "$PREP/polar_regions.tsv"
```

This creates one extraction interval for each polar-bear sequence reached by the SNPs.

### 4.2 — Extract the scores from the bigWig

```bash
while IFS=$'\t' read -r chromosome start end; do
  bigWigToBedGraph \
    -chrom="$chromosome" \
    -start="$start" \
    -end="$end" \
    "$GERP_BW" \
    "$PREP/GERP_${chromosome}.bedGraph"
done < "$PREP/polar_regions.tsv"
```

Combine the extracted score intervals:

```bash
cat "$PREP"/GERP_*.bedGraph |
  sort -k1,1 -k2,2n \
  > "$PREP/GERP_polar_regions.bedGraph"
```

A bedGraph row contains:

```text
chromosome    start    end    GERP_score
```

### 4.3 — Match SNPs with scores

```bash
bedtools intersect \
  -a "$PREP/polar_sites_unique.bed" \
  -b "$PREP/GERP_polar_regions.bedGraph" \
  -wa -wb \
  > "$PREP/sites_with_polar_scores.tsv"
```

Inspect the result:

```bash
head "$PREP/sites_with_polar_scores.tsv"
```

The first four columns describe the lifted SNP. The last four columns describe the matching GERP interval. The GERP score is in the final column.

## Step 5 — Return the scores to Apennine coordinates

**Purpose:** use the stored site IDs to create a GERP track in the original Scaffold 25 coordinates.

```bash
awk '
  BEGIN {
    OFS="\t"
  }
  $9 ~ /^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/ {
    count[$4]++
    score[$4]=$9
  }
  END {
    for (id in count) {
      if (count[id]==1) {
        split(id,position,":")
        print position[1],position[2]-1,position[2],score[id]
      }
    }
  }
' "$PREP/sites_with_polar_scores.tsv" |
  sort -k1,1 -k2,2n \
  > "$PREP/GERP_on_Apennine_unique.bed"
```

Compress and index the final track:

```bash
GERP="$OUTDIR/GERP_on_Apennine_unique.bed.gz"

bgzip -c "$PREP/GERP_on_Apennine_unique.bed" > "$GERP"
tabix -f -p bed "$GERP"
```

Check the result:

```bash
zcat "$GERP" | head

zcat "$GERP" |
  awk '
    $1!="Scaffold_25" || $3-$2!=1 {
      invalid++
    }
    END {
      print "Invalid rows:",invalid+0
    }
  '
```

The last command should report:

```text
Invalid rows: 0
```

The instructor provides this prepared file, so students do not need to download the approximately 7 GB bigWig or repeat Steps 1–5.

## Step 6 — Join GERP scores with GenoLoader sites

If starting from the supplied GERP track, define:

```bash
OUTDIR="$COURSE_DIR/results/gerp"
GERP="$COURSE_DIR/data/GERP_on_Apennine_unique.bed.gz"
GT="$COURSE_DIR/results/snpeff/Bears_4pops_s25.ann.POP_OUT.gt"

mkdir -p "$OUTDIR"
```

### 6.1 — Convert the GenoLoader table to BED

```bash
awk -F'\t' '
  BEGIN {
    OFS="\t"
  }
  NR>1 && $2~/^[0-9]+$/ {
    print $1,$2-1,$2,$3,$4,$5
  }
' "$GT" > "$OUTDIR/genoloader_sites.bed"
```

### 6.2 — Add the GERP scores

```bash
bedtools intersect \
  -a "$OUTDIR/genoloader_sites.bed" \
  -b "$GERP" \
  -wa -wb \
  > "$OUTDIR/sites_with_gerp.tsv"
```

Check the joined table:

```bash
wc -l \
  "$OUTDIR/genoloader_sites.bed" \
  "$OUTDIR/sites_with_gerp.tsv"

head "$OUTDIR/sites_with_gerp.tsv"
```

The final column contains the GERP score.

Not every GenoLoader site will necessarily have a GERP score. Missing scores should not be interpreted as zero.

## Step 7 — Plot the results

Define the focal-population sample lists:

```bash
ABB_LIST="$COURSE_DIR/data/ABB.samples"
SBB_LIST="$COURSE_DIR/data/SBB.samples"
```

Run the plotting script:

```bash
python scripts/plot_gerp_genotypes.py \
  "$GERP" \
  "$GT" \
  "$ABB_LIST" \
  "$SBB_LIST" \
  "$OUTDIR/GERP_derived_sites_GERP_gt2_by_sample.tsv" \
  "$OUTDIR/GERP_ABB_SBB_Scaffold_25"
```

The script produces:

```text
GERP_derived_sites_GERP_gt2_by_sample.tsv
GERP_ABB_SBB_Scaffold_25_score_distribution.pdf
GERP_ABB_SBB_Scaffold_25_derived_sites_GERP_gt2.pdf
```

Check the files:

```bash
ls -lh \
  "$OUTDIR/GERP_derived_sites_GERP_gt2_by_sample.tsv" \
  "$OUTDIR/GERP_ABB_SBB_Scaffold_25_score_distribution.pdf" \
  "$OUTDIR/GERP_ABB_SBB_Scaffold_25_derived_sites_GERP_gt2.pdf"
```

Inspect the table:

```bash
column -t \
  "$OUTDIR/GERP_derived_sites_GERP_gt2_by_sample.tsv" |
  head
```

## Understanding the plots

The score-distribution figure shows:

- the overall distribution of GERP scores;
- the threshold at GERP \(>2\);
- the distributions for the SnpEff impact classes.

The derived-site figure counts, for every bear:

- heterozygous derived sites with GERP \(>2\);
- homozygous-derived sites with GERP \(>2\);
- the total number of these derived sites.

A homozygous-derived genotype counts as **one site** in this figure, not two derived copies.

## Important interpretation

GERP \(>2\) identifies positions showing evidence of evolutionary constraint. It does not automatically identify deleterious variants.

The results should be interpreted cautiously because:

- the threshold of 2 is an analytical choice;
- some variants cannot be mapped between assemblies;
- some mapped sites do not have a GERP score;
- ancestral-state polarization can be incorrect;
- the analysis uses only Scaffold 25;
- GERP does not include dominance, gene expression or measured fitness effects.

### Questions for discussion

1. Why must the GERP bigWig and polar-bear FASTA use the same sequence names?
2. Why do we remove sites mapping to multiple locations?
3. Why is a missing GERP score different from a score of zero?
4. What is the difference between SnpEff impact and GERP conservation?
5. Why might homozygous-derived variants at conserved positions be important in an inbred population?
