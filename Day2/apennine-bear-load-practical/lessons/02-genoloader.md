# Part 2 — Which alleles are derived, and who carries them?

In this lesson, we will use **GenoLoader** (https://github.com/emitruc/genoloader) to:

1. infer the ancestral and derived allele at each SNP;
2. count how many derived alleles each bear carries;
3. compare the ABB and SBB populations.

The SnpEff-annotated VCF produced in Part 1 is used as input. Black bears and polar bears are used as outgroups to help identify the ancestral allele.

## Start your terminal

From your Day 2 directory, run:

```bash
conda activate bear-load-practical

COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$COURSE_DIR/software/genoloader:$PATH"
```

Define the input files and output directory:

```bash
ANNOTATED="$COURSE_DIR/results/snpeff/Bears_4pops_s25.ann.vcf"

ABB_LIST="$COURSE_DIR/data/ABB.samples"
SBB_LIST="$COURSE_DIR/data/SBB.samples"
BLB_LIST="$COURSE_DIR/data/BLB.samples"
POB_LIST="$COURSE_DIR/data/POB.samples"

GENOLOADER="$COURSE_DIR/software/genoloader/genoloader"

OUTDIR="$COURSE_DIR/results/genoloader"
mkdir -p "$OUTDIR"
```

Here:

- `ANNOTATED` is the VCF annotated with SnpEff;
- `ABB_LIST` and `SBB_LIST` contain the focal brown-bear samples;
- `BLB_LIST` and `POB_LIST` contain the black- and polar-bear outgroups;
- `GENOLOADER` is the program used for polarization;
- `OUTDIR` stores the summary tables and figures.

Check the input files:

```bash
ls -lh \
  "$ANNOTATED" \
  "$ABB_LIST" \
  "$SBB_LIST" \
  "$BLB_LIST" \
  "$POB_LIST"

test -x "$GENOLOADER" && echo "GenoLoader is ready"
```

## Step 1 — Prepare the outgroup samples

**Purpose:** combine black and polar bears into one outgroup list.

```bash
cat "$BLB_LIST" "$POB_LIST" | \
  sort -u > "$OUTDIR/outgroups.samples"
```

`cat` joins the two sample lists. `sort -u` sorts the IDs and removes any duplicate names.

Count the samples in each group:

```bash
N_ABB=$(wc -l < "$ABB_LIST")
N_SBB=$(wc -l < "$SBB_LIST")
N_OUT=$(wc -l < "$OUTDIR/outgroups.samples")

printf 'ABB=%s SBB=%s outgroups=%s\n' \
  "$N_ABB" "$N_SBB" "$N_OUT"
```

The expected focal-population counts are:

```text
ABB=10
SBB=8
```

Display the combined outgroup list:

```bash
cat "$OUTDIR/outgroups.samples"
```

### Why do we need outgroups?

The VCF `REF` allele is simply the allele present in the Apennine reference genome. It is not necessarily the ancestral allele.

GenoLoader uses the outgroup genotypes to determine which allele is likely ancestral. The other allele is then considered derived.

Using two outgroup species provides additional evidence, but polarization can still be affected by:

- missing or incorrect genotypes;
- ancestral polymorphism;
- mutations occurring independently in different lineages;
- introgression between species.

## Step 2 — Run GenoLoader

**Purpose:** polarize the SNPs and recode every genotype according to the number of derived alleles.

```bash
"$GENOLOADER" "$ANNOTATED" \
  --p1 "$ABB_LIST" \
  --p2 "$SBB_LIST" \
  --p0 "$OUTDIR/outgroups.samples" \
  --m1 "$N_ABB" \
  --m2 "$N_SBB" \
  --m0 "$N_OUT" \
  --polX POP_OUT \
  --low_cov NO
```

The main arguments are:

| Argument | Meaning |
|---|---|
| `--p1` | Sample list for ABB |
| `--p2` | Sample list for SBB |
| `--p0` | Combined outgroup sample list |
| `--m1`, `--m2`, `--m0` | Minimum numbers of called individuals required in each group |
| `--polX POP_OUT` | Infer the ancestral allele using the outgroup population |
| `--low_cov NO` | Use the diploid genotype calls from the VCF |

GenoLoader reports:

- the samples assigned to each population;
- the number of loci processed;
- the number of loci retained;
- how many loci were reoriented because the VCF `ALT` allele was inferred to be ancestral;
- how many loci failed the missing-data requirements.

## Step 3 — Find and inspect the GenoLoader output

GenoLoader writes its output beside the input VCF. It removes the final `.vcf` extension and adds `.POP_OUT.gt`.

Define the resulting filename:

```bash
GT="${ANNOTATED%.vcf}.POP_OUT.gt"
```

This corresponds to:

```text
results/snpeff/Bears_4pops_s25.ann.POP_OUT.gt
```

Check the file:

```bash
ls -lh "$GT"
wc -l "$GT"
head -n 3 "$GT"
```

The `.gt` file is a tab-separated table. It contains:

| Column | Meaning |
|---|---|
| Scaffold and position | Genomic location of the SNP |
| Effect | SnpEff impact category |
| Variant type | Simplified functional annotation |
| Flag | How the ancestral allele was inferred |
| Reference | Allele inferred by GenoLoader to be ancestral |
| Sample columns | Number of derived alleles carried by each bear |

The sample dosages are:

- `0`: ancestral homozygote;
- `1`: heterozygote;
- `2`: derived homozygote;
- `nan`: missing genotype.

Count the polarization flags:

```bash
awk -F'\t' '
  NR>1 {
    count[$5]++
  }
  END {
    for (flag in count)
      print flag, count[flag]
  }
' "$GT" | sort
```

These flags help distinguish confidently polarized sites from ambiguous or fallback cases.

## Step 4 — Summarize derived alleles

**Purpose:** calculate derived-allele summaries for every ABB and SBB individual.

```bash
python scripts/summarize_genoloader.py \
  "$GT" \
  "$ABB_LIST" \
  "$SBB_LIST" \
  "$OUTDIR/derived_burden_by_sample.tsv" \
  "$OUTDIR/derived_frequency_by_impact.tsv" \
  "$OUTDIR/GenoLoader_ABB_SBB_Scaffold_25"
```

The script reads the polarized `.gt` table and calculates:

- the number of heterozygous derived sites;
- the number of homozygous-derived sites;
- the total number of derived allele copies;
- values normalized by the number of called sites;
- mean derived-allele frequencies in ABB and SBB;
- directional Rxy values.


List the generated figures:

```bash
ls -lh "$OUTDIR"/GenoLoader_ABB_SBB_Scaffold_25_*.pdf
```

## How are the counts calculated?

At each site:

- a heterozygous genotype contributes **one** derived copy;
- a derived homozygous genotype contributes **two** derived copies.

Therefore:

```text
total derived copies =
heterozygous sites + 2 × homozygous-derived sites
```

The plots compare these values between ABB and SBB for the SnpEff impact categories:

- `HIGH`;
- `MODERATE`;
- `LOW`;
- `MODIFIER`.

Each point represents one bear, while the horizontal bar represents the population mean.

## Important interpretation

Remember that:

- SnpEff effects are computational predictions;
- derived does not automatically mean deleterious;
- the analysis uses only Scaffold 25;
- the sample contains 10 ABB and 8 SBB bears;
- polarization errors can affect derived-allele counts.

## Questions for discussion

1. Why is the VCF `REF` allele not necessarily the ancestral allele?
2. Why might black and polar bears disagree about the ancestral allele?
3. What is the difference between a heterozygous derived site and a homozygous-derived site?
4. Why could homozygous-derived variants be especially relevant for recessive deleterious effects?
