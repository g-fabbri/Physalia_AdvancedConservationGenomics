# Optional Part 3B — currentNe2 with a full autosomal VCF

This extension is **outside the timed practical** and requires a full autosomal VCF that is not included in the one-scaffold classroom dataset. GONE2 estimates a recent trajectory from LD at different recombination distances. currentNe2 also uses LD, but targets a **contemporary** effective population size. Agreement between them is therefore not independent confirmation.

This is an example that students or instructors can try later when a suitable genome-wide VCF is available. Do not run currentNe2 on Scaffold_25 alone: currentNe2 recognises chromosome information only when at least two chromosomes are present, and a one-scaffold estimate would not represent genome-wide contemporary Ne.

## Start your terminal

From the Day 1 course directory, activate the main environment and identify the executable:

~~~bash
conda activate bear-ne-practical
CURRENTNE2_BIN="$PWD/software/currentNe2/currentne2"
~~~

**Check:**

~~~bash
command -v bcftools
test -x "$CURRENTNE2_BIN" && echo "currentNe2 ready: $CURRENTNE2_BIN"
~~~

## Step 1 — Define the full autosomal VCF

**Purpose:** select a multisample VCF containing the autosomes, physical positions, and all ABB and SBB individuals.

~~~bash
FULL_VCF=data/UrArMa_18i_autosomes.vcf.gz
WORKDIR=results/currentne2/input
OUTDIR=results/currentne2
mkdir -p "$WORKDIR" "$OUTDIR"
~~~

The filename is an example. Replace it with the path to the prepared full autosomal VCF. It should already contain high-quality biallelic SNPs and exclude sex chromosomes, unplaced sequence, and loci that failed the project-level genotype filters.

**Check:**

~~~bash
bcftools query -l "$FULL_VCF" | wc -l
bcftools index -s "$FULL_VCF" | head
bcftools index -n "$FULL_VCF"
~~~

**Expected:** 18 samples, multiple autosomes, and a positive record count. CHROM and POS in the VCF provide the map information used by currentNe2.

## Step 2 — Create one VCF per population

**Purpose:** estimate ABB and SBB separately while retaining the same autosomal regions and variant-processing rules.

~~~bash
bcftools view \
  -S data/apennine.samples \
  -m2 -M2 -v snps \
  -Oz -o "$WORKDIR/ABB_autosomes.vcf.gz" \
  "$FULL_VCF"

bcftools view \
  -S data/slovak.samples \
  -m2 -M2 -v snps \
  -Oz -o "$WORKDIR/SBB_autosomes.vcf.gz" \
  "$FULL_VCF"

bcftools index -f "$WORKDIR/ABB_autosomes.vcf.gz"
bcftools index -f "$WORKDIR/SBB_autosomes.vcf.gz"
~~~

The sample list selects individuals, **-v snps** retains SNPs, and **-m2 -M2** retains biallelic records. Some sites can become monomorphic after separating populations; currentNe2 removes non-polymorphic sites during preprocessing.

**Check:**

~~~bash
for POPULATION in ABB SBB; do
  printf '%s samples: ' "$POPULATION"
  bcftools query -l "$WORKDIR/${POPULATION}_autosomes.vcf.gz" | wc -l

  printf '%s chromosomes: ' "$POPULATION"
  bcftools query -f '%CHROM\n' "$WORKDIR/${POPULATION}_autosomes.vcf.gz" |
    sort -u | wc -l

  printf '%s records: ' "$POPULATION"
  bcftools index -n "$WORKDIR/${POPULATION}_autosomes.vcf.gz"
done
~~~

**Expected:** 10 ABB samples, 8 SBB samples, more than one chromosome in both files, and positive record counts.

## Step 3 — Reduce very large inputs reproducibly

**Purpose:** remain below currentNe2's default compiled limit of 2,000,000 input loci and reduce memory and pairwise-comparison requirements.

If either population VCF approaches or exceeds two million records, thin both using the same rule:

~~~bash
for POPULATION in ABB SBB; do
  bcftools +prune \
    "$WORKDIR/${POPULATION}_autosomes.vcf.gz" \
    -Oz -o "$WORKDIR/${POPULATION}_autosomes_thinned.vcf.gz" \
    -- -n 1 -w 2kb -N rand --random-seed 1

  bcftools index -f "$WORKDIR/${POPULATION}_autosomes_thinned.vcf.gz"
done
~~~

This keeps at most one randomly selected SNP per 2 kb window with a fixed seed. It is a computational teaching choice, not a universally optimal biological filter. Use the same window and seed for both populations. The retained sites need not be identical because polymorphism differs between ABB and SBB.

**Check:**

~~~bash
for POPULATION in ABB SBB; do
  printf '%s thinned records: ' "$POPULATION"
  bcftools index -n "$WORKDIR/${POPULATION}_autosomes_thinned.vcf.gz"
done
~~~

Each count must be below 2,000,000. Increase the window if necessary. If the original VCFs are already comfortably below the limit, skip this step.

## Step 4 — Write plain VCF input

**Purpose:** create uncompressed VCF files with predictable filenames for currentNe2.

~~~bash
for POPULATION in ABB SBB; do
  bcftools view \
    -Ov -o "$WORKDIR/${POPULATION}_autosomes_currentNe2.vcf" \
    "$WORKDIR/${POPULATION}_autosomes_thinned.vcf.gz"
done
~~~

If Step 3 was skipped, replace the thinned input filename with **${POPULATION}_autosomes.vcf.gz**.

**Expected:** two nonempty text VCFs. currentNe2 reads chromosome assignments from CHROM, physical locations from POS, and population genotypes from the sample columns.

## Step 5 — Estimate contemporary Ne

**Purpose:** convert physical distance to genetic distance under the same approximation for both populations and estimate contemporary Ne.

~~~bash
for POPULATION in ABB SBB; do
  "$CURRENTNE2_BIN" \
    -r 1 \
    -t 2 \
    -o "$OUTDIR/${POPULATION}_currentNe2_OUTPUT.txt" \
    "$WORKDIR/${POPULATION}_autosomes_currentNe2.vcf"
done
~~~

| Option | Meaning |
|---|---|
| **-r 1** | Convert physical positions using a constant recombination rate of 1 cM/Mb. An empirical genetic map is preferable when available. |
| **-t 2** | Use two computational threads. |
| **-o** | Write each population's report to an explicit output filename. |

No genome-size argument is supplied: the full VCF provides multiple chromosome assignments and physical positions. We do not use **-x**, which fits a two-subpopulation metapopulation model and represents a different biological hypothesis.

**Expected:** two reports with a chromosome count, positive genome size, retained SNP count, diagnostics, and contemporary Ne estimates. The run may require substantial RAM and time because the number of possible SNP pairs grows rapidly with marker count.

## Step 6 — Inspect the reports

~~~bash
cat "$OUTDIR/ABB_currentNe2_OUTPUT.txt"
cat "$OUTDIR/SBB_currentNe2_OUTPUT.txt"
~~~

First inspect:

- total and effective sample sizes;
- input and retained polymorphic SNP counts;
- chromosome count and inferred genome size;
- missing-data proportion;
- warnings or convergence messages.

Then examine observed d², expected and observed heterozygosity, the F statistic, inferred full-sibling pairs, and the reported Ne estimates. Candidate sibling pairs are model-based inferences, not verified pedigrees. A positive F indicates a homozygote excess; a negative F indicates a heterozygote excess. Either may reflect biology, structure, relatives, sampling, or filtering and should be investigated.

Compare ABB and SBB only after confirming comparable autosomal coverage, filtering, thinning, recombination assumptions, and sample definitions. Ask whether currentNe2 and the recent portion of GONE2 point in the same direction, and whether shared LD assumptions could explain agreement or disagreement.

Even with a full VCF, Ne is not census size. Small samples, relatives, population structure, uneven recombination, genotyping error, missingness, and SNP ascertainment can strongly influence a contemporary LD estimate.


