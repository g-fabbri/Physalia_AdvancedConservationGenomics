# Part 0 — Meet the four-species dataset

Estimated terminal time: 10 minutes.

ABB and SBB are the **focal populations**. BLB and POB (polar bears) are **outgroups** used to infer ancestral states; they are not added to either focal population. All four groups were called against the Apennine reference assembly, so REF is the assembly allele, not necessarily the ancestral allele. This practical examines **Scaffold_34 only**.

## Step 1 — Name the inputs

**Purpose:** use the prepared Scaffold_34 VCF and four actual sample lists throughout the lesson, without extracting or re-filtering the full genome during class.

**Input:** indexed VCF and plain-text sample lists, one ID per line.

~~~bash
VCF=data/Bears_4pops_s34.vcf.gz
ABB_LIST=data/ABB.samples
SBB_LIST=data/SBB.samples
BLB_LIST=data/BLB.samples
POB_LIST=data/POB.samples
OUTDIR=results/genetic_load
mkdir -p "$OUTDIR"
~~~

**Expected:** no output. These variables exist only in the current terminal.

**Check:**

~~~bash
ls -lh "$VCF" "$VCF.csi" "$ABB_LIST" "$SBB_LIST" "$BLB_LIST" "$POB_LIST"
bcftools query -l "$VCF" | wc -l
for LIST in "$ABB_LIST" "$SBB_LIST" "$BLB_LIST" "$POB_LIST"; do
  printf '%s: ' "$LIST"
  wc -l < "$LIST"
done
~~~

Every file must exist; the four list counts should sum to the expected number of VCF samples if no samples are intentionally unassigned. A list entry absent from the VCF must be resolved before polarization.

## Step 2 — Inspect variant representation

**Purpose:** verify the scaffold and variant representation before interpreting per-individual counts. The filename suggests Scaffold_34, but we check its contents rather than assume.

~~~bash
bcftools view -H "$VCF" | head -n 3
bcftools index -n "$VCF"
bcftools query -f '%CHROM\n' "$VCF" | sort -u
~~~

**Expected:** variant rows, a record count, and exactly one chromosome name: **Scaffold_34**. The record count is **not** the number of callable bases. Check with the instructor that the input has been restricted to biallelic SNPs before GenoLoader; do not silently discard records during the lesson.

**Check:** chromosome names in the VCF must match the SnpEff database. Because this is one scaffold, differences in counts cannot establish a genome-wide difference in genetic load. Continue to [SnpEff](01-snpeff.md).
