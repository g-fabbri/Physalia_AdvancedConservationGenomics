# Part 0 — Meet the four-species dataset

Estimated terminal time: 10 minutes.

ABB and SBB are the **focal populations**. BLB and POB are **outgroups** used to infer ancestral states; they are not added to either focal population. All four groups were called against the Apennine reference assembly, so REF is the assembly allele, not necessarily the ancestral allele.

## Step 1 — Name the inputs

**Purpose:** use a single VCF and four explicit sample lists throughout the lesson. Replace the VCF placeholder with the instructor-provided path.

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
ls -lh "$VCF" "$ABB_LIST" "$SBB_LIST" "$BLB_LIST" "$POB_LIST"
bcftools query -l "$VCF" | wc -l
for LIST in "$ABB_LIST" "$SBB_LIST" "$BLB_LIST" "$POB_LIST"; do
  printf '%s: ' "$LIST"
  wc -l < "$LIST"
done
~~~

Every file must exist; the four list counts should sum to the expected number of VCF samples if no samples are intentionally unassigned. A list entry absent from the VCF must be resolved before polarization.

## Step 2 — Inspect variant representation

**Purpose:** confirm that the teaching input is a SNP VCF, not a collection of all callable sites. This distinction matters when interpreting per-individual counts.

~~~bash
bcftools view -H "$VCF" | head -n 3
bcftools index -n "$VCF"
~~~

**Expected:** variant rows and a record count. The count is **not** the number of callable bases. If a full-genome file is too large for class, the instructor should distribute a small, fixed autosomal subset and separately prepared genome-wide summaries.

**Check:** chromosome names in the VCF must match the SnpEff database. Continue to [SnpEff](01-snpeff.md).
