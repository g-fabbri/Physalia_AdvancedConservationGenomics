# Part 0 — Meet the four-species dataset

Estimated terminal time: 10 minutes.

## Start your terminal

From the Day 2 directory containing `data/`, `software/`, and `results/`, run this in each new terminal:

```bash
conda activate bear-load-practical
COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$COURSE_DIR/software/genoloader:$PATH"
```

Conda supplies the shared tools; standalone programs are kept under Day 2 `software/`. The instructor places command-line executables needed for GERP in `software/bin/`. Check your location with `pwd` before continuing.

ABB and SBB are the **focal populations**. BLB and POB (polar bears) are **outgroups** used to infer ancestral states; they are not added to either focal population. All four groups were called against the Apennine reference assembly, so REF is the assembly allele, not necessarily the ancestral allele. This practical examines **Scaffold_25 only**.

## Step 1 — Name the inputs

**Purpose:** use the prepared Scaffold_25 VCF and four actual sample lists throughout the lesson, without extracting or re-filtering the full genome during class.

**Input:** indexed VCF and plain-text sample lists, one ID per line.

### 1.1 — Set the input paths

**Purpose:** give the prepared VCF, sample lists, and results directory names that later commands can reuse.

~~~bash
VCF=data/Bears_4pops_s25.vcf.gz
ABB_LIST=data/ABB.samples
SBB_LIST=data/SBB.samples
BLB_LIST=data/BLB.samples
POB_LIST=data/POB.samples
OUTDIR=results/genetic_load
mkdir -p "$OUTDIR"
~~~

**Expected:** no output. These variables exist only in the current terminal.

### 1.2 — Confirm the files and sample counts

**Purpose:** catch missing inputs or incomplete sample lists before analysis.

~~~bash
ls -lh "$VCF" "$VCF.csi" "$ABB_LIST" "$SBB_LIST" "$BLB_LIST" "$POB_LIST"
~~~

~~~bash
bcftools query -l "$VCF" | wc -l
~~~

~~~bash
for LIST in "$ABB_LIST" "$SBB_LIST" "$BLB_LIST" "$POB_LIST"; do
  printf '%s: ' "$LIST"
  wc -l < "$LIST"
done
~~~

Every file must exist; the four list counts should sum to the expected number of VCF samples if no samples are intentionally unassigned. A list entry absent from the VCF must be resolved before polarization.

## Step 2 — Inspect variant representation

**Purpose:** verify the scaffold and variant representation before interpreting per-individual counts. The filename suggests Scaffold_25, but we check its contents rather than assume.

### 2.1 — Inspect records and scaffold names

**Purpose:** confirm that the indexed VCF contains the intended scaffold and count its variant records.

~~~bash
bcftools view -H "$VCF" | head -n 3
bcftools index -n "$VCF"
bcftools query -f '%CHROM\n' "$VCF" | sort -u
~~~

**Expected:** variant rows, a record count, and exactly one chromosome name: **Scaffold_25**. The record count is **not** the number of callable bases. Check with the instructor that the input has been restricted to biallelic SNPs before GenoLoader; do not silently discard records during the lesson.

**Check:** chromosome names in the VCF must match the SnpEff database. Because this is one scaffold, differences in counts cannot establish a genome-wide difference in genetic load. Continue to [SnpEff](01-snpeff.md).
