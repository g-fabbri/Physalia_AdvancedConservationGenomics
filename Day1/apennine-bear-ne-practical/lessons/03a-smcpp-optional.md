Optional Part 3 — SMC++ with all ABB and SBB individuals

This extension is **outside the 90-minute practical**. MSMC2 uses one diploid individual per population; SMC++ can incorporate the other individuals without phasing. Both use sequence patterns to infer a size trajectory, so this is a comparison within the coalescent-method family, not an independent LD validation of GONE2. We still use only **Scaffold_25**, and the resulting curves are demonstrations rather than genome-wide estimates.

## Start your terminal

From the Day 1 course directory, activate the optional SMC++ environment:

```bash
conda activate bear-ne-smcpp
```

This environment is separate from `bear-ne-practical` because the Bioconda SMC++ package requires a different Python version. It provides `smc++`, `bcftools`, and `bedtools`. Run the activation command whenever you open a new terminal.

**Check:**

```bash
command -v smc++
command -v bcftools
command -v bedtools
```

Each command should print an executable path inside the `bear-ne-smcpp` environment. If `smc++` is missing, follow the [software setup](../software/README.md) before continuing.

This lesson uses **Scaffold_25** throughout, matching the MSMC2 and GONE2 classroom analyses.

## Step 1 — Identify samples and sequence length

**Purpose:** confirm that the 18-individual VCF contains the 10 ABB and 8 SBB samples and obtain the reference length of Scaffold_25.

```bash
CHROM=Scaffold_25
VCF=data/UrArMa_18i_s25.vcf.gz
CALLABLE=data/UrArMa_callable.bed.gz
OUTDIR=results/smcpp
mkdir -p "$OUTDIR"

ABB_SAMPLES=$(paste -sd, data/apennine.samples)
SBB_SAMPLES=$(paste -sd, data/slovak.samples)
printf 'ABB: %s\nSBB: %s\n' "$ABB_SAMPLES" "$SBB_SAMPLES"
bcftools query -l "$VCF" | wc -l
bcftools index -s "$VCF" | awk -v c="$CHROM" '$1==c && $2~/^[0-9]+$/ {print $1"\t"$2}' > "$OUTDIR/${CHROM}.genome"
cat "$OUTDIR/${CHROM}.genome"
```

**Expected:** 18 sample IDs and one `Scaffold_25` line with its reference length. If the genome file is empty, obtain the verified Scaffold_25 length from the reference FASTA index (`.fai`) before continuing; do **not** substitute the position of the last variant for chromosome length.

## Step 2 — Convert callable intervals into an exclusion mask

**Purpose:** SMC++ `--mask` marks **uncallable** positions. The course BED lists *callable* positions, so we must take its complement across the full scaffold.

```bash
zcat "$CALLABLE" |
  awk -v c="$CHROM" '$1==c {print $1"\t"$2"\t"$3}' |
  sort -k1,1 -k2,2n |
  bedtools merge -i - > "$OUTDIR/${CHROM}.callable.bed"

bedtools complement \
  -i "$OUTDIR/${CHROM}.callable.bed" \
  -g "$OUTDIR/${CHROM}.genome" \
  > "$OUTDIR/${CHROM}.uncallable.bed"

head "$OUTDIR/${CHROM}.uncallable.bed"
```

**Expected:** BED intervals on Scaffold_25 that were **not** in the callable file. Do not give `UrArMa_callable.bed.gz` directly to `smc++ --mask`: that would hide the sequence we want to analyse. Confirm that this shared callable mask is appropriate for all 18 bears; common alignment to one reference alone does not establish equal callability.

## Step 3 — Convert each population

**Purpose:** use **all** individuals in each population while choosing one unphased diploid representative as SMC++'s distinguished pair. The remaining individuals contribute their allele counts. Do not put ABB and SBB into one population label.

```bash
smc++ vcf2smc \
  -d 4573 4573 \
  --mask "$OUTDIR/${CHROM}.uncallable.bed" \
  "$VCF" "$OUTDIR/ABB_${CHROM}.smc.gz" "$CHROM" "ABB:$ABB_SAMPLES"

smc++ vcf2smc \
  -d U1916 U1916 \
  --mask "$OUTDIR/${CHROM}.uncallable.bed" \
  "$VCF" "$OUTDIR/SBB_${CHROM}.smc.gz" "$CHROM" "SBB:$SBB_SAMPLES"

ls -lh "$OUTDIR"/*.smc.gz
```

**Expected:** two non-empty `.smc.gz` files. With unphased genotypes, SMC++ recommends specifying the **same individual twice** for `-d`. The joint VCF can supply both conversions; the population sample lists select which genotypes contribute to each result.

## Step 4 — Estimate trajectories

**Purpose:** fit each population separately using the same mutation rate as the MSMC2 lesson.

```bash
MU=1.82e-8
smc++ estimate -o "$OUTDIR/ABB" "$MU" "$OUTDIR/ABB_${CHROM}.smc.gz"
smc++ estimate -o "$OUTDIR/SBB" "$MU" "$OUTDIR/SBB_${CHROM}.smc.gz"
ls -lh "$OUTDIR"/{ABB,SBB}/model.final.json
```

**Expected:** one `model.final.json` per population. Fitting may take much longer than the classroom slot. SMC++ folds the frequency-spectrum information by default; do not use `--unfold` unless ancestral alleles have been independently established.

## Step 5 — Plot and discuss

```bash
smc++ plot -g 11 --logy \
  "$OUTDIR/SMCPP_ABB_SBB_${CHROM}.pdf" \
  "$OUTDIR/ABB/model.final.json" \
  "$OUTDIR/SBB/model.final.json"
```

**Expected:** a comparison PDF in `results/smcpp/`. Inspect the fit diagnostics before interpreting it. Ask: What changes when the same broad method uses 10 or 8 individuals rather than one? Which differences might instead come from the distinguished bear, the shared mask, or the single scaffold? SMC++ assumes sites absent from a variant-only VCF are homozygous ancestral unless masked, so mask construction is central to the interpretation. The classroom VCF was filtered to exclude sites with missing genotypes; if a genuinely variable, otherwise-callable site was removed by that filter, SMC++ would misread its absence as reference homozygosity. A full analysis needs input and callability definitions designed together.

Return to the [synthesis](04-synthesis.md).
