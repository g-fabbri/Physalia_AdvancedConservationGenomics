# Optional Part 3 — SMC++ with all ABB and SBB individuals

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
command -v bgzip
command -v tabix
```

Each command should print an executable path inside the `bear-ne-smcpp` environment. `bgzip` and `tabix` are needed because SMC++ requires an indexed mask. If `smc++` is missing, follow the [software setup](../software/README.md) before continuing.

This lesson uses **Scaffold_25** throughout, matching the MSMC2 and GONE2 classroom analyses.

The workflow has four transformations:

```text
population sample lists + joint VCF + callable regions
                         ↓
        population-specific SMC++ input files
                         ↓
             fitted demographic models
                         ↓
              comparative Ne figure
```

Unlike GONE2, SMC++ does not estimate Ne directly from pairwise LD bins. It combines the sequence of homozygous and heterozygous regions in a distinguished diploid individual with allele-frequency information from the remaining individuals.

## Step 1 — Identify samples and sequence length

**Purpose:** confirm that the 18-individual VCF contains the 10 ABB and 8 SBB samples and obtain the reference length of Scaffold_25.

```bash
CHROM=Scaffold_25
VCF=data/UrArMa_18i_s25.vcf.gz
CALLABLE=data/UrArMa_callable_s25.bed.gz
OUTDIR=results/smcpp
mkdir -p "$OUTDIR"

ABB_SAMPLES=$(paste -sd, data/apennine.samples)
SBB_SAMPLES=$(paste -sd, data/slovak.samples)
printf 'ABB: %s\nSBB: %s\n' "$ABB_SAMPLES" "$SBB_SAMPLES"
bcftools query -l "$VCF" | wc -l

bcftools index -s "$VCF" | awk -v c="$CHROM" '$1==c && $2~/^[0-9]+$/ {print $1"\t"$2}' > "$OUTDIR/${CHROM}.genome"
cat "$OUTDIR/${CHROM}.genome"
```

`paste -sd,` joins each population's sample IDs with commas, which is the format expected later by `vcf2smc`. `bcftools query -l` counts the samples in the joint VCF. The final command reads the VCF index and writes a two-column genome file containing the scaffold name and its complete reference length.

**Expected:** the population lines contain 10 comma-separated ABB IDs and 8 comma-separated SBB IDs; the count is `18`; and the genome file contains one line shaped like:

```text
Scaffold_25    <complete scaffold length in bp>
```

The second value is the full reference length, not the number of variants or the position of the final SNP. `bedtools complement` needs this boundary to identify uncallable sequence before the first callable interval and after the last one. If the genome file is empty, obtain the verified Scaffold_25 length from the reference FASTA index (`.fai`) before continuing.

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
  | bgzip -c > "$OUTDIR/${CHROM}.uncallable.bed.gz"

tabix -f -p bed "$OUTDIR/${CHROM}.uncallable.bed.gz"

zcat "$OUTDIR/${CHROM}.uncallable.bed.gz" | head
ls -lh "$OUTDIR/${CHROM}.uncallable.bed.gz" \
  "$OUTDIR/${CHROM}.uncallable.bed.gz.tbi"
```

`bedtools complement` creates intervals that are outside the callable regions. `bgzip` compresses the BED in a block-addressable format, and `tabix -p bed` creates the `.tbi` index required by SMC++. Ordinary `gzip` compression is not sufficient.

The intermediate files have different meanings:

| File | Meaning |
|---|---|
| **Scaffold_25.callable.bed** | Merged regions where sequence data are considered reliable. This is a positive mask. |
| **Scaffold_25.uncallable.bed.gz** | Complementary regions that SMC++ must ignore. This is the negative mask passed to `--mask`. |
| **Scaffold_25.uncallable.bed.gz.tbi** | Tabix index that lets SMC++ retrieve mask intervals for Scaffold_25 efficiently. |

**Expected:** the first command displays BED intervals on Scaffold_25 that were **not** in the callable file. The second check lists both the compressed mask and its Tabix index. Do not give `UrArMa_callable.bed.gz` directly to `smc++ --mask`: that would hide the sequence we want to analyse. Confirm that this shared callable mask is appropriate for all 18 bears; common alignment to one reference alone does not establish equal callability.

## Step 3 — Convert each population

**Purpose:** use **all** individuals in each population while choosing one unphased diploid representative as SMC++'s distinguished pair. The remaining individuals contribute their allele counts. Do not put ABB and SBB into one population label.

The distinguished individual contributes the sequential information used to model changes between hidden genealogical states along the scaffold. The other individuals are *undistinguished*: SMC++ does not trace their individual genealogies, but uses their derived/reference allele counts to strengthen the frequency-spectrum information. Here bear 4573 represents ABB and U1916 represents SBB.

```bash
smc++ vcf2smc \
  -d 4573 4573 \
  --mask "$OUTDIR/${CHROM}.uncallable.bed.gz" \
  "$VCF" "$OUTDIR/ABB_${CHROM}.smc.gz" "$CHROM" "ABB:$ABB_SAMPLES"

smc++ vcf2smc \
  -d U1916 U1916 \
  --mask "$OUTDIR/${CHROM}.uncallable.bed.gz" \
  "$VCF" "$OUTDIR/SBB_${CHROM}.smc.gz" "$CHROM" "SBB:$SBB_SAMPLES"

ls -lh "$OUTDIR"/*.smc.gz
```

Important arguments:

| Argument | Meaning |
|---|---|
| **-d 4573 4573** | Use both unphased chromosome copies of ABB bear 4573 as the distinguished pair. The SBB command does the same for U1916. |
| **--mask** | Exclude uncallable intervals rather than treating absence of variants as evidence of homozygosity. |
| **ABB:$ABB_SAMPLES** | Name the population ABB and select only the IDs stored in the ABB sample list. |
| **CHROM** | Restrict conversion to Scaffold_25. |

During conversion, SMC++ should report one population, two distinguished lineages, and the remaining population chromosomes as undistinguished lineages. For ABB, 4573 should appear under `Distinguished lineages`; the other nine ABB bears should appear under `Undistinguished lineages`. For SBB, these positions should be occupied by U1916 and the other seven SBB bears.

**Expected:** two non-empty `.smc.gz` files. These are compressed SMC++ observation files, not ordinary VCFs and not final Ne estimates. They encode runs of sequence with their callable length, distinguished genotype state, and allele-count information. The joint VCF can supply both conversions because the population expressions select the relevant genotypes. The `pkg_resources is deprecated` message comes from the packaged SMC++ version and is a warning, not the cause of a failed conversion.

**Check:** a traceback means conversion failed even if some informational lines were printed. Successful conversion returns to the prompt and leaves both files with nonzero sizes.

## Step 4 — Estimate trajectories

**Purpose:** fit each population separately using the same mutation rate as the MSMC2 lesson. Because this teaching dataset contains only one scaffold, we use a deliberately simple model over an explicit time interval rather than relying on SMC++ to infer very broad model limits automatically.

```bash
MU=1.82e-8

smc++ estimate \
  --timepoints 10 100000 \
  --knots 6 \
  --cores 2 \
  -o "$OUTDIR/ABB_bounded" \
  "$MU" "$OUTDIR/ABB_${CHROM}.smc.gz"

smc++ estimate \
  --timepoints 10 100000 \
  --knots 6 \
  --cores 2 \
  -o "$OUTDIR/SBB_bounded" \
  "$MU" "$OUTDIR/SBB_${CHROM}.smc.gz"

ls -lh "$OUTDIR"/{ABB_bounded,SBB_bounded}/model.final.json
```

| Option | Meaning |
|---|---|
| **--timepoints 10 100000** | Fit the model from 10 to 100,000 generations before present instead of using automatically selected limits. These are model-fitting bounds, not proof that the complete interval is well resolved. |
| **--knots 6** | Use a relatively small number of change points, reducing flexibility and numerical instability in a one-scaffold exercise. |
| **--cores 2** | Use two processor cores. |
| **-o** | Write each population to a separate, newly named output directory. |

**Expected:** one `model.final.json` per population. Messages about EM iterations, log likelihood, and the current model are normal. The repeated `pkg_resources is deprecated` message is a packaging warning and does not itself indicate failure.

The main output and progress messages mean:

| Output | Interpretation |
|---|---|
| **Loading data / Gb of data** | SMC++ successfully opened the `.smc.gz` file and reports the amount of sequence represented. |
| **theta and rho** | Internal scaled mutation and recombination quantities used to initialise the model; they are not direct Ne estimates. |
| **E-step / EM iteration** | SMC++ is alternating between calculating genealogical expectations and updating demographic parameters. |
| **Loglik** | Model log likelihood. It is mainly useful for monitoring optimisation; its absolute value should not be compared casually between datasets of different sizes. |
| **Current model / ASCII plot** | An intermediate trajectory printed during optimisation, not the final result. |
| **model.iter*.json** | Intermediate models that can reveal instability or overfitting across iterations. |
| **model.final.json** | Final fitted model used by `smc++ plot`. Its presence indicates that the estimation command completed. |

Do not interpret the numerical rows printed below `Current model` as a ready-to-use results table. SMC++ stores the fitted time and population-size parameters in JSON and performs the biological scaling during plotting using the supplied mutation rate and generation time.

If the command ends with `RuntimeError: erroneous average coalescence time`, no final model has been produced. This means the fitted model became numerically invalid, often because the automatic or highly flexible trajectory is poorly supported by the limited data. Confirm that the mask and scaffold length are correct, then try the bounded commands above. Failure even with the bounded model is itself an informative result: this one scaffold is insufficient for a stable SMC++ fit, and more independent scaffolds should be supplied rather than repeatedly tuning parameters until a curve appears.

Fitting may take much longer than the classroom slot. SMC++ folds the frequency-spectrum information by default; do not use `--unfold` unless ancestral alleles have been independently established.

## Step 5 — Plot and discuss

**Purpose:** convert the two final model files into one figure on a common time scale. `-g 11` changes generations into years using an 11-year brown-bear generation time, while `--logy` places Ne on a logarithmic axis.

```bash
smc++ plot -g 11 \
  "$OUTDIR/SMCPP_ABB_SBB_${CHROM}.pdf" \
  "$OUTDIR/ABB_bounded/model.final.json" \
  "$OUTDIR/SBB_bounded/model.final.json"
```

**Expected:** **results/smcpp/SMCPP_ABB_SBB_Scaffold_25.pdf** containing one trajectory for each final model. The x-axis is years before present and the y-axis is effective population size. Both represent model estimates, not census sizes. Features near the chosen time boundaries and sharp changes supported by only one scaffold deserve particular caution.

**Check:**

```bash
ls -lh "$OUTDIR/SMCPP_ABB_SBB_${CHROM}.pdf"
```

The PDF should have a nonzero size and show both populations. Inspect the iteration history and input limitations before interpreting it. Ask: What changes when the same broad method uses 10 or 8 individuals rather than one? Which differences might instead come from the distinguished bear, the shared mask, or the single scaffold?

SMC++ assumes sites absent from a variant-only VCF are homozygous reference unless masked, so mask construction is central to the interpretation. The classroom VCF was filtered to exclude sites with missing genotypes; if a genuinely variable, otherwise-callable site was removed by that filter, SMC++ would misread its absence as reference homozygosity. A full analysis needs input and callability definitions designed together.

Return to the [synthesis](04-synthesis.md).
