# Part 2 — Which alleles are derived, and who carries them?


## Start your terminal

From the Day 2 directory containing `data/`, `software/`, and `results/`, run this in each new terminal:

```bash
conda activate bear-load-practical
COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$COURSE_DIR/software/genoloader:$PATH"
```

```bash
OUTDIR="$COURSE_DIR/results/genoloader"
ANNOTATED="$COURSE_DIR/results/snpeff/Bears_4pops_s25.ann.vcf"
ABB_LIST="$COURSE_DIR/data/ABB.samples"
SBB_LIST="$COURSE_DIR/data/SBB.samples"
BLB_LIST="$COURSE_DIR/data/BLB.samples"
POB_LIST="$COURSE_DIR/data/POB.samples"
GENOLOADER="$COURSE_DIR/software/genoloader/genoloader"
mkdir -p "$OUTDIR"
```

Conda supplies the shared tools, and GenoLoader lives in Day 2 `software/genoloader/`. `ANNOTATED` is the result produced in Part 1, while `OUTDIR` stores the sample lists and summary tables made in this lesson. Defining all paths here makes the lesson safe to start in a new terminal.

GenoLoader has no separate output-directory option. It removes the final `.vcf` extension, adds the polarization mode and `.gt`, and writes the result beside the input. Therefore, the annotated VCF remains in `results/snpeff/`, and GenoLoader creates `results/snpeff/Bears_4pops_s25.ann.POP_OUT.gt`. Nothing is copied or moved.

**Check:**

```bash
ls -lh "$ANNOTATED" "$ABB_LIST" "$SBB_LIST" "$BLB_LIST" "$POB_LIST"
```

```bash
test -x "$GENOLOADER" && echo "GenoLoader ready: $GENOLOADER"
```

Every input must exist, and the last command should print the executable path.

The VCF `REF` allele is defined by the Apennine assembly. It need not be ancestral. GenoLoader uses genotypes in an outgroup to polarize each annotated biallelic SNP, recoding every individual as `0`, `1`, or `2` derived copies. We use black and polar bears together as a **strict-consensus outgroup**: disagreement or missing calls should be excluded from the conservative summary, not silently converted into ancestry.

The workflow is:

~~~text
SnpEff-annotated VCF + ABB/SBB/outgroup sample lists
                         ↓ GenoLoader
polarized per-individual derived-allele dosage table (.gt)
                         ↓ conservative flag filtering
HIGH, missense, and synonymous burden proxies by individual
~~~

## Step 1 — Build the combined outgroup list

**Purpose:** ask whether black and polar bears support the same ancestral allele.

**Input:** the two outgroup sample lists from Part 0.

### 1.1 — Combine and count outgroups

**Purpose:** create one nonredundant list and set the group sizes required by GenoLoader.

~~~bash
cat "$BLB_LIST" "$POB_LIST" | sort -u > "$OUTDIR/outgroups.samples"
~~~

~~~bash
N_ABB=$(wc -l < "$ABB_LIST")
N_SBB=$(wc -l < "$SBB_LIST")
N_OUT=$(wc -l < "$OUTDIR/outgroups.samples")
~~~

~~~bash
printf 'ABB=%s SBB=%s outgroups=%s\n' "$N_ABB" "$N_SBB" "$N_OUT"
~~~

`cat` joins the black-bear and polar-bear IDs, while `sort -u` removes any accidental duplicate. The three `N_*` variables are later used as missing-data thresholds: GenoLoader will require that many called individuals in each group.

**Expected:** the outgroup count equals the BLB count plus the POB count if no IDs overlap. The ABB and SBB values should match the sample counts recorded in Part 0.

**Check:**

~~~bash
cat "$OUTDIR/outgroups.samples"
~~~

~~~bash
comm -12 <(sort "$BLB_LIST") <(sort "$POB_LIST")
~~~

The first command displays the combined list. The second should print nothing; any output would identify an ID assigned to both outgroup lists.

Outgroup support is a **hypothesis**, not a guarantee. Shared ancestral polymorphism, sequencing or genotype error, and historical introgression among bear lineages can cause the two outgroup species to disagree.

## Step 2 — Run GenoLoader

**Purpose:** make a derived-allele genotype table from the SnpEff-annotated VCF.

**Input:** `ANN`-annotated uncompressed VCF, ABB/SBB sample lists, and the combined outgroup list.

### 2.1 — Polarize annotated SNPs

**Purpose:** run GenoLoader with ABB and SBB as focal groups and the combined outgroup list as ancestral-state evidence.

~~~bash
"$GENOLOADER" "$ANNOTATED" \
  --p1 "$ABB_LIST" --p2 "$SBB_LIST" \
  --p0 "$OUTDIR/outgroups.samples" \
  --m1 "$N_ABB" --m2 "$N_SBB" --m0 "$N_OUT" \
  --polX POP_OUT --low_cov NO
~~~

For each biallelic annotated SNP, GenoLoader reads the focal and outgroup genotypes, chooses an ancestral allele according to `POP_OUT`, and recodes all VCF samples relative to that allele. If the inferred ancestral allele is the VCF ALT allele, GenoLoader reverses the orientation so the original ALT copies become ancestral and original REF copies become derived.

| Argument | Meaning |
|---|---|
| **$ANNOTATED** | SnpEff-annotated, uncompressed input VCF; records without `ANN` are skipped |
| **--p1 / --p2** | ABB and SBB sample lists used as focal populations |
| **--p0** | Combined black- and polar-bear list used for outgroup polarization |
| **--m1 / --m2 / --m0** | Minimum numbers of called individuals required in the three groups; here they equal the full group sizes |
| **--polX POP_OUT** | Infer ancestry from an allele fixed in the outgroup when possible and record alternative cases with diagnostic flags |
| **--low_cov NO** | Use diploid GT dosages; do not perform one-read pseudo-haploid resampling |

The population lists guide polarization, but GenoLoader writes dosage columns for **all samples present in the VCF**, including outgroups. We later summarize only ABB and SBB.

**Expected terminal output:** GenoLoader reports the number of loci written, the number re-polarized relative to VCF REF, and the number failing the requested missingness thresholds. A zero or unexpectedly small retained count should trigger checks of sample IDs, group sizes, `ANN`, and genotype completeness.

Define the path of the table created by GenoLoader:

~~~bash
GT="${ANNOTATED%.vcf}.POP_OUT.gt"
~~~

`${ANNOTATED%.vcf}` means “use the value of `ANNOTATED` after removing the final `.vcf`”. Therefore, `GT` becomes **results/snpeff/Bears_4pops_s25.ann.POP_OUT.gt**. It is a tab-delimited text table, not a VCF, and no original genotypes are modified.

**Check:**

~~~bash
ls -lh "$GT"
~~~

~~~bash
wc -l "$GT"
~~~

The file must have a nonzero size and more than one line. One line is the header; the remaining lines are annotated loci written by GenoLoader.

### 2.2 — Inspect dosages and polarization flags

**Purpose:** verify the table structure and count reliable versus fallback or ambiguous polarization states.

~~~bash
head -n 3 "$GT"
~~~

~~~bash
awk -F'\t' 'NR>1 {n[$5]++} END {for (flag in n) print flag,n[flag]}' "$GT" | sort
~~~

The output columns mean:

| Column | Meaning |
|---|---|
| **scaffold / position** | Genomic location copied from VCF CHROM and POS |
| **effect** | SnpEff impact category, such as HIGH, MODERATE, LOW, or MODIFIER |
| **vartype** | Simplified annotation class, such as missense, synonymous, intron, or intergenic |
| **flag** | Trace of how the ancestral allele was chosen and whether fallback information was needed |
| **ref** | GenoLoader's inferred ancestral allele; it is not necessarily the original VCF REF allele |
| **sample columns** | Derived-allele dosage: 0 for ancestral homozygote, 1 for heterozygote, 2 for derived homozygote, and `nan` for missing |

The flag counts show how many loci were supported by a monomorphic outgroup and how many used ambiguous or fallback polarization. Flags beginning with `unfolded` indicate a monomorphic outgroup, but the suffix still records the ingroup configuration—for example fixation, segregation, missingness, or possible incomplete lineage sorting. `allFold`, `inFold`, and other folded/fallback flags should not be silently treated as secure outgroup polarization.

**Check:** the header should include every expected focal individual. Confirm explicitly:

~~~bash
for SAMPLE in $(cat "$ABB_LIST" "$SBB_LIST"); do
  head -n 1 "$GT" | tr '\t' '\n' | grep -Fx "$SAMPLE" >/dev/null || \
    echo "Missing sample column: $SAMPLE"
done
~~~

No output means that every ABB and SBB sample was found.

## Step 3 — Summarize putative burden proxies

**Purpose:** compare a predicted damaging class with a synonymous comparator while respecting genotype state.

### 3.1 — Count derived alleles by individual and class

**Purpose:** summarize called sites, derived copies, and homozygous-derived sites for ABB and SBB.

~~~bash
python scripts/summarize_genoloader.py \
  "$GT" "$ABB_LIST" "$SBB_LIST" \
  "$OUTDIR/derived_burden_by_sample.tsv" \
  "$OUTDIR/derived_frequency_by_impact.tsv" \
  "$OUTDIR/GenoLoader_ABB_SBB_Scaffold_25"
~~~

~~~bash
column -t "$OUTDIR/derived_burden_by_sample.tsv" | head -n 16
~~~

~~~bash
column -t "$OUTDIR/derived_frequency_by_impact.tsv"
~~~

The helper script does not re-run polarization. It reads the `.gt` table, keeps conservative outgroup-supported flags, assigns each retained locus to a comparison category, and sums each focal individual's dosages. It also calculates population mean derived-allele frequencies and `Rxy` for the four SnpEff impact classes.

| Output column | Meaning |
|---|---|
| **population / sample** | Focal group and individual ID |
| **category** | HIGH impact, missense, or synonymous comparator |
| **called_sites** | Retained category sites with a non-missing dosage for that individual |
| **derived_copies** | Sum of dosage values: heterozygote contributes 1 and derived homozygote contributes 2 |
| **homozygous_derived** | Number of sites with dosage 2; a proxy relevant to recessive effects |
| **derived_copies_per_called_site** | Derived copies divided by called sites, helping account for different denominators |
| **homozygous_derived_per_called_site** | Homozygous-derived sites divided by called sites |

**Expected:** three rows per ABB/SBB individual: `HIGH`, `missense`, and `synonymous`. The script retains flags beginning with `unfolded`—a monomorphic combined outgroup—or `allFix`; other fallback and polymorphic-outgroup flags are excluded.

The second TSV has one row for each SnpEff impact class:

| Output column | Meaning |
|---|---|
| **impact** | SnpEff `MODIFIER`, `LOW`, `MODERATE`, or `HIGH` class |
| **shared_called_sites** | Retained sites with at least one called individual in both ABB and SBB |
| **mean_DAF_ABB / mean_DAF_SBB** | Mean derived-allele frequency across those sites in each population |
| **Rxy_ABB_over_SBB** | Sum of ABB-derived/SBB-ancestral contributions divided by the reverse contributions |

`Rxy > 1` means the derived allele is relatively more frequent in ABB than SBB for that class; `Rxy < 1` means the reverse. This is directional and depends on which population is placed in the numerator.

The final argument is an output **prefix**, not a complete filename. The script adds an informative suffix and `.pdf` for each separate figure.

**Check:** the script prints how many annotated loci it read and retained after conservative flag filtering, followed by the paths of both TSVs and all PDFs. If `called_sites` is zero for a category, do not compare its ratio. A transcript's impact label is not a direct estimate of deleteriousness. Every result here is a **Scaffold_25 burden proxy**, not a whole-genome load estimate.

This summary deliberately reports several proxies instead of one number called “genetic load.” Derived copies are closer to an additive count, while homozygous-derived sites are informative for completely recessive models. Neither incorporates selection coefficients, dominance values, expression, or validated phenotypic effects.

### 3.2 — Inspect the population comparison plot

**Purpose:** visualize individual variation and population means without confusing different numbers of annotated sites with different burdens.

~~~bash
ls -lh "$OUTDIR"/GenoLoader_ABB_SBB_Scaffold_25_*.pdf
~~~

The script creates seven separate figures:

| Figure | Contents |
|---|---|
| **..._normalized_burden.pdf** | Derived copies per called site and homozygous-derived sites per called site for HIGH, missense, and synonymous variants |
| **..._MODIFIER.pdf** | Total derived copies, heterozygous sites, and homozygous-derived sites for MODIFIER variants |
| **..._LOW.pdf** | The same three measures for LOW-impact variants |
| **..._MODERATE.pdf** | The same three measures for MODERATE-impact variants |
| **..._HIGH.pdf** | The same three measures for HIGH-impact variants |
| **..._total_derived_by_impact.pdf** | Four side-by-side panels comparing total derived copies for HIGH, MODERATE, LOW, and MODIFIER variants |
| **..._frequency_Rxy.pdf** | Mean derived-allele frequency for ABB and SBB plus directional `Rxy` for every impact class |

In every individual-level panel, each point is one bear and the short horizontal line is the population mean. **Total derived copies** is calculated as `heterozygous sites + 2 × homozygous-derived sites`. Therefore, a difference in total derived copies can be decomposed into its heterozygous and homozygous-derived contributions in the adjacent panels.

The combined **total-derived-by-impact** figure follows the faceted style shown in class. ABB and SBB are on the x-axis of every panel, individual bears are the colored points, and black horizontal bars mark population means. The panels have separate y-axis scales because the four impact classes contain very different numbers of variants; compare ABB with SBB **within** a panel rather than comparing the absolute heights among panels.

Start with the normalized-burden figure because it accounts for each individual's number of called sites. The four impact-specific figures explain whether the pattern arises through more heterozygous sites, more homozygous-derived sites, or both, but their raw counts can change with the number of retained and callable sites. The frequency/`Rxy` figure summarizes allele-frequency shifts; it does not describe dominance or individual genotype state.

Do not infer statistical significance from overlapping or separated points in these teaching plots. There are only 10 ABB and 8 SBB individuals, all loci come from one scaffold, individuals may be related, and sites are not statistically independent. Neither SnpEff impact, derived-allele frequency, nor `Rxy` is a direct estimate of selection coefficients, fitness, or realized genetic load.

## Step 4 — Interpret, then challenge, the comparison

### 4.1 — Compare the populations cautiously

**Purpose:** separate observed Scaffold_25 patterns from claims about genome-wide genetic load.







