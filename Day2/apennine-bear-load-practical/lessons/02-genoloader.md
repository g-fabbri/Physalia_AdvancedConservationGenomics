# Part 2 — Which alleles are derived, and who carries them?

Estimated terminal time: 35 minutes.

The VCF `REF` allele is defined by the Apennine assembly. It need not be ancestral. GenoLoader uses genotypes in an outgroup to polarize each annotated biallelic SNP, recoding every individual as `0`, `1`, or `2` derived copies. We use black and polar bears together as a **strict-consensus outgroup**: disagreement or missing calls should be excluded from the conservative summary, not silently converted into ancestry.

## Step 1 — Build the combined outgroup list

**Purpose:** ask whether black and polar bears support the same ancestral allele.

**Input:** the two outgroup sample lists from Part 0.

~~~bash
cat "$BLB_LIST" "$POB_LIST" | sort -u > "$OUTDIR/outgroups.samples"
N_ABB=$(wc -l < "$ABB_LIST")
N_SBB=$(wc -l < "$SBB_LIST")
N_OUT=$(wc -l < "$OUTDIR/outgroups.samples")
printf 'ABB=%s SBB=%s outgroups=%s\n' "$N_ABB" "$N_SBB" "$N_OUT"
~~~

**Expected:** the outgroup count equals BLB plus POB counts if no IDs overlap. Outgroup support is a **hypothesis**, not a guarantee: shared ancestral polymorphism, sequencing error, and bear introgression can cause disagreement.

## Step 2 — Run GenoLoader

**Purpose:** make a derived-allele genotype table from the SnpEff-annotated VCF.

**Input:** `ANN`-annotated uncompressed VCF, ABB/SBB sample lists, and the combined outgroup list.

~~~bash
GENOLOADER=software/genoloader/genoloader
"$GENOLOADER" "$ANNOTATED" \
  --p1 "$ABB_LIST" --p2 "$SBB_LIST" \
  --p0 "$OUTDIR/outgroups.samples" \
  --m1 "$N_ABB" --m2 "$N_SBB" --m0 "$N_OUT" \
  --polX POP_OUT --low_cov NO
~~~

`--p1` and `--p2` define the focal groups; `--p0` supplies outgroup genotypes. The `--m*` values request called individuals in each group. `POP_OUT` uses outgroup information, but its output flags reveal whether a locus was truly outgroup-polarized or used a fallback. `--low_cov NO` keeps diploid dosage rather than one-read resampling.

**Expected:** a table named like `results/genetic_load/Bears_4pops_s25.ann.vcf.POP_OUT.gt`. The exact path is printed by GenoLoader; confirm it before the next command.

**Check:**

~~~bash
GT="$ANNOTATED.POP_OUT.gt"
head -n 3 "$GT"
awk -F'\t' 'NR>1 {n[$5]++} END {for (flag in n) print flag,n[flag]}' "$GT" | sort
~~~

The header should contain `scaffold`, `position`, `effect`, `vartype`, `flag`, `ref`, and sample IDs. Dosages are `0`, `1`, `2`, or missing. The flag counts show how many loci have reliable outgroup polarization versus fallback or ambiguous states.

## Step 3 — Summarize putative burden proxies

**Purpose:** compare a predicted damaging class with a synonymous comparator while respecting genotype state.

~~~bash
python scripts/summarize_genoloader.py \
  "$GT" "$ABB_LIST" "$SBB_LIST" \
  "$OUTDIR/derived_burden_by_sample.tsv"

column -t "$OUTDIR/derived_burden_by_sample.tsv" | head -n 16
~~~

**Expected:** three rows per ABB/SBB individual: `HIGH`, `missense`, and `synonymous`. The table reports called sites, derived copies (`heterozygote=1`, derived homozygote=2), homozygous-derived sites, and copies per called site. The script retains only flags indicating a monomorphic outgroup or an `allFix` site; all other fallback and ambiguous flags are excluded.

**Check:** the script prints how many annotated loci it read and retained after conservative flag filtering. If `called_sites` is zero for a category, do not compare its ratio. A transcript's impact label is not a direct estimate of deleteriousness. Every result here is a **Scaffold_25 burden proxy**, not a whole-genome load estimate.

## Step 4 — Interpret, then challenge, the comparison

1. Does ABB have more homozygous-derived `HIGH` sites per individual than SBB? Is the same true for total derived copies?
2. Are any differences also present at synonymous sites? What would that imply about ancestry, sampling, or technical bias?
3. Why can a small, inbred population expose recessive alleles while also losing some strongly deleterious variants through drift and purging?
4. How many loci were excluded because the outgroups disagreed or were not confidently callable?
5. What changes if black bears and polar bears are used **separately** as outgroups? Treat that as a sensitivity analysis, not a way to choose the preferred answer.

Write a cautious conclusion in the [answer sheet](../answers/student_answers.md), then continue to [GERP](03-gerp.md).
