# Part 3 — Where does evolutionary constraint add information?

Estimated terminal time: 15 minutes if a lifted score track is supplied. The liftOver itself is instructor preparation.

SnpEff predicts a variant's consequence relative to gene models. GERP measures evolutionary constraint from a multi-species alignment. A high positive rejected-substitution score suggests stronger constraint at that aligned position; it does not identify the derived allele or prove the mutation is damaging in bears.

## Before using the track

The instructor must identify the score track's actual coordinate reference. Polar bear appearing among 92 aligned mammals does not mean that the released score coordinates are polar-bear coordinates. If the track is on a polar-bear assembly, a polar-to-Apennine chain is needed. If it is on another reference, that assembly is the liftOver source instead.

Provide a one-base-per-row BED score track already mapped uniquely to the Apennine assembly: **data/GERP_on_Apennine_unique.bed.gz**. Its first four columns must be chromosome, zero-based start, end, and numeric GERP score. Document source assembly, chain direction, one-to-one mapping rate, rejected loci, and score sign convention.

## Step 1 — Verify the supplied coordinates

**Purpose:** confirm that score and VCF chromosome names match and that each score covers one base.

~~~bash
GERP=data/GERP_on_Apennine_unique.bed.gz
zcat "$GERP" | head -n 5
zcat "$GERP" | awk 'NR<=1000 && $3-$2!=1 {bad++} END {print "Non-single-base intervals in first 1000:",bad+0}'
~~~

**Expected:** BED rows such as **Scaffold_34  12344  12345  2.7**; the width check should print zero. BED coordinate 12344–12345 corresponds to VCF position **12345**.

## Step 2 — Compare scores across annotated sites

**Purpose:** attach a positional conservation score to each polarizable SNP. This is a descriptive comparison, not a load estimator.

~~~bash
awk -F'\t' 'BEGIN {OFS="\t"} NR>1 && $2~/^[0-9]+$/ {print $1,$2-1,$2,$3,$4,$5}' "$GT" \
  > "$OUTDIR/genoloader_sites.bed"

bedtools intersect \
  -a "$OUTDIR/genoloader_sites.bed" -b "$GERP" \
  -wa -wb > "$OUTDIR/sites_with_gerp.tsv"

wc -l "$OUTDIR/genoloader_sites.bed" "$OUTDIR/sites_with_gerp.tsv"
head "$OUTDIR/sites_with_gerp.tsv"
~~~

**Expected:** joined rows; the last column is the GERP score if the supplied BED has exactly four columns. The joined count should be no larger than the site count if scores are unique. A larger count signals duplicate or overlapping score intervals.

**Check:** what fraction of polarizable sites has a mapped score? Missing scores are not zero scores. Ask whether score availability differs among coding categories or genomic regions.

## Synthesis questions

1. Are HIGH or missense sites generally more constrained than synonymous sites in the mapped subset?
2. Can a synonymous site have a high GERP score? Can a missense site have a low one?
3. Why must score-track coverage and liftOver failures be reported before comparing populations?
4. Why does the same GERP score apply to an ABB and SBB allele at the same genomic position, while derived dosage can differ?
5. What would be needed to move from a putative burden proxy toward a defensible estimate of realized genetic load?

Record your final interpretation in the [answer sheet](../answers/student_answers.md).
