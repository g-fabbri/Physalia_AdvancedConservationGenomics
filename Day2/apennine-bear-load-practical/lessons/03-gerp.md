# Part 3 — Add GERP constraint scores to Scaffold_25

Estimated terminal time: 15 minutes **if the instructor supplies the chain and score files**. Building the assembly alignment is instructor preparation; the full process is shown so the liftOver is transparent and can be repeated outside class.

SnpEff predicts consequences from gene models; GERP measures evolutionary constraint at an alignment column. Neither identifies a bear allele's fitness effect. We will score the **Scaffold_25 SNPs** used in Parts 1–2.

The GERP bigWig's **coordinate assembly must be verified first**. Polar bear being among the 92 aligned mammals does *not* prove the bigWig uses polar coordinates. The commands below assume the score track uses the exact polar assembly in `POLAR_FA`. If it does not, stop and identify the correct assembly. We lift Apennine SNPs to the score assembly and carry the score back using a stable site ID; there is no need to lift a whole bigWig or make separate files for each SnpEff category.

## Step 1 — Build the one-scaffold chain (instructor preparation)

**Purpose:** map Apennine Scaffold_25 to the polar coordinate system of the verified score track.

**Input:** the exact Apennine VCF-reference FASTA, the polar score-reference FASTA, minimap2, and Transanno. [Part 0](00-inputs.md) checks the matching `data/Bears_4pops_s25.vcf.gz`. If that file is absent, extract Scaffold_25 from the verified, indexed **genome-wide four-population VCF** before class; do not extract from a VCF containing a different scaffold.

Only when the Scaffold_25 VCF has **not** already been prepared, and after confirming the source contains the intended four populations, create it once:

```bash
SOURCE_VCF=/jarvis/scratch/usr/biello/bear/snpeff/VCF/marpolblk.sorted.merged.alignable.final.SNP.vcf.gz
bcftools view -r Scaffold_25 -Oz \
  -o data/Bears_4pops_s25.vcf.gz "$SOURCE_VCF"
bcftools index -f data/Bears_4pops_s25.vcf.gz
bcftools query -l data/Bears_4pops_s25.vcf.gz
```

**Expected:** only Scaffold_25 records and the expected four-population sample IDs. This is extraction, not a new allele-frequency filter; Part 0 checks site representation. Do not overwrite an existing verified teaching VCF.

```bash
CHROM=Scaffold_25
APP_FA=/jarvis/data/refgenomes/Uarcmar/mUrsArc1.1.primarysoftmask.fasta
POLAR_FA=/path/to/verified/polar_score_assembly.fasta
TRANSANNO=/jarvis/scratch/usr/benazzo/endemixit/ursusarctos/liftover/transanno-0.2.4/transanno
PREP=results/gerp/preparation
mkdir -p "$PREP"

samtools faidx "$APP_FA" "$CHROM" > "$PREP/${CHROM}.fa"
minimap2 -cx asm20 --cs -t 8 "$POLAR_FA" "$PREP/${CHROM}.fa" \
  > "$PREP/Apennine_to_polar.paf"
"$TRANSANNO" minimap2chain "$PREP/Apennine_to_polar.paf" \
  --output "$PREP/Apennine_to_polar.chain"
head -n 1 "$PREP/Apennine_to_polar.chain"
```

**Expected:** non-empty one-scaffold FASTA, PAF, and chain files. Transanno's documented convention puts the destination/query FASTA first in minimap2 and the source/reference FASTA second. Its current documentation calls the conversion command `minimap2-to-chain`; check `"$TRANSANNO" --help` if the installed 0.2.4 `minimap2chain` spelling differs. `asm20` is an alignment preset, **not** a hard 5% divergence cutoff. [Transanno](https://github.com/informationsea/transanno), [minimap2](https://github.com/lh3/minimap2/blob/master/cookbook.md).

**Check:** compare PAF names and chain header with the two FASTAs. Test a few known positions: the chain must accept **Apennine** coordinates and produce **polar** coordinates. Do not swap file labels simply to make the command finish.

## Step 2 — Make a named BED file of VCF SNPs

**Purpose:** assign each Apennine SNP an ID that survives liftOver.

**Input:** the same indexed Scaffold_25 VCF used for SnpEff and GenoLoader.

```bash
VCF=data/Bears_4pops_s25.vcf.gz
bcftools query -f '%CHROM\t%POS\n' "$VCF" |
  awk -v c="$CHROM" 'BEGIN{OFS="\t"} $1==c {print $1,$2-1,$2,$1":"$2}' |
  sort -k1,1 -k2,2n -u > "$PREP/Apennine_sites.bed"
head "$PREP/Apennine_sites.bed"
wc -l "$PREP/Apennine_sites.bed"
```

**Expected:** BED4 rows such as `Scaffold_25  12344  12345  Scaffold_25:12345`. BED starts are zero-based; the VCF position is one-based.

**Check:** compare the site count with the VCF record count from Part 0. A large difference suggests duplicated positions or a chromosome-name mismatch.

## Step 3 — Lift sites to polar and discard ambiguous mappings

**Purpose:** use the chain to locate the SNPs on the GERP score assembly. `-multiple` exposes alternative mappings; keep only IDs with exactly one one-base result.

**Input:** the BED from Step 2 and validated chain from Step 1.

```bash
liftOver -multiple "$PREP/Apennine_sites.bed" \
  "$PREP/Apennine_to_polar.chain" \
  "$PREP/polar_sites_all.bed" "$PREP/unmapped.bed"

awk 'BEGIN{OFS="\t"} {n[$4]++; line[$4]=$0}
     END{for(id in n) if(n[id]==1) print line[id]}' \
  "$PREP/polar_sites_all.bed" |
  awk 'BEGIN{OFS="\t"} $3-$2==1 {print}' |
  sort -k1,1 -k2,2n > "$PREP/polar_sites_unique.bed"

wc -l "$PREP/Apennine_sites.bed" "$PREP/polar_sites_all.bed" \
  "$PREP/polar_sites_unique.bed"
head "$PREP/polar_sites_unique.bed"
```

**Expected:** columns 1–3 now contain **polar** coordinates; column 4 still contains the original Apennine `Scaffold_25:position` ID. The unique count is no greater than the input count.

**Check:** inspect `unmapped.bed`, count ambiguous and length-changing loci, and report the uniquely mapped fraction. Unmapped sites have *missing* scores, not score zero.

## Step 4 — Read GERP scores at the polar positions

**Purpose:** extract bigWig scores only from polar regions reached by the sites, then overlap them with the one-base polar BED. bedGraph intervals already carry scores; `bedtools makewindows -w 1` is unnecessary.

**Input:** verified polar-coordinate GERP bigWig and unique polar sites. The instructor may provide the resulting bedGraph checkpoint if this extraction is too slow in class.

```bash
GERP_BW=/path/to/verified/polar_coordinate_GERP.bigWig
awk 'BEGIN{OFS="\t"}
     {if(!($1 in min) || $2<min[$1]) min[$1]=$2;
      if(!($1 in max) || $3>max[$1]) max[$1]=$3}
     END{for(c in min) print c,min[c],max[c]}' \
  "$PREP/polar_sites_unique.bed" |
while IFS=$'\t' read -r chr start end; do
  bigWigToBedGraph -chrom="$chr" -start="$start" -end="$end" \
    "$GERP_BW" "$PREP/GERP_${chr}.bedGraph"
done

cat "$PREP"/GERP_*.bedGraph |
  sort -k1,1 -k2,2n > "$PREP/GERP_polar_regions.bedGraph"
bedtools intersect \
  -a "$PREP/polar_sites_unique.bed" \
  -b "$PREP/GERP_polar_regions.bedGraph" \
  -wa -wb > "$PREP/sites_with_polar_scores.tsv"
head "$PREP/sites_with_polar_scores.tsv"
```

**Expected:** eight columns: the four mapped-site columns followed by four polar bedGraph columns. Column 4 is the original Apennine site ID and column 8 the GERP score.

**Check:** `bigWigToBedGraph` reporting an unknown chromosome indicates a likely assembly or contig-name mismatch. Missing values are not zero. Check for sites matching more than one score interval.

## Step 5 — Restore Apennine coordinates and validate

**Purpose:** write a four-column, one-base score table in the coordinates used by SnpEff and GenoLoader. Retain only IDs with one numeric score match.

**Input:** the polar-site/score join from Step 4.

```bash
awk 'BEGIN{OFS="\t"}
     $8 ~ /^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/ {
       n[$4]++; score[$4]=$8
     }
     END {for(id in n) if(n[id]==1) {
       split(id,a,":"); print a[1],a[2]-1,a[2],score[id]
     }}' "$PREP/sites_with_polar_scores.tsv" |
  sort -k1,1 -k2,2n > "$PREP/GERP_on_Apennine_unique.bed"

GERP=results/gerp/GERP_on_Apennine_unique.bed.gz
bgzip -c "$PREP/GERP_on_Apennine_unique.bed" > "$GERP"
tabix -f -p bed "$GERP"
zcat "$GERP" | head
zcat "$GERP" |
  awk '$1!="Scaffold_25" || $3-$2!=1 {bad++} END{print "Invalid rows:",bad+0}'
wc -l "$PREP/Apennine_sites.bed" "$PREP/GERP_on_Apennine_unique.bed"
```

**Expected:** `Scaffold_25`, zero-based start, end, numeric score. `Invalid rows: 0`; fewer scored sites than input sites is normal.

**Check:** spot-check positions against both FASTAs and the bigWig for one-base shifts and contig mismatches. Report the fraction of VCF sites with a unique score. This is a **site-only** track, not a whole-scaffold GERP profile. liftOver transfers coordinates; it does not recalculate GERP or infer ancestral bear alleles.

## Step 6 — Compare constraint with annotated sites

**Purpose:** attach scores to the polarizable SNP table from Part 2. This is a descriptive comparison, not a direct estimate of realized genetic load.

**Input:** `GT` from [GenoLoader](02-genoloader.md) and `GERP` from Step 5. With a precomputed track, point `GERP` to its actual path and start here.

```bash
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $2~/^[0-9]+$/ {print $1,$2-1,$2,$3,$4,$5}' "$GT" \
  > "$OUTDIR/genoloader_sites.bed"
bedtools intersect \
  -a "$OUTDIR/genoloader_sites.bed" -b "$GERP" \
  -wa -wb > "$OUTDIR/sites_with_gerp.tsv"
wc -l "$OUTDIR/genoloader_sites.bed" "$OUTDIR/sites_with_gerp.tsv"
head "$OUTDIR/sites_with_gerp.tsv"
```

**Expected:** joined rows with the GERP score in the last column. If scores are unique, the joined count cannot exceed the input site count.

**Check:** a larger joined count signals duplicate/overlapping score intervals. Ask whether score availability differs among HIGH, missense, and synonymous categories. The same positional score applies to ABB and SBB at one site; their derived dosages may differ.

## Stop and discuss

1. Are HIGH or missense sites more constrained than synonymous sites in the uniquely mapped subset?
2. Why are unscored or ambiguously mapped positions not assigned score zero?
3. How could chain direction, assembly mismatch, or a one-base BED error produce convincing but wrong results?
4. Why are SnpEff impact, GERP constraint, derived status, and realized fitness load distinct?

Record your interpretation in the [answer sheet](../answers/student_answers.md).
