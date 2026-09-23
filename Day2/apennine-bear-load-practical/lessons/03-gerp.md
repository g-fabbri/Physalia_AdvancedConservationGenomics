# Optional Part 3 — Add GERP constraint scores to Scaffold_25

Estimated terminal time: 15 minutes **if the instructor supplies the chain and score files**. Building the assembly alignment is instructor preparation; the full process is shown so the liftOver is transparent and can be repeated outside class.

This lesson is an **optional extension**. The core Day 2 practical is complete after GenoLoader; the instructor may teach GERP if time allows or leave it for independent work. No GERP output is needed for the core conclusion.

## Start your terminal

From the Day 2 directory containing `data/`, `software/`, and `results/`, run this in each new terminal:

```bash
conda activate bear-load-practical
COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$PATH"
```

Conda supplies the shared tools, including `minimap2`, `transanno`, `liftOver`, `bigWigInfo`, and `bigWigToBedGraph`. During setup, `bash software/link_conda_tools.sh` creates course-local links in Day 2 `software/bin/` where needed. Check the active commands before continuing:

```bash
command -v minimap2 transanno liftOver bigWigInfo bigWigToBedGraph
```

`transanno` is installed from Bioconda by `environment.yml`; students do **not** need to download it separately. If `command -v transanno` prints nothing, update the Conda environment rather than adding an unrelated binary manually. The [Bioconda Transanno recipe](https://bioconda.github.io/recipes/transanno/README.html) documents `conda install transanno`; the [Transanno repository](https://github.com/informationsea/transanno) provides releases and source-build instructions as alternatives when Conda is unavailable.

SnpEff predicts consequences from gene models; GERP measures evolutionary constraint at an alignment column. Neither identifies a bear allele's fitness effect. We will score the **Scaffold_25 SNPs** used in Parts 1–2.

### What is a bigWig?

A **bigWig** (`.bw`) is an indexed, compressed binary track of numerical values along a genome—for example, one GERP conservation score at a genomic position. Unlike a text BED or bedGraph file, it is not meant to be read with `head`. Its index lets `bigWigToBedGraph` retrieve a selected chromosome interval without converting the entire track. The extracted **bedGraph** is a small, readable table with `chrom start end score` columns. [UCSC bigWig guide](https://genome.ucsc.edu/goldenPath/help/bigWig).

The score track for this exercise is Ensembl release 114's [91-mammal GERP bigWig for polar bear](https://ftp.ensembl.org/pub/release-114/compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.ursus_maritimus.UrsMar_1.0.bw), named `gerp_conservation_scores.ursus_maritimus.UrsMar_1.0.bw`. Its filename identifies the intended **UrsMar_1.0 polar-bear coordinates**.

Use the polar-bear FASTA from the **same Ensembl release** as the conservation track: [Ensembl release 114 polar-bear DNA files](https://ftp.ensembl.org/pub/release-114/fasta/ursus_maritimus/dna/) and `Ursus_maritimus.UrsMar_1.0.dna.toplevel.fa.gz`. This toplevel FASTA uses the same GenBank-style `AVOR...` sequence identifiers as the bigWig. The `GCF_...` RefSeq version represents the same UrsMar_1.0 assembly sequence but uses `NW_...` headers; `bigWigToBedGraph` would therefore return empty output unless those accessions were translated.

The polar FASTA and bigWig must describe the same assembly and contig names. **Scaffold_25 is an Apennine scaffold, not a name to search for in the polar bigWig.** We first map that Apennine scaffold to its corresponding polar region(s), lift the SNP positions, and retrieve polar-coordinate GERP values. Stable site IDs then return those values to Apennine coordinates.

## Step 1 — Build the one-scaffold chain (instructor preparation)

**Purpose:** map Apennine Scaffold_25 to the polar coordinate system of the verified score track.

**Input:** the prepared Apennine Scaffold_25 VCF-reference FASTA, the complete polar score-reference FASTA, minimap2, and Transanno. [Part 0](00-inputs.md) has already checked `data/Bears_4pops_s25.vcf.gz`; this lesson does not recreate it.

### 1.1 — Name the assemblies and output directories

**Purpose:** keep the source and destination assemblies explicit throughout the liftOver. Replace `POLAR_FA` only after verifying the bigWig's assembly.

```bash
CHROM=Scaffold_25
APP_FA="$COURSE_DIR/data/mUrsArc1.1.genome.s25.fasta"
POLAR_FA="$COURSE_DIR/data/Ursus_maritimus.UrsMar_1.0.dna.toplevel.fa"
GERP_BW="$COURSE_DIR/data/gerp_conservation_scores.ursus_maritimus.UrsMar_1.0.bw"
VCF="$COURSE_DIR/data/Bears_4pops_s25.vcf.gz"
GT="$COURSE_DIR/results/snpeff/Bears_4pops_s25.ann.POP_OUT.gt"
TRANSANNO=$(command -v transanno)
OUTDIR="$COURSE_DIR/results/gerp"
PREP="$OUTDIR/preparation"
mkdir -p "$PREP"
```

`APP_FA` is already restricted to Scaffold_25, so no additional FASTA extraction is necessary. `POLAR_FA` is the Ensembl release 114 toplevel UrsMar_1.0 assembly matching the GERP coordinate system and names. `GERP_BW` is the large indexed score track; students need it only when reproducing instructor-preparation Steps 4–5.

If the polar FASTA has not already been prepared, the instructor can download and decompress the matching NCBI file once:

```bash
wget -c \
  https://ftp.ensembl.org/pub/release-114/fasta/ursus_maritimus/dna/Ursus_maritimus.UrsMar_1.0.dna.toplevel.fa.gz \
  -O "$POLAR_FA.gz"
gunzip -k "$POLAR_FA.gz"
```

Do not repeat the download for every student. Confirm the inputs and software:

```bash
ls -lh "$APP_FA" "$POLAR_FA" "$VCF" "$GT" "$GERP_BW"
test -n "$TRANSANNO" && "$TRANSANNO" minimap2chain --help | head
grep '^>' "$APP_FA" | head
grep '^>' "$POLAR_FA" | head
```

**Expected:** the Apennine FASTA contains `>Scaffold_25`; the Ensembl polar FASTA contains `>AVOR...` accessions; all files are nonempty; and Transanno prints help text. If the bigWig is not present, the class can still discuss Steps 1–3 and then use the instructor's precomputed Apennine-coordinate score track from Step 5.

### 1.2 — Verify the bigWig assembly and sequence names

**Purpose:** confirm that the polar FASTA and GERP bigWig use identical sequence identifiers before performing an expensive alignment. The assembly name alone is insufficient because RefSeq and GenBank releases of the same assembly can use different accessions.

```bash
bigWigInfo -chroms "$GERP_BW" | head -n 30
```

**Expected:** the summary reports thousands of sequences and the chromosome list begins with names such as `AVOR01003489.1`. These are GenBank accessions. `bigWigInfo` also confirms that the 7 GB file is a readable indexed bigWig rather than a truncated download.

Create sorted name lists and verify that the two resources overlap:

```bash
bigWigInfo -chroms "$GERP_BW" |
  awk '$2~/^[0-9]+$/ && $3~/^[0-9]+$/ {print $1}' |
  sort -u > "$PREP/GERP_contigs.txt"

grep '^>' "$POLAR_FA" |
  sed 's/^>//; s/[[:space:]].*$//' |
  sort -u > "$PREP/polar_FASTA_contigs.txt"

printf 'Shared FASTA/bigWig contig names: '
comm -12 "$PREP/polar_FASTA_contigs.txt" "$PREP/GERP_contigs.txt" |
  wc -l

comm -12 "$PREP/polar_FASTA_contigs.txt" "$PREP/GERP_contigs.txt" |
  head
```

**Expected:** a positive shared count and names beginning with `AVOR`. A count of `0` means the headers are incompatible. Do not continue with alignment: use the Ensembl release 114 toplevel FASTA above or explicitly translate accessions with an assembly report.

### 1.3 — Align the two assemblies

**Purpose:** find corresponding segments between Apennine Scaffold_25 and the verified polar assembly. In minimap2, the **first FASTA is the target/destination** and the second is the query/source. The output is a PAF alignment, **not** yet a liftOver chain.

```bash
minimap2 -cx asm20 --cs -t 2 "$POLAR_FA" "$APP_FA" \
  > "$PREP/Apennine_to_polar.paf"
```

```bash
head -n 2 "$PREP/Apennine_to_polar.paf"
```

**Expected:** PAF rows whose first field is `Scaffold_25` and whose sixth field names a polar contig. `-c` requests base-level alignment and a CIGAR-like `cg` tag; `--cs` writes detailed substitutions and gaps; `-t 2` uses eight threads.

`asm20` is a bundle of assembly-alignment parameters, not an instruction that the genomes differ by exactly 20% and not a hard divergence filter. Among the assembly presets, `asm5` is the strictest for highly similar assemblies, `asm10` is intermediate, and `asm20` is the most permissive. We use `asm20` to maintain sensitivity in a cross-species brown-bear-to-polar-bear alignment. Reasonable alternatives are:

- try `asm10` when the assemblies are expected to be highly similar and compare aligned coverage, mapping uniqueness, and chain quality;
- try `asm5` mainly for very closely related or within-species assemblies;
- retain `asm20` when stricter presets fragment or lose valid orthologous alignment.

Do not select a preset simply because it maps the most bases: permissive settings can also increase paralogous or repetitive mappings. Compare the uniquely lifted fraction and spot-check loci. The [minimap2 documentation](https://github.com/lh3/minimap2) describes `asm5` for intra-species assembly alignment and recommends tuning assembly presets to cross-species divergence.

**Check:** an empty PAF means there is no usable alignment.

### 1.4 — Convert the alignment to a chain

**Purpose:** convert the pairwise assembly alignment into the coordinate-mapping format read by `liftOver`.

```bash
"$TRANSANNO" minimap2chain "$PREP/Apennine_to_polar.paf" \
  --output "$PREP/Apennine_to_polar.chain"
```

```bash
head -n 8 "$PREP/Apennine_to_polar.chain"
```

A **UCSC chain file** describes how continuous blocks in one assembly correspond to blocks in another assembly. It does not contain DNA sequences or GERP scores. Instead, it records chromosome names, coordinate ranges, strand orientation, aligned-block sizes, and the gaps between successive blocks. `liftOver` follows this map to translate an interval from one coordinate system to the other.

Each alignment starts with a header shaped like this:

```text
chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
```

The fields mean:

- `score`: an alignment score used to rank chains; it is not a GERP score;
- `tName`, `tSize`, `tStart`, `tEnd`: target sequence name, total length, and covered range;
- `qName`, `qSize`, `qStart`, `qEnd`: query sequence name, total length, and covered range;
- `tStrand` and `qStrand`: alignment orientation; a minus query strand indicates a reverse-complement mapping;
- `id`: an identifier for that chain.

Coordinates in the chain format are zero-based and half-open. After the header come one or more block lines:

```text
size  dt  dq
size  dt  dq
size
```

`size` is the length of an aligned block; `dt` is the gap before the next block in the target; and `dq` is the corresponding gap in the query. The final block contains only its size. A chain file can contain many chains because one scaffold may align in several segments or to several destination contigs.



**Expected:** a non-empty file containing one or more headers beginning with `chain`, followed by numeric alignment-block lines. The header should contain an Apennine scaffold name and a polar-bear contig accession from the two FASTAs.

Count the chains and inspect their sequence names:

```bash
grep -c '^chain[[:space:]]' "$PREP/Apennine_to_polar.chain"
grep '^chain[[:space:]]' "$PREP/Apennine_to_polar.chain" | head
```

`[[:space:]]` deliberately accepts either a tab or a space after `chain`. Transanno commonly writes tab-delimited chain headers, so `grep '^chain '` can incorrectly report zero even when the chain file is valid.

**Check:** the existence of a chain does not prove that its direction is correct. Before processing all SNPs, test a few Apennine BED intervals. The chain used in Step 3 must accept **Apennine** coordinates and emit **polar** coordinates. If every test interval is unmapped, inspect the header and the minimap2/Transanno source–destination convention rather than reversing labels blindly. Also confirm that sequence names and sizes agree with the corresponding FASTA indexes. Confirm the installed syntax with `"$TRANSANNO" minimap2chain --help`. [Transanno documentation](https://github.com/informationsea/transanno).

## Step 2 — Make a named BED file of VCF SNPs

**Purpose:** assign each Apennine SNP an ID that survives liftOver.

**Input:** the same indexed Scaffold_25 VCF used for SnpEff and GenoLoader.

### 2.1 — Convert VCF positions to BED

**Purpose:** represent each one-based VCF position as a zero-based, one-base BED interval while retaining the original site ID.

```bash
bcftools query -f '%CHROM\t%POS\n' "$VCF" |
  awk -v c="$CHROM" 'BEGIN{OFS="\t"} $1==c {print $1,$2-1,$2,$1":"$2}' |
  sort -k1,1 -k2,2n -u > "$PREP/Apennine_sites.bed"
```

`bcftools query` reads chromosome and one-based position from the VCF; `awk` converts each SNP to a one-base BED interval and stores the original position in column 4. `sort -u` removes any exact duplicate site rows.

### 2.2 — Inspect the BED input

**Purpose:** catch coordinate or record-count mistakes before liftOver.

```bash
head "$PREP/Apennine_sites.bed"
wc -l "$PREP/Apennine_sites.bed"
```

**Expected:** BED4 rows such as `Scaffold_25  12344  12345  Scaffold_25:12345`. BED starts are zero-based; the VCF position is one-based.

**Check:** compare the site count with the VCF record count from Part 0. A large difference suggests duplicated positions or a chromosome-name mismatch.

## Step 3 — Lift sites to polar and discard ambiguous mappings

**Purpose:** use the chain to locate the SNPs on the GERP score assembly. `-multiple` exposes alternative mappings; keep only IDs with exactly one one-base result.

**Input:** the BED from Step 2 and validated chain from Step 1.

### 3.1 — Lift the candidate sites

**Purpose:** produce all candidate polar mappings and a separate list of unmapped sites.

```bash
liftOver -multiple "$PREP/Apennine_sites.bed" \
  "$PREP/Apennine_to_polar.chain" \
  "$PREP/polar_sites_all.bed" "$PREP/unmapped.bed"
```

The two outputs have different meanings: `polar_sites_all.bed` contains mapped coordinates, while `unmapped.bed` records failures and reasons. A single input site may have multiple mapped rows, which is why we have **not** called this file unique.

Multiple mappings can occur when an Apennine interval aligns to more than one polar region. Common causes include repetitive sequence, recent segmental duplication, paralogous regions, retained haplotigs or alternate assembly sequence, and overlapping primary/secondary chains produced by the assembly alignment. The permissive `asm20` preset can also retain more weak or repeated matches. These alternatives are useful to expose during quality control, but they do not tell us which polar position carries the orthologous GERP score. We therefore exclude them from this conservative exercise.

Check that both output files were created and inspect their structure:

```bash
ls -lh "$PREP/polar_sites_all.bed" "$PREP/unmapped.bed"
head "$PREP/polar_sites_all.bed"
grep -v '^#' "$PREP/unmapped.bed" | head
```

The mapped file should contain polar contig names in columns 1–3 and the original Apennine site ID in column 4. The unmapped file may contain comment lines beginning with `#`, followed by the original BED record.

Summarize the result before filtering:

```bash
printf 'Input sites: '
wc -l < "$PREP/Apennine_sites.bed"

printf 'Mapped output rows: '
wc -l < "$PREP/polar_sites_all.bed"

printf 'Distinct mapped input IDs: '
cut -f4 "$PREP/polar_sites_all.bed" | sort -u | wc -l

printf 'Distinct unmapped input IDs: '
grep -v '^#' "$PREP/unmapped.bed" | cut -f4 | sort -u | wc -l

printf 'Input IDs with multiple mapped rows: '
awk '{n[$4]++} END{for(id in n) if(n[id]>1) multiple++;
     print multiple+0}' "$PREP/polar_sites_all.bed"

printf 'Mapped rows not one base wide: '
awk '$3-$2!=1 {n++} END{print n+0}' "$PREP/polar_sites_all.bed"
```

**Expected:** mapped plus unmapped distinct input IDs should approximately account for the original input sites. The number of mapped rows can exceed the number of distinct mapped IDs when `-multiple` finds alternatives. If nearly everything is unmapped, first suspect chain direction or sequence-name incompatibility rather than a biological absence of conservation.

### 3.2 — Retain unambiguous one-base mappings

**Purpose:** remove sites with multiple mappings or altered interval length. Count output rows by original site ID (column 4) and keep IDs occurring once.

```bash
awk 'BEGIN{OFS="\t"} {n[$4]++; line[$4]=$0}
     END{for(id in n) if(n[id]==1) print line[id]}' \
  "$PREP/polar_sites_all.bed" |
  awk 'BEGIN{OFS="\t"} $3-$2==1 {print}' |
  sort -k1,1 -k2,2n > "$PREP/polar_sites_unique.bed"
```

The first `awk` retains only IDs represented by exactly one mapped row. The second requires the lifted interval to remain exactly one base wide. The final `sort` orders the retained sites by polar contig and coordinate.

Check the filtered file immediately:

```bash
ls -lh "$PREP/polar_sites_unique.bed"
head "$PREP/polar_sites_unique.bed"

printf 'Retained unique one-base sites: '
wc -l < "$PREP/polar_sites_unique.bed"

printf 'Duplicate IDs remaining: '
cut -f4 "$PREP/polar_sites_unique.bed" | sort | uniq -d | wc -l

printf 'Non-one-base intervals remaining: '
awk '$3-$2!=1 {n++} END{print n+0}' "$PREP/polar_sites_unique.bed"
```

**Expected:** the last two checks both print `0`. The retained count should be no greater than the number of distinct mapped input IDs.

### 3.3 — Measure mapping losses

**Purpose:** quantify how many VCF sites remain usable before looking at GERP scores.

```bash
wc -l "$PREP/Apennine_sites.bed" "$PREP/polar_sites_all.bed" \
  "$PREP/polar_sites_unique.bed"
head "$PREP/polar_sites_unique.bed"
```

**Expected:** columns 1–3 now contain **polar** coordinates; column 4 still contains the original Apennine `Scaffold_25:position` ID. The unique count is no greater than the input count.

**Check:** inspect `unmapped.bed`, count ambiguous and length-changing loci, and report the uniquely mapped fraction. Unmapped sites have *missing* scores, not score zero.

## Step 4 — Read GERP scores at the polar positions

**Purpose:** extract bigWig scores only from polar regions reached by the sites, then overlap them with the one-base polar BED. bedGraph intervals already carry scores; `bedtools makewindows -w 1` is unnecessary.

**Input:** verified polar-coordinate GERP bigWig and unique polar sites. **The full approximately 7 GB bigWig is needed only once, for instructor preparation.** Students need not download or copy it: the instructor can run Steps 4–5 once and distribute the much smaller `GERP_on_Apennine_unique.bed.gz` (with its index), then students start at Step 6. The commands remain here to show exactly how that teaching file was made.

### 4.1 — Find the polar regions to query

**Purpose:** define one score-extraction span per polar contig reached by the SNPs, avoiding whole-genome conversion. If running this substep, point `GERP_BW` to the instructor's local copy; do not download a copy for every student.

```bash
ls -lh "$GERP_BW"
awk 'BEGIN{OFS="\t"}
     {if(!($1 in min) || $2<min[$1]) min[$1]=$2;
      if(!($1 in max) || $3>max[$1]) max[$1]=$3}
     END{for(c in min) print c,min[c],max[c]}' \
  "$PREP/polar_sites_unique.bed" > "$PREP/polar_regions.tsv"
cat "$PREP/polar_regions.tsv"
```

This first reduces the mapped positions to one covered span per polar contig. It does **not** yet extract scores. **Expected:** rows with a polar contig, minimum start, and maximum end. If the table is empty, return to Step 3.

### 4.2 — Extract the GERP score intervals

**Purpose:** convert only those polar spans from bigWig to scored bedGraph intervals.

```bash
while IFS=$'\t' read -r chr start end; do
  bigWigToBedGraph -chrom="$chr" -start="$start" -end="$end" \
    "$GERP_BW" "$PREP/GERP_${chr}.bedGraph"
done < "$PREP/polar_regions.tsv"
```

**Expected:** one bedGraph file per reached polar contig. Each row has `chrom start end score`; the intervals need not be one base wide. **Check:** an unknown-chromosome error means the bigWig and polar FASTA names or assemblies do not match.

### 4.3 — Join mapped SNPs to GERP scores

**Purpose:** collect the interval scores and match them to the uniquely mapped, one-base SNP positions.

```bash
cat "$PREP"/GERP_*.bedGraph |
  sort -k1,1 -k2,2n > "$PREP/GERP_polar_regions.bedGraph"
bedtools intersect \
  -a "$PREP/polar_sites_unique.bed" \
  -b "$PREP/GERP_polar_regions.bedGraph" \
  -wa -wb > "$PREP/sites_with_polar_scores.tsv"
```

### 4.4 — Inspect the scored join

**Purpose:** verify the original Apennine site ID and polar GERP score are both present before restoring coordinates.

```bash
head "$PREP/sites_with_polar_scores.tsv"
```

**Expected:** eight columns: the four mapped-site columns followed by four polar bedGraph columns. Column 4 is the original Apennine site ID and column 9 the GERP score.

**Check:** `bigWigToBedGraph` reporting an unknown chromosome indicates a likely assembly or contig-name mismatch. Missing values are not zero. Check for sites matching more than one score interval.

## Step 5 — Restore Apennine coordinates and validate

**Purpose:** write a four-column, one-base score table in the coordinates used by SnpEff and GenoLoader. Retain only IDs with one numeric score match.

**Input:** the polar-site/score join from Step 4.

### 5.1 — Restore original Apennine coordinates

**Purpose:** use each stored site ID to return the score to its original VCF coordinate, retaining only a single numeric score per site.

```bash
awk 'BEGIN{OFS="\t"}
     $8 ~ /^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/ {
       n[$4]++; score[$4]=$9
     }
     END {for(id in n) if(n[id]==1) {
       split(id,a,":"); print a[1],a[2]-1,a[2],score[id]
     }}' "$PREP/sites_with_polar_scores.tsv" |
  sort -k1,1 -k2,2n > "$PREP/GERP_on_Apennine_unique.bed"
```

```bash
head "$PREP/GERP_on_Apennine_unique.bed"
```

Column 4 of the scored join still holds the **original** Apennine `Scaffold_25:position` ID. This command converts that ID back to a one-base Apennine BED row and keeps the polar GERP score in column 4. It rejects non-numeric scores and IDs with more than one score match.

### 5.2 — Save an indexed score track

**Purpose:** make the Apennine-coordinate score table reusable in another session.

```bash
GERP="$OUTDIR/GERP_on_Apennine_unique.bed.gz"
```

```bash
bgzip -c "$PREP/GERP_on_Apennine_unique.bed" > "$GERP"
tabix -f -p bed "$GERP"
```

### 5.3 — Check the final score track

**Purpose:** catch scaffold or BED-width mistakes and report how many sites retained a unique score.

```bash
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

If the instructor supplied the precomputed score track and you skipped Steps 1–5, restore the required variables first:

```bash
OUTDIR="$COURSE_DIR/results/gerp"
GT="$COURSE_DIR/results/snpeff/Bears_4pops_s25.ann.POP_OUT.gt"
GERP="$OUTDIR/GERP_on_Apennine_unique.bed.gz"
mkdir -p "$OUTDIR"
```

### 6.1 — Convert GenoLoader sites to BED

**Purpose:** put annotated GenoLoader sites into the same one-base Apennine coordinate system as the GERP track.

```bash
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $2~/^[0-9]+$/ {print $1,$2-1,$2,$3,$4,$5}' "$GT" \
  > "$OUTDIR/genoloader_sites.bed"
head "$OUTDIR/genoloader_sites.bed"
```

This converts the GenoLoader table to one-base Apennine BED rows, keeping its site categories and flags after the coordinates. Check that the first column is `Scaffold_25` before joining.

### 6.2 — Join annotations and scores

**Purpose:** append the positional GERP score to each annotated site with a score available.

```bash
bedtools intersect \
  -a "$OUTDIR/genoloader_sites.bed" -b "$GERP" \
  -wa -wb > "$OUTDIR/sites_with_gerp.tsv"
```

`-wa -wb` writes both the GenoLoader row and the matching GERP row, so the **last column** is the score.

### 6.3 — Inspect coverage and examples

**Purpose:** check whether the join retained a plausible number of sites without duplicating rows.

```bash
wc -l "$OUTDIR/genoloader_sites.bed" "$OUTDIR/sites_with_gerp.tsv"
head "$OUTDIR/sites_with_gerp.tsv"
```

**Expected:** joined rows with the GERP score in the last column. If scores are unique, the joined count cannot exceed the input site count.

**Check:** a larger joined count signals duplicate/overlapping score intervals. Ask whether score availability differs among HIGH, missense, and synonymous categories. The same positional score applies to ABB and SBB at one site; their derived dosages may differ.

## Step 7 — Plot score distributions and derived-genotype contributions

**Purpose:** first examine the GERP scores independently of population genotypes, then combine positive GERP scores with each bear's derived-allele dosage. These are descriptive summaries for Scaffold_25, not estimates of realized fitness load.

**Input:** the Apennine-coordinate GERP track, the polarized GenoLoader table, and the ABB/SBB sample lists.

### 7.1 — Run the plotting script

**Purpose:** produce a per-sample table and two separate PDF figures.

```bash
ABB_LIST="$COURSE_DIR/data/ABB.samples"
SBB_LIST="$COURSE_DIR/data/SBB.samples"

python scripts/plot_gerp_genotypes.py \
  "$GERP" \
  "$GT" \
  "$ABB_LIST" \
  "$SBB_LIST" \
  "$OUTDIR/GERP_derived_scores_by_sample.tsv" \
  "$OUTDIR/GERP_ABB_SBB_Scaffold_25"
```

The script keeps GenoLoader rows whose polarization flag begins with `unfolded` or is `allFix`, matches them to unique GERP-scored positions, and separates the SnpEff classes `MODIFIER`, `LOW`, `MODERATE`, and `HIGH`.

For the genotype-weighted summaries, a positive GERP score is treated as evidence of constraint:

- a heterozygous derived genotype contributes `1 × GERP`;
- a homozygous derived genotype contributes `2 × GERP`;
- total contribution is the sum of the heterozygous and homozygous-derived components.

Negative scores remain visible in the distribution plot but are set to zero in the genotype-weighted summary. This prevents unconstrained scores from cancelling positive constraint, but it is an explicit teaching choice rather than a universal definition of genetic load.

### 7.2 — Inspect the outputs

**Purpose:** confirm that every expected table and figure was created before interpretation.

```bash
ls -lh \
  "$OUTDIR/GERP_derived_scores_by_sample.tsv" \
  "$OUTDIR/GERP_ABB_SBB_Scaffold_25_score_distribution.pdf" \
  "$OUTDIR/GERP_ABB_SBB_Scaffold_25_derived_genotype_scores.pdf"

head "$OUTDIR/GERP_derived_scores_by_sample.tsv"
```

**Expected outputs:**

- `GERP_derived_scores_by_sample.tsv` contains population, individual, impact class, number of called/scored sites, derived-copy count, heterozygous contribution, homozygous-derived contribution, total contribution, and total per called/scored site;
- `..._score_distribution.pdf` shows the raw GERP distribution at retained SNPs and compares scores among SnpEff impact classes;
- `..._derived_genotype_scores.pdf` shows individual ABB and SBB values for heterozygous, homozygous-derived, and total positive-GERP contributions. Points are individuals and horizontal lines are population means.

### 7.3 — Interpret the figures carefully

**Purpose:** distinguish evolutionary constraint from predicted molecular consequence and from genetic load.

A high positive GERP score means that a position is more conserved across the mammal alignment than expected under neutrality. It does **not** prove that the derived bear allele is deleterious. SnpEff impact and GERP answer different questions: SnpEff predicts how a variant changes an annotated transcript, while GERP measures long-term evolutionary constraint at that genomic position.

In the derived-genotype figure, compare both the populations and the two genotype components. A larger homozygous-derived contribution can be especially informative in a small, inbred population because recessive alleles are more often exposed in homozygous form. However, the plotted value does not incorporate selection coefficients, dominance, gene expression, or actual fitness. It is best called a **constraint-weighted derived-allele summary**, not realized genetic load.

**Check:** population differences can also arise from unequal numbers of called/scored sites. Use the `called_scored_sites` and `total_per_called_scored_site` columns before interpreting raw totals. Results from one scaffold are an illustration and must not be generalized automatically to the whole genome.

## Stop and discuss

1. Do GERP distributions differ among MODIFIER, LOW, MODERATE, and HIGH sites? Should they necessarily follow the SnpEff ranking?
2. Does an ABB–SBB contrast arise mainly from heterozygous or homozygous-derived contributions?
3. Why are negative scores retained in the distribution but not subtracted from the positive-GERP derived summary?
4. Why are unscored or ambiguously mapped positions not assigned score zero?
5. How could chain direction, assembly mismatch, or a one-base BED error produce convincing but wrong results?
6. Why should these figures not be described as direct estimates of realized genetic load?

Record your interpretation in the [answer sheet](../answers/student_answers.md).
