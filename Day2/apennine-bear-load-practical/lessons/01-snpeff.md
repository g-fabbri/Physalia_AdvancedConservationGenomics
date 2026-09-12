# Part 1 — What might each SNP do?

Estimated terminal time: 20 minutes.

SnpEff compares each VCF allele with gene models and writes predicted consequences to the VCF `ANN` field. The same SNP can have multiple transcript annotations. `HIGH`, `MODERATE`, `LOW`, and `MODIFIER` are **predicted impact categories**, not measurements of selection coefficients or fitness.

## Step 1 — Annotate the teaching VCF

**Purpose:** add consequences on the same reference assembly used for variant calling.

**Input:** the teaching VCF and a prebuilt, instructor-verified SnpEff database. Replace `BEAR_DB` with its actual database ID.

~~~bash
SNPEFF_DB=BEAR_DB
SNPEFF_CONFIG=software/snpEff/snpEff.config
SNPEFF_JAR=software/snpEff/snpEff.jar
ANNOTATED="$OUTDIR/four_species.ann.vcf"

java -Xmx4g -jar "$SNPEFF_JAR" \
  -c "$SNPEFF_CONFIG" \
  -noStats \
  "$SNPEFF_DB" "$VCF" > "$ANNOTATED"
~~~

**Expected:** an uncompressed annotated VCF. GenoLoader's documented C++ command accepts a `.vcf` input; the uncompressed output avoids assuming gzip support.

**Check:**

~~~bash
bcftools view -h "$ANNOTATED" | grep 'ID=ANN'
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$ANNOTATED" | head -n 3
~~~

The first command should find an `ANN` header. The second should show consequence strings separated by `|`. If annotations are unexpectedly absent, stop and check chromosome names, genome build, and database provenance. GenoLoader skips loci without `ANN`.

## Step 2 — Compare consequence classes

**Purpose:** see what the annotations actually contain before treating any class as putatively deleterious.

~~~bash
bcftools query -f '%INFO/ANN\n' "$ANNOTATED" | \
  awk -F'[,|]' '$1!="." {impact[$3]++} END {for (x in impact) print x, impact[x]}' | \
  sort
~~~

**Expected:** counts labelled `HIGH`, `MODERATE`, `LOW`, and/or `MODIFIER`. This quick check uses the **first** transcript annotation per VCF record; it is not a complete transcript-aware summary.

**Check:** record which categories are common. Discuss how transcript choice and reference gene-model quality could change the labels. Continue to [polarization and burden](02-genoloader.md).
