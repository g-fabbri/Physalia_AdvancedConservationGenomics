# Part 1 — What might each SNP do?

Estimated terminal time: 20 minutes with a prebuilt database; allow extra time if students prepare database files themselves.

SnpEff compares each VCF allele with gene models and writes predicted consequences to the VCF `ANN` field. The same SNP can have multiple transcript annotations. `HIGH`, `MODERATE`, `LOW`, and `MODIFIER` are **predicted impact categories**, not measurements of selection coefficients or fitness.

The **UrArMar_mUrsArc2** database uses the frozen Apennine-bear BRAKER3/TSEBRA gene annotation and a reference FASTA. We first show how those biological inputs are assembled. The software itself should already be available; its installation is described separately in [software setup](../software/README.md).

## Step 1 — Prepare the bear reference and gene annotation

**Purpose:** understand exactly which genome sequence and gene models SnpEff uses to predict variant consequences. This is data preparation, not software installation. The instructor can demonstrate it or provide the prepared database as a checkpoint if copying a whole genome would take too long in class.

**Input:** the frozen GFF3 and the **verified reference FASTA used to call the teaching VCF**. The original recipe called the database mUrsArc2 but copied a FASTA named mUrsArc1.1; those names must be reconciled before entering a real path below. A matching scaffold name alone does not prove that two assemblies are identical.

From your **private course directory**, set:

~~~bash
SNPEFF_HOME="$PWD/software/snpEff"
DB=UrArMar_mUrsArc2
GFF=/jarvis/scratch/usr/biello/bear/annotation/annotation_versions/frozen/UrArMar.braker3.tsebra.gff3
REF_FASTA=/path/to/verified/VCF_reference.fasta

bcftools view -h "$VCF" | grep '^##contig' | head
grep '^>' "$REF_FASTA" | head
awk '$0 !~ /^#/ {print $1; if (++n==5) exit}' "$GFF"
~~~

**Expected:** scaffold IDs from the VCF, FASTA, and GFF3 that can be matched to the same assembly. The instructor must also check contig lengths and sampled VCF REF alleles against the FASTA; this short inspection is only an orientation check.

After the reference has been verified, place files where SnpEff expects them:

~~~bash
mkdir -p "$SNPEFF_HOME/data/$DB"
cp "$GFF" "$SNPEFF_HOME/data/$DB/genes.gff"
cp "$REF_FASTA" "$SNPEFF_HOME/data/$DB/sequences.fa"
gzip "$SNPEFF_HOME/data/$DB/genes.gff"
gzip "$SNPEFF_HOME/data/$DB/sequences.fa"
~~~

**Expected:** **genes.gff.gz** and **sequences.fa.gz** under **software/snpEff/data/UrArMar_mUrsArc2/**. The previous Jarvis preparation also kept a GTF, but the supplied build used **-gff3**, so that build reads the GFF3. These files prepare the database inputs; the instructor's [database-build instructions](../software/README.md) explain the config entry and build command. Do not rebuild a shared installation from every student account.

**Check:**

~~~bash
ls -lh "$SNPEFF_HOME/data/$DB/genes.gff.gz" \
       "$SNPEFF_HOME/data/$DB/sequences.fa.gz"
~~~

## Step 2 — Annotate the teaching VCF

**Purpose:** add consequences on the same reference assembly used for variant calling.

**Input:** the teaching VCF and the prebuilt, instructor-verified **UrArMar_mUrsArc2** SnpEff database.

~~~bash
SNPEFF_DB="$DB"
SNPEFF_CONFIG="$SNPEFF_HOME/snpEff.config"
SNPEFF_JAR="$SNPEFF_HOME/snpEff.jar"
ANNOTATED="$OUTDIR/Bears_4pops_s25.ann.vcf"

java -Xmx4g -jar "$SNPEFF_JAR" \
  -c "$SNPEFF_CONFIG" \
  -dataDir "$SNPEFF_HOME/data" \
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

## Step 3 — Compare consequence classes

**Purpose:** see what the annotations actually contain before treating any class as putatively deleterious.

~~~bash
bcftools query -f '%INFO/ANN\n' "$ANNOTATED" | \
  awk -F'[,|]' '$1!="." {impact[$3]++} END {for (x in impact) print x, impact[x]}' | \
  sort
~~~

**Expected:** counts labelled `HIGH`, `MODERATE`, `LOW`, and/or `MODIFIER`. This quick check uses the **first** transcript annotation per VCF record; it is not a complete transcript-aware summary.

**Check:** record which categories are common. Discuss how transcript choice and reference gene-model quality could change the labels. Continue to [polarization and burden](02-genoloader.md).
