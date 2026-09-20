# Part 1 — What might each SNP do?

Estimated terminal time: 20 minutes with a prebuilt database; allow extra time if students prepare database files themselves.

## Start your terminal

From the Day 2 directory containing `data/`, `software/`, and `results/`, run this in each new terminal:

```bash
conda activate bear-load-practical
COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$COURSE_DIR/software/genoloader:$PATH"
```

Conda supplies the shared tools; standalone programs are kept under Day 2 `software/`. SnpEff is called through its JAR at `software/snpEff/snpEff.jar`, so it does not need a PATH entry. The variables from Part 0 must also be set if you opened a new terminal.

SnpEff compares each VCF allele with gene models and writes predicted consequences to the VCF `ANN` field. The same SNP can have multiple transcript annotations. `HIGH`, `MODERATE`, `LOW`, and `MODIFIER` are **predicted impact categories**, not measurements of selection coefficients or fitness.

The **mUrsArc1.1** database uses the Apennine-bear gene annotation and the matching reference FASTA available from [Zenodo](https://zenodo.org/records/15349716). For this exercise, both files have already been restricted to Scaffold_25 and given consistent names. We first show how those biological inputs are assembled. The software itself should already be available; its installation is described separately in [software setup](../software/README.md).

## Step 1 — Prepare the bear reference and gene annotation

**Purpose:** understand exactly which genome sequence and gene models SnpEff uses to predict variant consequences and prepare them in the structure expected by SnpEff.

**Input:** the Scaffold_25 GFF3 and the **reference FASTA used to call the teaching VCF**.

### 1.1 — Identify and compare the reference inputs

**Purpose:** verify that the VCF, gene annotation, and FASTA use the same assembly before building the custom database.

Set the definitive database name and input paths:

~~~bash
SNPEFF_HOME="$COURSE_DIR/software/snpEff"
DB=mUrsArc1.1
GFF="$COURSE_DIR/data/mUrsArc1.1.annotation.s25.gff3"
REF_FASTA="$COURSE_DIR/data/mUrsArc1.1.genome.s25.fasta"
~~~

Inspect the assembly identifiers independently:

~~~bash
bcftools query -f '%CHROM\n' "$VCF" | sort -u
~~~

~~~bash
grep '^>' "$REF_FASTA" | head
~~~

~~~bash
awk '$0 !~ /^#/ {print $1; if (++n==5) exit}' "$GFF"
~~~

**Expected:** the VCF, FASTA, and GFF3 all use **Scaffold_25** from the same assembly. Matching names are necessary, although database provenance and coordinate compatibility must also have been verified during course preparation.

### 1.2 — Understand the SnpEff configuration entry

**Purpose:** connect the database name used on the command line with the custom genome stored under the SnpEff data directory.

The instructor has already added the following entry to **software/snpEff/snpEff.config**:

~~~text
#---
# Non-standard Databases
#---

# Ursus arctos marsicanus genome, version mUrsArc1.1
mUrsArc1.1.genome : Ursus arctos marsicanus
mUrsArc1.1.codonTable : Standard
~~~

The text before `.genome` is the database identifier. It must match both the value of **DB** and the directory name under **software/snpEff/data/**. The `.genome` line supplies a human-readable description, while `.codonTable` tells SnpEff to interpret coding sequences using the standard genetic code.

Students do **not** need to edit the shared config file. Check that the prepared entry is present:

~~~bash
grep '^mUrsArc1.1\.' "$SNPEFF_HOME/snpEff.config"
~~~

**Expected:** the two database properties shown above. If nothing is printed, stop and ask the instructor rather than modifying the shared installation during class.

### 1.3 — Place the verified database inputs

**Purpose:** give SnpEff the GFF3 and FASTA filenames expected by its GFF3 database build. Do this only after the reference has been verified.

~~~bash
mkdir -p "$SNPEFF_HOME/data/$DB"
~~~

~~~bash
cp "$GFF" "$SNPEFF_HOME/data/$DB/genes.gff"
cp "$REF_FASTA" "$SNPEFF_HOME/data/$DB/sequences.fa"
~~~

~~~bash
gzip "$SNPEFF_HOME/data/$DB/genes.gff"
gzip "$SNPEFF_HOME/data/$DB/sequences.fa"
~~~

**Expected:** **genes.gff.gz** and **sequences.fa.gz** under **software/snpEff/data/mUrsArc1.1/**. The database build uses the GFF3 file. Students receive these files and the built database already prepared; they do not rebuild a shared installation from every account.

### 1.4 — Check the prepared database files

**Purpose:** ensure both compressed inputs exist before the instructor builds or supplies the database.

~~~bash
ls -lh "$SNPEFF_HOME/data/$DB/genes.gff.gz" \
       "$SNPEFF_HOME/data/$DB/sequences.fa.gz"
~~~

## Step 2 — Build the custom SnpEff database

**Purpose:** combine the Scaffold_25 reference sequence and gene models into the binary database that SnpEff uses during annotation.

### 2.1 — Define the SnpEff program and config paths

~~~bash
SNPEFF_JAR="$SNPEFF_HOME/snpEff.jar"
SNPEFF_CONFIG="$SNPEFF_HOME/snpEff.config"
DATA="$SNPEFF_HOME/data"
~~~

**Check:**

~~~bash
ls -lh "$SNPEFF_JAR" "$SNPEFF_CONFIG"
~~~

Both files must exist before continuing.

### 2.2 — Build the database

The `build` command does not annotate the bear variants yet. It prepares a reusable representation of the reference genome and its gene annotation. In broad terms, SnpEff:

1. reads **sequences.fa.gz** and creates the Scaffold_25 reference sequence in its internal genome model;
2. parses **genes.gff.gz** and connects features through their GFF3 IDs and parent relationships;
3. reconstructs the hierarchy of genes, transcripts, exons, CDS segments, and UTRs, including their coordinates and strands;
4. assigns the configured codon table and derives transcript/coding information from the reference sequence;
5. checks the feature structure and, unless disabled, can compare reconstructed CDS and proteins with independent validation FASTAs;
6. serializes the resulting predictor so it can be loaded quickly when variants are annotated.

The essential point is that SnpEff is compiling a coordinate-aware gene model. Later, when it receives a VCF position and allele change, it can determine which genomic features overlap that position and predict how the alternate allele changes each affected transcript.

~~~bash
java -Xmx4g -jar "$SNPEFF_JAR" build \
  -gff3 \
  -c "$SNPEFF_CONFIG" \
  -nodownload \
  -dataDir "$DATA" \
  -noCheckCds \
  -noCheckProtein \
  -v \
  "$DB"
~~~

| Argument | Meaning |
|---|---|
| **-Xmx4g** | Allow Java to use up to 4 GB of memory |
| **build** | Build a SnpEff genome database rather than annotate a VCF |
| **-gff3** | Read gene models from **genes.gff.gz** as GFF3 |
| **-c** | Use the prepared course config file |
| **-nodownload** | Use local files and do not try to download a public database |
| **-dataDir** | Use the database directory inside the Day 2 course folder |
| **-noCheckCds** | Skip comparison against a separate CDS validation file |
| **-noCheckProtein** | Skip comparison against a separate protein validation file |
| **-v** | Print detailed progress messages |
| **mUrsArc1.1** | Database identifier; it must match the config entry and directory name |

The two `-noCheck...` options are used because the course supplies the reference FASTA and GFF3 but not independent CDS and protein FASTA files. They allow the build to proceed, but they also remove two useful validation checks. In a production annotation workflow, provide matching CDS and protein sequences and keep those checks enabled whenever possible.

**Expected terminal output:** SnpEff reports the config and database paths, reads the GFF3 and reference sequence, and prints summaries of the genes, transcripts, exons, and coding features found on Scaffold_25. It then announces that it is saving the database. The command should return to the prompt without a Java exception or fatal error.

Warnings should be read rather than ignored, particularly warnings about missing chromosome sequences, inconsistent scaffold names, features outside the reference boundaries, missing parents, or malformed transcripts. The build can sometimes finish while dropping or modifying problematic features, so “a file was created” is necessary but not sufficient evidence of a high-quality annotation database.

The database directory contains both source inputs and compiled output:

| File | Created when? | What it contains |
|---|---|---|
| **genes.gff.gz** | Prepared in Step 1 | Source gene annotation: genomic features, coordinates, strands, IDs, and parent-child relationships |
| **sequences.fa.gz** | Prepared in Step 1 | Source reference nucleotide sequence for Scaffold_25 |
| **snpEffectPredictor.bin** | Created by `build` | SnpEff's serialized predictor: the parsed genome and feature hierarchy used to locate genes/transcripts and predict allele effects efficiently |
| **sequence.Scaffold_25.bin** | Created by `build` | Serialized Scaffold_25 nucleotide sequence. SnpEff uses it to recover reference bases, construct transcript/CDS sequences, and evaluate codon and amino-acid changes |
| **Other `.bin` files, if present** | Version- or option-dependent | Auxiliary regulation, motif, or other feature data; they are not expected from every simple GFF3 build |

Recent SnpEff versions store reference sequence in chromosome-specific binary files by default. Therefore, the approximately 1 MB predictor and the larger Scaffold_25 sequence file serve different purposes: the predictor describes **where and how genomic features are organized**, while the sequence file provides the **nucleotides underlying those features**.

The build does **not** create an annotated VCF, a table of variants, or genetic-load results. Those are produced only when the compiled predictor is applied to the teaching VCF in Step 3.

**Check:** confirm that the compiled database file was created:

~~~bash
ls -lh "$SNPEFF_HOME/data/$DB/snpEffectPredictor.bin" \
  "$SNPEFF_HOME/data/$DB/sequence.${CHROM}.bin"
find "$SNPEFF_HOME/data/$DB" -maxdepth 1 -type f -printf '%f\n' | sort
~~~

Both binary files should exist and have nonzero sizes. The file listing should also show the two compressed source inputs. If either `snpEffectPredictor.bin` or `sequence.Scaffold_25.bin` is missing, inspect the final build messages before attempting annotation.


## Step 3 — Annotate the teaching VCF

**Purpose:** predict and record how each alternate allele may affect annotated genes and transcripts—for example, whether it is synonymous, missense, stop-gained, intronic, or intergenic. The SnpEff database must use the same reference assembly and coordinates as the VCF.

**Input:** the teaching VCF and the **mUrsArc1.1** SnpEff database built in Step 2.

During annotation, SnpEff processes each VCF record as follows:

1. it confirms that the VCF chromosome and position can be located in the database;
2. it compares the VCF REF allele with the stored Scaffold_25 reference sequence;
3. it finds genes, transcripts, exons, CDS regions, splice regions, UTRs, or intergenic regions overlapping the variant;
4. for coding variants, it reconstructs the affected codon and predicts any amino-acid change;
5. it creates one or more transcript-specific annotations and adds them to the VCF INFO field named **ANN**.

SnpEff annotates variants; it does not change the sample genotypes, infer ancestral alleles, calculate allele frequencies, or measure genetic load. The original VCF columns and sample data remain present in the output.

### 3.1 — Add predicted functional annotations

**Purpose:** write an ANN-annotated, uncompressed VCF that GenoLoader can read.

~~~bash
SNPEFF_DB="$DB"
ANNOTATED="$OUTDIR/Bears_4pops_s25.ann.vcf"

java -Xmx4g -jar "$SNPEFF_JAR" \
  -c "$SNPEFF_CONFIG" \
  -dataDir "$SNPEFF_HOME/data" \
  -noStats \
  "$SNPEFF_DB" "$VCF" > "$ANNOTATED"
~~~

| Argument | Meaning |
|---|---|
| **-Xmx4g** | Allow Java to use up to 4 GB of memory |
| **-c** | Read database definitions from the prepared course config |
| **-dataDir** | Load the compiled database from the Day 2 course directory |
| **-noStats** | Do not create the optional SnpEff HTML and gene-summary statistics reports |
| **mUrsArc1.1** | Load the database built in Step 2 |
| **$VCF** | Read the original four-population Scaffold_25 variants and genotypes |
| **> $ANNOTATED** | Redirect the annotated VCF written to standard output into the selected file |

SnpEff normally writes the annotated VCF to standard output and progress or warning messages to standard error. The `>` redirection captures only the VCF stream. Do not omit it: otherwise the annotated VCF will be printed to the terminal rather than saved.

**Expected terminal output:** progress messages may identify the SnpEff version, config, database, and input file. There should be no fatal error about a missing database, chromosome, or reference mismatch.

**Expected file:** **results/genetic_load/Bears_4pops_s25.ann.vcf**, an uncompressed VCF containing the original records and sample genotypes plus SnpEff header definitions and `ANN` values. GenoLoader's documented C++ command accepts a `.vcf` input; keeping this file uncompressed avoids assuming gzip support.

**Check:** confirm that the file is nonempty and that annotation did not change the number of VCF records:

~~~bash
ls -lh "$ANNOTATED"
printf 'Input records: '; bcftools view -H "$VCF" | wc -l
printf 'Annotated records: '; bcftools view -H "$ANNOTATED" | wc -l
~~~

The two record counts should match. Equal counts do not mean every record has a useful gene annotation: intergenic records can be annotated as intergenic, and records that cannot be matched correctly require investigation.

### 3.2 — Inspect the ANN field

**Purpose:** confirm that the annotated VCF has both the ANN definition and per-site consequences.

~~~bash
bcftools view -h "$ANNOTATED" | grep 'ID=ANN'
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$ANNOTATED" | head -n 3
~~~

The first command should find an `ANN` header. The second should show consequence strings separated by `|`. If annotations are unexpectedly absent, stop and check chromosome names, genome build, and database provenance. GenoLoader skips loci without `ANN`.

The ANN field is comma-separated when one allele has annotations for multiple transcripts. Within each annotation, pipe-separated fields describe items such as:

| ANN component | Meaning |
|---|---|
| **Allele** | Alternate allele being annotated |
| **Annotation** | Sequence Ontology effect term, such as `missense_variant` or `synonymous_variant` |
| **Impact** | Broad predicted category: `HIGH`, `MODERATE`, `LOW`, or `MODIFIER` |
| **Gene name / gene ID** | Affected annotated gene, when applicable |
| **Feature type / feature ID** | Usually the affected transcript and its identifier |
| **Biotype** | Transcript class, such as protein coding, when provided by the GFF3 |
| **Rank** | Exon or intron number within the transcript |
| **HGVS.c / HGVS.p** | Predicted coding-DNA and protein-level change, when applicable |

One VCF variant may therefore receive several ANN entries with different predicted effects because transcripts can use different exons or reading frames. The annotation is a model-based prediction determined by the reference allele, alternate allele, transcript model, and selected database; it is not direct experimental evidence that the variant changes fitness.
## Step 4 — Compare predicted-effect classes

**Purpose:** see what the annotations actually contain before treating any class as putatively deleterious.

### 4.1 — Count first-transcript impact labels

**Purpose:** get a quick orientation to the predicted-impact categories; this is not a transcript-aware final summary.

~~~bash
bcftools query -f '%INFO/ANN\n' "$ANNOTATED" | \
  awk -F'[,|]' '$1!="." {impact[$3]++} END {for (x in impact) print x, impact[x]}' | \
  sort
~~~

**Expected:** counts labelled `HIGH`, `MODERATE`, `LOW`, and/or `MODIFIER`. This quick check uses the **first** transcript annotation per VCF record; it is not a complete transcript-aware summary.

**Check:** record which categories are common. Discuss how transcript choice and reference gene-model quality could change the labels. Continue to [polarization and burden](02-genoloader.md).
