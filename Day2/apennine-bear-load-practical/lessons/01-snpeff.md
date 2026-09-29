# Part 1 — What might each SNP do?

In this lesson, we will use **SnpEff** to predict the possible functional effects of SNPs on annotated genes and transcripts.

SnpEff provides many ready-made databases for commonly studied species and reference genomes. When one of these databases matches the assembly used to produce the VCF, it can be downloaded and used directly.

Our Apennine brown bear assembly is not included among the prepared SnpEff databases, so we need a **custom database**. 

To run, SnpEff requires:

- a reference genome sequence;
- a matching gene annotation;
- a configuration entry identifying the custom genome;
- a VCF produced using the same reference assembly.

It is important that the FASTA, annotation and VCF all use the same assembly, scaffold names and coordinate system. Otherwise, SnpEff may assign variants to the wrong features or fail to annotate them.

For this practical, all the required files have already been prepared for **Scaffold 25**. Students will build the custom SnpEff database, but they do not need to modify the configuration or prepare the input files.

## Start your terminal

From your Day 2 directory, run:

```bash
conda activate bear-load-practical

COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$PATH"
```


Define the files and directories used during the lesson:

```bash
VCF="$COURSE_DIR/data/Bears_4pops_s25.vcf.gz"

SNPEFF_HOME="$COURSE_DIR/software/snpEff"
SNPEFF_JAR="$SNPEFF_HOME/snpEff.jar"
SNPEFF_CONFIG="$SNPEFF_HOME/snpEff.config"
SNPEFF_DATA="$SNPEFF_HOME/data"

DB=mUrsArc1.1
CHROM=Scaffold_25

OUTDIR="$COURSE_DIR/results/snpeff"
mkdir -p "$OUTDIR"
```

Here:

- `VCF` is the four-population VCF containing variants from Scaffold 25;
- `SNPEFF_HOME` is the SnpEff installation directory;
- `SNPEFF_JAR` is the SnpEff program;
- `SNPEFF_CONFIG` is the prepared configuration file;
- `SNPEFF_DATA` contains the custom bear database;
- `DB` is the identifier of the custom database;
- `OUTDIR` is where the annotated VCF will be written.

## Step 1 — Understand the custom SnpEff database

To build a custom database, SnpEff needs:

- a **reference FASTA** containing the nucleotide sequence;
- a matching **GFF3 annotation** describing genes, transcripts, exons and coding regions;
- an entry in `snpEff.config` identifying the custom genome;
- a database directory inside the SnpEff `data/` directory.

The reference and annotation used here contain only **Scaffold 25** and use the same assembly and coordinates as the VCF.

### 1.1 — Configuration file (ALREADY PREPARED)

The following entry was added to `software/snpEff/snpEff.config`:

```text
#---
# Non-standard Databases
#---

# Ursus arctos marsicanus genome, version mUrsArc1.1
mUrsArc1.1.genome : Ursus arctos marsicanus
mUrsArc1.1.codonTable : Standard
```

Here:

- `mUrsArc1.1` is the database identifier used by SnpEff;
- `.genome` provides a description of the genome;
- `.codonTable` tells SnpEff to use the standard genetic code.

The database identifier must match the name of the corresponding directory inside `software/snpEff/data/`.

### 1.2 — Database directory (ALREADY PREPARED)

A directory called `data` was created inside the SnpEff installation. Inside it, a database directory was created using the database identifier:

```text
software/snpEff/data/mUrsArc1.1/
```


The prepared database inputs are therefore:

```text
software/snpEff/data/mUrsArc1.1/genes.gff.gz
software/snpEff/data/mUrsArc1.1/sequences.fa.gz
```

The two files have different roles:

- `genes.gff.gz` describes the locations and structures of genes, transcripts, exons and coding regions;
- `sequences.fa.gz` contains the reference nucleotide sequence used to predict codon and amino-acid changes.

These files and the modified configuration are already provided. Students only need to confirm that they are available:

```bash
grep '^mUrsArc1.1\.' "$SNPEFF_CONFIG"
```

```bash
ls -lh \
  "$SNPEFF_DATA/$DB/genes.gff.gz" \
  "$SNPEFF_DATA/$DB/sequences.fa.gz"
```

If the configuration entry and both files are present, the custom database is ready to be built.

## Step 2 — Build the custom SnpEff database

**Purpose:** compile the reference sequence and gene annotation into a database that SnpEff can use efficiently.

```bash
java -Xmx4g -jar "$SNPEFF_JAR" build \
  -gff3 \
  -c "$SNPEFF_CONFIG" \
  -nodownload \
  -dataDir "$SNPEFF_DATA" \
  -noCheckCds \
  -noCheckProtein \
  -v \
  "$DB"
```

The most important options are:

| Option | Meaning |
|---|---|
| `build` | Build a SnpEff database |
| `-gff3` | Read the annotation in GFF3 format |
| `-c` | Use the supplied configuration file |
| `-nodownload` | Use the local course files |
| `-dataDir` | Store and read the database from your results directory |
| `-noCheckCds` | Skip validation against a separate CDS FASTA |
| `-noCheckProtein` | Skip validation against a separate protein FASTA |
| `-v` | Print detailed progress information |
| `-Xmx4g` | Allow Java to use up to 4 GB of memory |

The command creates two important binary files:
- snpEffectPredictor.bin contains the compiled annotation model. It stores the positions, structures and relationships of genes, transcripts, exons, coding regions and other annotated features. SnpEff uses it to determine which features overlap each variant.
- sequence.Scaffold_25.bin contains the Scaffold 25 reference sequence in SnpEff’s internal format. SnpEff uses the nucleotide sequence to reconstruct codons and predict whether a variant causes effects such as a synonymous change, missense change or premature stop codon.
These files allow SnpEff to load the custom database quickly without reading and rebuilding the original FASTA and GFF3 files every time.

Check that they were created:

```bash
ls -lh \
  "$SNPEFF_DATA/$DB/snpEffectPredictor.bin" \
  "$SNPEFF_DATA/$DB/sequence.${CHROM}.bin"
```

## Step 3 — Annotate the VCF

**Purpose:** predict how each alternate allele may affect genes and transcripts.

```bash
ANNOTATED="$OUTDIR/Bears_4pops_s25.ann.vcf"
```

```bash
java -Xmx4g -jar "$SNPEFF_JAR" \
  -c "$SNPEFF_CONFIG" \
  -dataDir "$SNPEFF_DATA" \
  -noStats \
  "$DB" "$VCF" > "$ANNOTATED"
```

The `>` symbol redirects the annotated VCF produced by SnpEff into the file stored in `$ANNOTATED`.

SnpEff preserves the original variants, genotypes and samples. It adds its predictions to the VCF `INFO` field named `ANN`.

Check the output:

```bash
ls -lh "$ANNOTATED"
```

Confirm that the input and output contain the same number of variant records:

```bash
printf 'Input records: '
bcftools view -H "$VCF" | wc -l

printf 'Annotated records: '
bcftools view -H "$ANNOTATED" | wc -l
```

The two counts should be identical.

## Step 4 — Inspect the annotations

Check that the `ANN` field has been added to the VCF header:

```bash
bcftools view -h "$ANNOTATED" | grep 'ID=ANN'
```

Inspect the first three annotated variants:

```bash
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
  "$ANNOTATED" | head -n 3
```

Each `ANN` entry can contain:

- the alternate allele;
- the predicted consequence, such as `missense_variant`;
- the predicted impact category;
- the affected gene and transcript;
- the predicted coding or protein change.

A variant can have several annotations because it may affect multiple transcripts.

The four broad impact categories are:

| Category | General interpretation |
|---|---|
| `HIGH` | Predicted major disruption, such as a stop-gained variant |
| `MODERATE` | Predicted protein-changing effect, such as a missense variant |
| `LOW` | Usually limited protein effect, such as a synonymous variant |
| `MODIFIER` | Usually non-coding or difficult to interpret directly |

These are predicted functional categories. They are **not** direct measurements of fitness or selection.

## Step 5 — Count the predicted-impact categories

**Purpose:** obtain a quick overview of the annotations.

```bash
bcftools query -f '%INFO/ANN\n' "$ANNOTATED" | \
  awk -F'[,|]' '
    $1!="." {
      impact[$3]++
    }
    END {
      for (x in impact)
        print x, impact[x]
    }
  ' | sort
```

This quick summary counts the impact category from the first transcript annotation reported for each variant. It is useful for checking the results, but it is not a complete transcript-aware analysis.

### Question for discussion

Why should we avoid interpreting every `HIGH`-impact variant as a deleterious mutation?

Consider annotation quality, alternative transcripts, genotype state, allele frequency and whether the affected gene is biologically important.
