# Day 2 software and downloads

Run these commands from the Day 2 directory containing `environment.yml`, `software/`, `data/`, and `lessons/`.

## 1. Install the Conda environment

For a new installation:

```bash
conda env create -f software/environment.yml
conda activate bear-load-practical
```

If the environment already exists:

```bash
conda env update -n bear-load-practical -f software/environment.yml --prune
conda activate bear-load-practical
```

This installs:

- BCFtools, BEDTools, and SAMtools;
- Python with Matplotlib, Java, Git, `wget`, and `unzip`;
- a C++ compiler;
- Minimap2 and Transanno;
- UCSC `liftOver`, `bigWigInfo`, and `bigWigToBedGraph`.

Transanno is installed directly from [Bioconda](https://bioconda.github.io/recipes/transanno/README.html) by `environment.yml`; it does not require a separate manual download. If Conda is unavailable, releases and source-build instructions are provided by the [official Transanno repository](https://github.com/informationsea/transanno).

Create the paths used in the lessons:

```bash
COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$COURSE_DIR/software/genoloader:$PATH"
```

Check the commands:

```bash
command -v bcftools bedtools samtools python java
command -v minimap2 transanno liftOver bigWigInfo bigWigToBedGraph
```

## 2. Download SnpEff

```bash
cd software
wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd ..

ls -lh software/snpEff/snpEff.jar software/snpEff/snpEff.config
```

The instructor must also prepare the custom `mUrsArc1.1` database from:

```text
data/mUrsArc1.1.annotation.s25.gff3
data/mUrsArc1.1.genome.s25.fasta
```

Add this block to `software/snpEff/snpEff.config` before building the database:

```text
#---
# Non-standard Databases
#---

# Ursus arctos marsicanus genome, version mUrsArc1.1
mUrsArc1.1.genome : Ursus arctos marsicanus
mUrsArc1.1.codonTable : Standard
```

Students should receive the config file already modified and the two prepared database inputs under:

```text
software/snpEff/data/mUrsArc1.1/genes.gff.gz
software/snpEff/data/mUrsArc1.1/sequences.fa.gz
```

They build `snpEffectPredictor.bin` in their own course directory by following the [SnpEff lesson](../lessons/01-snpeff.md).

## 3. Download and compile GenoLoader

```bash
git clone https://github.com/emitruc/genoloader.git software/genoloader

g++ -O2 -std=c++17 \
  -o software/genoloader/genoloader \
  software/genoloader/GenoLoader.v3.2.cpp
```

Check the filename of the C++ source after cloning because a later GenoLoader release may use a different filename.

## 4. Files required for the optional GERP lesson

The GERP analysis requires more than the four command-line programs. The instructor needs:

1. The Ensembl release 114 polar-bear GERP score file:

   [gerp_conservation_scores.ursus_maritimus.UrsMar_1.0.bw](https://ftp.ensembl.org/pub/release-114/compara/conservation_scores/91_mammals.gerp_conservation_score)

2. The **Ensembl release 114 UrsMar_1.0 toplevel FASTA** matching the `AVOR...` sequence names in that bigWig, from the [polar-bear DNA directory](https://ftp.ensembl.org/pub/release-114/fasta/ursus_maritimus/dna/):

   ```text
   data/Ursus_maritimus.UrsMar_1.0.dna.toplevel.fa
   ```

3. The prepared Apennine Scaffold_25 FASTA from the reference used to call the teaching VCF:

   ```text
   data/mUrsArc1.1.genome.s25.fasta
   ```
4. The indexed Scaffold_25 teaching VCF:

   ```text
   data/Bears_4pops_s25.vcf.gz
   data/Bears_4pops_s25.vcf.gz.csi
   ```

5. Minimap2, Transanno, UCSC `liftOver`, `bigWigInfo`, and `bigWigToBedGraph`, already installed through `environment.yml`.

The bigWig is approximately 7 GB and should be downloaded only once. Students do not need their own copy if the instructor prepares and distributes:

```text
results/gerp/GERP_on_Apennine_unique.bed.gz
results/gerp/GERP_on_Apennine_unique.bed.gz.tbi
```

This smaller file contains only the GERP scores successfully mapped back to the Apennine Scaffold_25 positions. Students can then begin at Step 6 of the optional GERP lesson.

## 5. Final installation check

```bash
conda activate bear-load-practical
COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$COURSE_DIR/software/genoloader:$PATH"

command -v bcftools bedtools samtools minimap2 transanno liftOver bigWigInfo bigWigToBedGraph
java -version
ls -lh software/snpEff/snpEff.jar
ls -lh software/genoloader/genoloader
```
