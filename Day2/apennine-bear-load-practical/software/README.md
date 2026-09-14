# Instructor software and data preparation

Complete this outside the timed practical. Test every command on the teaching compute node. The sample VCF and SnpEff database are not supplied in this draft.

## One Conda environment for shared dependencies

~~~bash
conda env create -f environment.yml
conda activate bear-load-practical
~~~

The environment provides BCFtools, BEDTools, SAMtools, Java, a C++ compiler, Python, wget, and unzip. It does not automatically provide the bear-specific SnpEff database, GenoLoader binary, or liftOver chain.

## SnpEff

For a **fresh installation**, download and unpack SnpEff from the course directory (do not unpack over an existing customized installation):

~~~bash
cd software
wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd ..
ls -lh software/snpEff/snpEff.jar software/snpEff/snpEff.config
~~~

The ZIP extracts to **software/snpEff/**. Because `latest` can change, record the downloaded SnpEff version for reproducibility. The custom bear database still needs to be prepared and built as described below.

The existing Jarvis analysis used the custom database ID **UrArMar_mUrsArc2** with the frozen BRAKER3/TSEBRA annotation. The [SnpEff lesson](../lessons/01-snpeff.md) now contains the reference/GFF3 inspection and file-preparation commands. In particular, resolve the mUrsArc2 database-label versus mUrsArc1.1 FASTA-name discrepancy before building anything.

For the following instructor-only commands, use the prepared, writable SnpEff directory and the same database ID as the lesson:

~~~bash
SNPEFF_HOME=software/snpEff
DB=UrArMar_mUrsArc2
~~~

The previous workflow also copied **UrArMar.braker3.tsebra.gtf** to **genes.gtf.gz**. That file is optional when building with **-gff3**; SnpEff uses **genes.gff.gz** for the explicit GFF3 build. Keep GTF only if you plan to compare builds. The [official database documentation](https://pcingola.github.io/SnpEff/snpeff/build_db_gff_gtf/) describes the expected filenames and notes that GTF is generally preferred when both formats are valid.

Add these entries **once** to **$SNPEFF_HOME/snpEff.config** (inspect the file first to avoid duplicate definitions):

~~~text
UrArMar_mUrsArc2.genome : Ursus arctos marsicanus
UrArMar_mUrsArc2.codonTable : Standard
~~~

Build the database outside the timed practical:

~~~bash
java -Xmx4g -jar "$SNPEFF_HOME/snpEff.jar" build \
  -gff3 -c "$SNPEFF_HOME/snpEff.config" \
  -nodownload -dataDir "$SNPEFF_HOME/data" \
  -noCheckCds -noCheckProtein -v "$DB" \
  > "$SNPEFF_HOME/${DB}.build.log" 2>&1
~~~

The last two flags reproduce the supplied run but **skip CDS and protein validation**. Check the build log for errors and unexpected gene/transcript counts, and run independent reference/annotation checks; a successful exit alone is insufficient. If CDS and protein FASTAs are available, rebuild with checks enabled.

The historical whole-dataset annotation was:

~~~bash
VCF=/jarvis/scratch/usr/biello/bear/snpeff/VCF/marpolblk.sorted.merged.alignable.final.SNP.vcf.gz
java -Xmx4g -jar "$SNPEFF_HOME/snpEff.jar" \
  -c "$SNPEFF_HOME/snpEff.config" "$DB" "$VCF" |
  gzip -c > marpolblk.sorted.merged.alignable.final.SNP.snpeff.vcf.gz
~~~

For class, distribute a small indexed four-species VCF and make the prepared SnpEff installation available at **software/snpEff/** (or adjust the paths in Part 1). Record the exact DB build, source VCF, annotation rate, and software version. Installation and the database build are instructor tasks; students annotate only the small teaching VCF.

## GenoLoader

The [GenoLoader repository](https://github.com/emitruc/genoloader) documents a C++ implementation of the VCF-to-dosage step:

~~~bash
git clone https://github.com/emitruc/genoloader.git software/genoloader
g++ -O2 -std=c++17 \
  -o software/genoloader/genoloader \
  software/genoloader/GenoLoader.v3.2.cpp
~~~

Test the executable on a small, uncompressed, SnpEff-annotated **biallelic SNP** VCF with GT and ANN. Confirm its output filename and flag meanings against the installed version; GenoLoader's documented CLI and outputs may evolve. The Python notebook provides additional counting and plotting functions, but the supplied Day 2 summary script uses only Python's standard library.

## GERP liftOver — optional instructor preparation

First confirm the **source assembly of the published GERP score track**, its score sign, and whether it is a one-base bedGraph, bigWig, or another format. The fact that polar bear is one of the aligned species is not sufficient to identify the score coordinate reference.

The [GERP lesson](../lessons/03-gerp.md) gives the full **Scaffold_25-only** workflow: align the assemblies once, build and validate an Apennine-to-score-assembly chain, lift SNP positions, extract scores, and restore Apennine coordinates by stable site ID. It avoids converting or lifting the full-genome score track. The instructor prepares the chain and may provide scored checkpoints before the timed practical. This requires minimap2, Transanno, UCSC `liftOver`, and `bigWigToBedGraph` in addition to the Conda tools; check versions and chain direction on a few loci. Reject unmapped, multi-mapped, length-changing, and multiply scored sites, and document the fraction retained. UCSC documents [chain orientation](https://genome.ucsc.edu/goldenpath/help/chain.html) and [zero-based, half-open bedGraph coordinates](https://genome.ucsc.edu/goldenpath/help/bedgraph).
