# Instructor software and data preparation

Complete this outside the timed practical. Test every command on the teaching compute node. The sample VCF and SnpEff database are not supplied in this draft.

## One Conda environment for shared dependencies

~~~bash
conda env create -f environment.yml
conda activate bear-load-practical
~~~

The environment provides BCFtools, BEDTools, Java, a C++ compiler, and Python. It does not automatically provide the bear-specific SnpEff database, GenoLoader binary, or liftOver chain.

## SnpEff

Place a tested SnpEff distribution under **software/snpEff/** so that **snpEff.jar** and **snpEff.config** match the commands in Part 1. Build or obtain a database for the **exact Apennine reference assembly and gene annotation** used for the VCF. Record database ID, assembly accession, annotation source, and version. Follow the [official custom-database instructions](https://pcingola.github.io/SnpEff/snpeff/build_db/). Do not use a polar-bear or another brown-bear assembly merely because its database is available.

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

For one-base source BED records, a validated source-to-Apennine chain permits liftOver of coordinates. Keep a unique site identifier and score through transfer; reject unmapped, multi-mapped, length-changing, and non-1:1 sites. Check chromosome names and a sample of mapped positions against both assemblies, and document the fraction retained. Do not reuse the score track as if liftOver regenerated a 92-mammal alignment; it transfers existing scores to new coordinates. The student extension expects a verified four-column file **data/GERP_on_Apennine_unique.bed.gz**. UCSC documents [chain orientation](https://genome.ucsc.edu/goldenpath/help/chain.html) and [zero-based, half-open bedGraph coordinates](https://genome.ucsc.edu/goldenpath/help/bedgraph).
