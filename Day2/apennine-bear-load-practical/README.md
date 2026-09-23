# Day 2 — Genetic load in two bear populations

## The question

Day 1 asked whether Apennine (ABB) and Slovak (SBB) brown bears have different demographic histories. Today we ask whether they carry different **burdens of putatively deleterious derived alleles**. Black bears (BLB) and polar bears (POB) provide outgroup information for deciding which allele is likely ancestral. Functional effects come from SnpEff; GenoLoader polarizes genotypes and records derived-allele dosage. A GERP conservation-score comparison is an optional extension.

An allele annotated `HIGH` or `MODERATE` is **not** a measured fitness effect. A larger derived-allele count is not automatically a larger realized load. Keep prediction, ancestry, genotype state, and fitness distinct throughout.

The timed analysis uses **Scaffold_25 only**, from `data/Bears_4pops_s25.vcf.gz`. All site counts and burdens refer to that scaffold, not to the whole genome.


## Preparation and missing local details

Before class, supply:

- `data/Bears_4pops_s25.vcf.gz` and its `.csi` index; confirm it contains only biallelic SNPs on Scaffold_25 before GenoLoader;
- the four actual lists `data/ABB.samples`, `data/SBB.samples`, `data/BLB.samples`, and `data/POB.samples`;
- the custom **UrArMar_mUrsArc1.1** SnpEff database, built from the matching Scaffold_25 reference and annotation files supplied for the course;
- a precomputed annotation checkpoint for the same Scaffold_25 VCF if live SnpEff annotation takes too long;
- GenoLoader compiled and tested against that VCF;
- for the optional GERP extension: the score track's **source assembly**, score definition and format, a validated Scaffold_25-to-score-assembly chain, and a precomputed site-score checkpoint with mapping QC.

The four sample-list filenames follow the supplied Jarvis `day02/data/` listing. The new Scaffold_25 VCF and all downstream checkpoints must be prepared from the same verified four-population source; the earlier directory listing supplied only a different scaffold. The SnpEff database ID follows the supplied recipe, subject to the assembly verification above. Installation and database building belong in [software setup](software/README.md), not in the timed lessons.

## Begin

Before class, follow the [software installation and verification guide](software/README.md). It installs the core tools and the four commands used by the optional GERP lesson for every student. Then start with [the input and sample check](lessons/00-inputs.md).

