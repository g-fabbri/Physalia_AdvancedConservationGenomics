# Day 2 — Genetic load in two bear populations

## The question

Day 1 asked whether Apennine (ABB) and Slovak (SBB) brown bears have different demographic histories. Today we ask whether they carry different **burdens of putatively damaging derived alleles**. Black bears (BLB) and polar bears (PBB) provide outgroup information for deciding which allele is likely ancestral. Functional effects come from SnpEff; GenoLoader polarizes genotypes and records derived-allele dosage. A GERP conservation-score comparison is an optional extension.

An allele annotated `HIGH` or `MODERATE` is **not** a measured fitness effect. A larger derived-allele count is not automatically a larger realized load. Keep prediction, ancestry, genotype state, and fitness distinct throughout.

## Schedule (3 hours; about 90 minutes at the terminal)

| Time | Theory and practical, alternated |
|---|---|
| 00:00–00:25 | Genetic load, recessive versus additive effects, and why bottlenecks matter |
| 00:25–00:35 | [Orient to the four-species VCF](lessons/00-inputs.md) |
| 00:35–01:00 | Variant consequences and annotation uncertainty |
| 01:00–01:20 | [Inspect SnpEff annotations](lessons/01-snpeff.md) |
| 01:20–01:30 | Break |
| 01:30–01:55 | Ancestral-state polarization, outgroup disagreement, and introgression |
| 01:55–02:30 | [Polarize and compare burdens with GenoLoader](lessons/02-genoloader.md) |
| 02:30–02:45 | Conservation scores versus coding consequences |
| 02:45–03:00 | [GERP extension and synthesis](lessons/03-gerp.md) |

## Instructor preparation and missing local details

This is a **draft with explicit input placeholders**, not yet a tested runnable lesson. Before class, supply:

- a **biallelic SNP VCF** on the Apennine brown-bear reference, containing ABB, SBB, BLB, and PBB individuals, plus its index;
- four exact sample-ID lists in `data/` and the sample counts for each;
- the custom **UrArMar_mUrsArc2** SnpEff database, after resolving the mUrsArc1.1 FASTA-name discrepancy and verifying that its reference and gene model match the VCF;
- a small annotated teaching VCF, ideally on one well-covered autosome, so the exercise fits the time;
- GenoLoader compiled and tested against that VCF;
- for the optional GERP extension: the score track's **source assembly**, score definition and format, a validated source-to-Apennine chain, and a pre-lifted one-base score track with mapping QC.

The commands deliberately use `data/FOUR_SPECIES_TEACHING.vcf.gz` until the real classroom file name is confirmed. The SnpEff database ID now follows the supplied Jarvis recipe, subject to the assembly verification above. Installation and database building belong in [software setup](software/README.md), not in the timed lessons.

## Begin

Start with [the input and sample check](lessons/00-inputs.md). Record outputs in the [answer sheet](answers/student_answers.md).
