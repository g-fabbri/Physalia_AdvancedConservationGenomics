# Data preparation and provenance

No empirical genotypes are bundled in this template. The instructor should stage approved derived files under `data/`.

## Candidate source

The Apennine brown bear data described by Benazzo et al. (2017) are linked to [NCBI BioProject PRJNA395974](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA395974). The paper reports six Apennine samples, with one excluded for very low coverage, and heterogeneous sequencing depth among the retained genomes.

Raw reads are not appropriate to download or process during class. Generate the teaching VCF and callable masks beforehand, or obtain author-approved derived data.

## Minimum QC record

Document:

- reference assembly and contig name;
- sample IDs and population labels;
- sequencing depth and genotype missingness;
- variant-calling workflow;
- SNP, genotype-quality, depth, and missingness filters;
- treatment of repeats and low-mappability regions;
- relatedness filtering;
- callable-mask construction;
- mutation rate, generation time, and recombination assumptions;
- software versions and command history.

## Files expected by the exercises

```text
data/UrArMa_4573_s34.vcf.gz
data/UrArMa_4573_s34.vcf.gz.tbi
data/UrArMa_callable.bed.gz
data/UrArMa_18i_s34.vcf.gz
data/UrArMa_18i_s34.vcf.gz.tbi
```

The VCFs should contain both invariant callable sites or have an accompanying mask appropriate to the downstream method. A variant-only VCF is insufficient for constructing a defensible MSMC2 callable mask.

For the comparative tutorial, sample `4573` represents the Apennine population and `U1916` represents the Slovak population. `UrArMa_callable.bed.gz` is treated as a common teaching mask. Document whether it is a reference/mappability mask, a cohort-wide mask, or an intersection of individual callable masks. Alignment to the same reference alone does not justify reusing a sample-specific depth mask.
