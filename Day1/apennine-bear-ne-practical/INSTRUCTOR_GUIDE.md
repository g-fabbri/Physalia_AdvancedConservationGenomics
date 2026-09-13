# Instructor guide

## Recommended teaching configuration

Use Scaffold_34 for speed and distribute precomputed genome-wide outputs for comparison. Frame the practical around Apennine versus Slovak bears, while calling the one-scaffold analyses demonstrations rather than biological reconstructions.

### Required inputs

- the 18-individual joint VCF containing 10 Apennine and 8 Slovak bears;
- representative diploid genomes 4573 (Apennine) and U1916 (Slovak);
- a common mask valid for both representatives, plus sample-specific exclusion of missing genotypes;
- chromosome names consistent across VCF, BED, reference, and genetic map;
- no sex chromosome unless sex and ploidy are handled explicitly.

Do not manufacture a callable mask from variant-only positions. Verify how `UrArMa_callable.bed.gz` was constructed. Sharing an alignment reference permits a shared mappability mask but does not by itself justify sharing a sample-specific depth mask.

## Before class

1. Check permission to redistribute derived genotype data.
2. Choose an autosome with adequate callable sequence and SNP density.
3. Remove or flag close relatives and samples with excessive missingness.
4. Prepare all inputs and copy/paste every command from the lesson pages on the classroom platform.
5. Save intermediate and final checkpoints.
6. Test with the exact student account and available compute quota.
7. Decide and document mutation rate, generation time, and recombination assumptions.

## Checkpoints to distribute

```text
results/msmc2/precomputed/msmc2.final.txt
results/gone2/precomputed/GONE2_ABB_SBB.pdf
results/smcpp/precomputed/SMCPP_ABB_SBB_Scaffold_34.pdf
results/currentne2/precomputed/comparison.txt
```

Keep a provenance record containing source accession, sample identifiers, filtering commands, software versions, and checksums. If public redistribution is not allowed, provide a download/preparation script or use a simulated dataset.

## Facilitation notes

### MSMC2

Emphasize that the callable mask contributes information: long distances between heterozygous sites mean something only when the intervening sequence was observable. Ask students to distinguish a lack of variants from a lack of data.

Expected interpretation: the middle of the curve is generally more defensible than either extreme. One chromosome increases stochastic variation and reduces independent genealogical information.

### GONE2

Ask students why the full LD pattern—not an LD-pruned panel—is needed. Recent migrants, Wahlund effects, and relatives can create LD that is not due solely to historical population size.

GONE2 is the primary recent-demography practical. Students create the ABB and SBB population VCFs and the shared PED/MAP inputs, then run GONE2 with identical settings. Emphasize the sample size, SNP count, Fis, inferred genome length, and structure warnings in each **GONE2_STATS** file before discussing the trajectories.

The lesson uses GONE2's default upper recombination fraction of 0.05. Testing alternative `-u` values can be an instructor sensitivity analysis.

### Optional SMC++

This is an extension, not part of the timed practical. It uses the same one-scaffold 18-individual VCF but fits ABB and SBB separately. Verify that the VCF contains all expected sample IDs, that the reference chromosome length is available, and that the shared callable BED is appropriate for both populations. SMC++ needs an **uncallable** BED mask, so the lesson complements the provided callable BED. Test both fits and the plot on the classroom platform; distribute a precomputed plot if fitting takes too long.

### Optional currentNe2

This extension reuses the GONE2 PED/MAP inputs to estimate one contemporary value per population. It is LD-based like GONE2, but targets a different time scale. Treat its one-scaffold numbers as exploratory; prepare a multi-chromosome dataset for a substantive comparison, especially if using currentNe2's between-chromosome estimate.

## Suggested assessment

Award credit for reasoning, not whether a student's estimate matches an instructor value:

- 30% connects each method to its genomic signal;
- 25% distinguishes time scales;
- 25% identifies assumption violations;
- 20% communicates a cautious conservation interpretation.

## Ethical and conservation context

Avoid presenting a sensitive population merely as a convenient dataset. Discuss data sovereignty, sample scarcity, and why genetic results do not directly prescribe management. A low inferred `Ne` does not by itself determine whether genetic rescue, habitat intervention, or another action is appropriate.
