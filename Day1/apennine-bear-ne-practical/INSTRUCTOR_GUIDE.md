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
results/gone/precomputed/Output_Ne_<name>
results/neestimator/precomputed/result.txt
```

Keep a provenance record containing source accession, sample identifiers, filtering commands, software versions, and checksums. If public redistribution is not allowed, provide a download/preparation script or use a simulated dataset.

## Facilitation notes

### MSMC2

Emphasize that the callable mask contributes information: long distances between heterozygous sites mean something only when the intervening sequence was observable. Ask students to distinguish a lack of variants from a lack of data.

Expected interpretation: the middle of the curve is generally more defensible than either extreme. One chromosome increases stochastic variation and reduces independent genealogical information.

### GONE2

Ask students why the full LD pattern—not an LD-pruned panel—is needed. Recent migrants, Wahlund effects, and relatives can create LD that is not due solely to historical population size.

GONE2 is the primary recent-demography practical. Students create the ABB and SBB population VCFs and the shared PED/MAP inputs, then run GONE2 with identical settings. Emphasize the sample size, SNP count, Fis, inferred genome length, and structure warnings in each **GONE2_STATS** file before discussing the trajectories.

The lesson uses GONE2's default upper recombination fraction of 0.05. The value 0.0101 is reserved for instructor sensitivity testing against the original GONE setting **hc=0.01**.

### Optional original GONE

The optional lesson reuses the PED/MAP files prepared by GONE2. The supplied driver assumes chromosomes are consecutively numbered beginning at 1, so Part 2 creates a derived one-chromosome MAP with chromosome code **1** while preserving marker order and physical coordinates.

Each student works in a private course directory, so no group-specific run copies are required. The original driver uses generic intermediate files and ABB and SBB runs must be performed sequentially. Keep compilation and installation outside the timed practical; all setup is consolidated in **software/README.md**.

### NeEstimator

Ask why GONE and NeEstimator should not receive identical SNP treatment automatically. Prepare a valid GENEPOP file and test it with the exact NeEstimator V2.1 package distributed to students.

## Suggested assessment

Award credit for reasoning, not whether a student's estimate matches an instructor value:

- 30% connects each method to its genomic signal;
- 25% distinguishes time scales;
- 25% identifies assumption violations;
- 20% communicates a cautious conservation interpretation.

## Ethical and conservation context

Avoid presenting a sensitive population merely as a convenient dataset. Discuss data sovereignty, sample scarcity, and why genetic results do not directly prescribe management. A low inferred `Ne` does not by itself determine whether genetic rescue, habitat intervention, or another action is appropriate.
