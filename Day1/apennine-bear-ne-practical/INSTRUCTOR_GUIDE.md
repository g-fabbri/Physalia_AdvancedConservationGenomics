# Instructor guide

## Recommended teaching configuration

Use one long autosome for speed and distribute precomputed outputs from a genome-wide analysis for comparison. Call the chromosome analysis a demonstration, not a biological reconstruction.

### Required inputs

- one high-coverage diploid Apennine bear with a sample-specific callable mask;
- a jointly called, filtered VCF from multiple unrelated Apennine bears;
- chromosome names consistent across VCF, BED, reference, and genetic map;
- no sex chromosome unless sex and ploidy are handled explicitly.

Do not manufacture a callable mask from variant-only positions. The mask must represent sites where a genotype could reliably have been called.

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

### GONE

Ask students why the full LD pattern—not an LD-pruned panel—is needed. Recent migrants, Wahlund effects, and relatives can create LD that is not due solely to historical population size.

The GONE workflow is platform-sensitive. Confirm the input basename expected by `script_GONE.sh` and provide the exact tested directory. Never make compilation the objective of a 30-minute exercise.

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
