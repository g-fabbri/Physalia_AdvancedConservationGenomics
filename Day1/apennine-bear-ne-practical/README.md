# Two bear populations, three views of effective population size

## A conservation-genomics practical

The Apennine brown bear is a small and isolated population in central Italy. Slovak brown bears belong to a larger and more connected European population. Their contrasting histories give us a biological question:

> Do genomes from Apennine brown bears (ABB) and Slovak brown bears (SBB) record different histories of effective population size?

We will not search for one definitive value. Instead, we will use three genomic signals:

| Method | Genomic signal | Main time scale |
|---|---|---|
| MSMC2 | Coalescence along diploid genomes | Older and intermediate history |
| GONE | LD at different recombination distances | Recent generations |
| GONE2 (optional test) | Updated LD modelling and diagnostics | Recent generations |
| NeEstimator | LD in a population sample | Contemporary effective size |

The objective is to understand why these methods may produce different—but not necessarily contradictory—answers.

## The bears

Our chromosome-level dataset contains 18 individuals:

| Population | Identification | Sample size |
|---|---|---:|
| ABB — Apennine brown bears | IDs not beginning with U | 10 |
| SBB — Slovak brown bears | IDs beginning with U | 8 |

For the individual-genome comparison, we use ABB individual **4573** and SBB individual **U1916**. These abbreviations are used throughout the practical.

All analyses use **Scaffold_34** to keep computation short. A single scaffold is appropriate for learning the workflow, but provides less information and greater stochastic variation than a genome-wide analysis.

## The investigation

~~~text
Meet and inspect the samples
             ↓
Ask what one genome records
             ↓
Compare historical trajectories with MSMC2
             ↓
Ask what population-level LD records
             ↓
Compare recent histories with GONE
             ↓
Optionally repeat with GONE2
             ↓
Estimate contemporary Ne with NeEstimator
             ↓
Combine evidence and identify limitations
~~~

## Learning outcomes

By the end, you should be able to connect each estimator to its genomic signal and time scale, prepare and validate inputs, justify important parameters, and compare populations without ignoring sampling limitations.

## Schedule

| Time | Activity |
|---:|---|
| 00:00–00:15 | Effective population size and the bear case study |
| 00:15–00:25 | [Meet and inspect the data](lessons/00-data-and-qc.md) |
| 00:25–00:50 | Coalescence and historical demography |
| 00:50–01:20 | [Compare individual histories with MSMC2](lessons/01-msmc2.md) |
| 01:20–01:30 | Break |
| 01:30–01:50 | LD and recent demography |
| 01:50–02:20 | [Compare population histories with GONE](lessons/02-gone.md), or test [GONE2](lessons/02b-gone2.md) |
| 02:20–02:35 | Contemporary effective size and sampling |
| 02:35–02:55 | [Compare contemporary estimates](lessons/03-neestimator.md) |
| 02:55–03:00 | [Combine the evidence](lessons/04-synthesis.md) |

## Input files

~~~text
data/UrArMa_18i_s34.vcf.gz
data/UrArMa_18i_s34.vcf.gz.csi
data/UrArMa_4573_s34.vcf.gz
data/UrArMa_4573_s34.vcf.gz.csi
data/UrArMa_U1916_s34.vcf.gz
data/UrArMa_U1916_s34.vcf.gz.csi
data/UrArMa_callable.bed.gz
~~~

The 18-individual VCF supplies population data. The two already-filtered single-individual VCFs provide the ABB and SBB representatives for the historical comparison.

### About the callable mask

All individuals were aligned to the Apennine brown bear reference assembly. The course uses **UrArMa_callable.bed.gz** as a common teaching mask.

Sharing a reference assembly does not automatically make a sample-specific callable mask transferable. A shared mask is appropriate when it describes reference mappability or regions callable in all relevant samples. In a complete analysis, depth- and genotype-quality-based callability should be assessed separately for each individual and then combined explicitly.

## Before starting

~~~bash
bcftools --version
plink --version
msmc2 --help | head
Rscript --version
~~~

GONE, GONE2, and NeEstimator require separate software setup. Installation is not part of the timed practical. GONE2 is retained as an optional comparison while both LD workflows are being tested.

## Interpretation limits

- Ten Apennine and eight Slovak individuals are small samples for LD-based estimation.
- One scaffold provides much less independent information than a whole genome.
- Coverage, callability, missingness, relatedness, and population structure can imitate demographic differences.
- Classroom runs use fewer replicates to reduce runtime.
- Results demonstrate methods and hypotheses; they should not be used directly for management decisions.

## Begin

Start with [Terminal orientation and data QC](lessons/00-data-and-qc.md).

## Background and inspiration

- [Benazzo et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5692547/)
- [NCBI BioProject PRJNA395974](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA395974)
- Tutorial structure inspired by the narrative, command-first approach of the [Speciation & Population Genomics guide](https://speciationgenomics.github.io/pca/).

Additional method references are listed in [REFERENCES.md](REFERENCES.md).
