# Two bear populations, two complementary views of effective population size

## A conservation-genomics practical

The Apennine brown bear is a small and isolated population in central Italy. Slovak brown bears belong to a larger and more connected European population. Their contrasting histories give us a biological question:

> Do genomes from Apennine brown bears (ABB) and Slovak brown bears (SBB) record different histories of effective population size?

We will not search for one definitive value. Instead, we will use two complementary approaches, with optional method comparisons:

| Method | Genomic signal | Main time scale |
|---|---|---|
| MSMC2 | Coalescence along diploid genomes | Older and intermediate history |
| GONE2 | Updated LD modelling and diagnostics | Recent generations |
| SMC++ (optional extension) | Coalescent and frequency information from multiple unphased individuals | Historical trajectory |
| currentNe2 (optional extension) | Genome-wide LD in a population sample | Contemporary effective size |

The objective is to understand why these methods may produce different—but not necessarily contradictory—answers.

## The bears

Our chromosome-level dataset contains 18 individuals:

| Population | Identification | Sample size |
|---|---|---:|
| ABB — Apennine brown bears | IDs not beginning with U | 10 |
| SBB — Slovak brown bears | IDs beginning with U | 8 |

For the individual-genome comparison, we use ABB individual **4573** and SBB individual **U1916**. These abbreviations are used throughout the practical.

All analyses use **Scaffold_25** to keep computation short. A single scaffold is appropriate for learning the workflow, but provides less information and greater stochastic variation than a genome-wide analysis.

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
Compare recent histories with GONE2
             ↓
Optional extensions: SMC++ or currentNe2
             ↓
Combine evidence and identify limitations
~~~

## Learning outcomes

By the end, you should be able to connect each estimator to its genomic signal and time scale, prepare and validate inputs, justify important parameters, and compare populations without ignoring sampling limitations.


## Input files

~~~text
data/UrArMa_18i_s25.vcf.gz
data/UrArMa_18i_s25.vcf.gz.csi
data/UrArMa_4573_s25.vcf.gz
data/UrArMa_4573_s25.vcf.gz.csi
data/UrArMa_U1916_s25.vcf.gz
data/UrArMa_U1916_s25.vcf.gz.csi
data/UrArMa_callable_s25.bed.gz
~~~

The 18-individual VCF supplies population data. The two already-filtered single-individual VCFs provide the ABB and SBB representatives for the historical comparison.

### About the callable mask

All individuals were aligned to the Apennine brown bear reference assembly. The course uses **UrArMa_callable_s25.bed.gz** as a common teaching mask.

Sharing a reference assembly does not automatically make a sample-specific callable mask transferable. A shared mask is appropriate when it describes reference mappability or regions callable in all relevant samples. In a complete analysis, depth- and genotype-quality-based callability should be assessed separately for each individual and then combined explicitly.

## Before starting

From the course root, activate the prepared Conda environment and add **msmc-tools** to your current terminal session:

~~~bash
conda activate bear-ne-practical
export PATH="$PWD/software/msmc-tools:$PATH"
~~~

Repeat these two commands whenever you open a new terminal. Then check the programs:

~~~bash
bcftools --version
plink --version
msmc2 --help | head
generate_multihetsep.py --help | head
Rscript --version
~~~

Installation is not part of the timed practical. Create the shared environment and install the remaining programs using the consolidated [software setup](software/README.md).

## Interpretation limits

- Ten Apennine and eight Slovak individuals are small samples for LD-based estimation.
- One scaffold provides much less independent information than a whole genome.
- Coverage, callability, missingness, relatedness, and population structure can imitate demographic differences.
- Classroom runs use fewer replicates to reduce runtime.
  
## Begin

Start with [Terminal orientation and data QC](lessons/00-data-and-qc.md).

## Background

- [Benazzo et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5692547/)
- [NCBI BioProject PRJNA395974](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA395974)

Additional method references are listed in [REFERENCES.md](REFERENCES.md).
