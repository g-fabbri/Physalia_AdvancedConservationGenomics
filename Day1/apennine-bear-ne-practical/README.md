# Effective population size from genomes

## A terminal-based practical using the Apennine brown bear

This repository replaces a slide deck. Work through it from top to bottom: read a short concept, answer a prediction question, copy a command into the terminal, inspect its output, and interpret the result.

The central lesson is that **effective population size is not one universal number**. Different methods observe different genomic signals over different time windows.

> One diploid genome can be used for MSMC2. GONE and NeEstimator require genotypes from multiple individuals.

## Learning outcomes

By the end, you should be able to connect each method to its genomic signal and time scale, justify important parameters, and recognize when a conservation interpretation is unsupported.

## Three-hour route

| Time | Activity |
|---:|---|
| 00:00–00:15 | Introduction: what is `Ne`? |
| 00:15–00:25 | [Terminal orientation and data QC](lessons/00-data-and-qc.md) |
| 00:25–00:50 | MSMC2 concepts |
| 00:50–01:20 | [MSMC2 practical](lessons/01-msmc2.md) |
| 01:20–01:30 | Break |
| 01:30–01:50 | LD and recent demography |
| 01:50–02:20 | [GONE practical](lessons/02-gone.md) |
| 02:20–02:35 | Contemporary `Ne` and sampling |
| 02:35–02:55 | [NeEstimator practical](lessons/03-neestimator.md) |
| 02:55–03:00 | [Synthesis](lessons/04-synthesis.md) |

## How to use the commands

Copy commands inside code blocks. Do not copy a shell prompt such as `$`. Replace placeholders such as `chrN` with the value supplied by the instructor. Lines beginning with `#` are comments.

Start in the repository directory:

```bash
pwd
ls
```

Check the principal programs:

```bash
bcftools --version
plink --version
msmc2 --help | head
Rscript --version
```

If a command is missing, ask the instructor. Installation is not part of the timed practical.

## Dataset

```text
data/teaching/single_bear.chrN.vcf.gz
data/teaching/single_bear.chrN.callable.bed.gz
data/teaching/population.chrN.vcf.gz
```

The first two files represent one high-coverage diploid individual. The population VCF contains multiple individuals. See [data preparation and provenance](data/README.md).

## Case-study reading

- [Genomic history of the Apennine brown bear](https://pmc.ncbi.nlm.nih.gov/articles/PMC5692547/)
- [NCBI BioProject PRJNA395974](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA395974)

## Take-home message

> An `Ne` estimate belongs to a method, sample, genomic signal, time window, and set of assumptions—not simply to a species.

## Instructor resources

- [Instructor guide](INSTRUCTOR_GUIDE.md)
- [Student answer sheet](answers/student_answers.md)
- [Software environment](environment.yml)
- [References](REFERENCES.md)

