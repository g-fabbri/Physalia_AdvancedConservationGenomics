# GONE — recent demographic history from LD

Estimated practical time: 30 minutes.

Drift creates associations between alleles and recombination breaks them down. LD at different recombination distances therefore contains information about different recent generations. GONE requires **multiple individuals**.

## Predict

Which non-demographic processes can generate LD? Consider migrants, pooled subpopulations, relatives, selection, and genotype error.

## Prepare PED/MAP input

```bash
CHROM=chrN
VCF=data/teaching/population.${CHROM}.vcf.gz
OUTDIR=results/gone
PREFIX="$OUTDIR/population_${CHROM}"
mkdir -p "$OUTDIR"
```

Count samples and biallelic SNPs:

```bash
bcftools query -l "$VCF" | wc -l
bcftools view -m2 -M2 -v snps "$VCF" -Ou | bcftools view -H | wc -l
```

```bash
plink \
  --vcf "$VCF" \
  --double-id \
  --allow-extra-chr \
  --snps-only just-acgt \
  --biallelic-only strict \
  --geno 0.10 \
  --mac 2 \
  --recode \
  --out "$PREFIX"
```

| Flag | Effect | Decision question |
|---|---|---|
| `--double-id` | Uses sample ID as family and individual ID | Are real relationships known? |
| `--allow-extra-chr` | Accepts non-human chromosome labels | Are small contigs mixed with autosomes? |
| `--snps-only just-acgt` | Keeps canonical SNP alleles | Why exclude indels here? |
| `--biallelic-only strict` | Keeps exactly two alleles | Is this required downstream? |
| `--geno 0.10` | Removes loci with >10% missing genotypes | Is 10% sensible for this sample size? |
| `--mac 2` | Requires two minor-allele copies | Would three copies be safer but discard too much? |
| `--recode` | Produces PED/MAP files | Does your GONE version expect these? |

### MAF or MAC?

For 10 diploid individuals, `--maf 0.02` corresponds to less than one chromosome copy and therefore does almost nothing. Re-run with `--mac 3` and compare logs:

```bash
plink \
  --vcf "$VCF" --double-id --allow-extra-chr \
  --snps-only just-acgt --biallelic-only strict \
  --geno 0.10 --mac 3 --recode \
  --out "${PREFIX}_mac3"
```

```bash
grep -E 'variants.*remaining|variants loaded' "$OUTDIR"/*.log
```

Which threshold would you choose, and why?

## Inspect the genetic map

```bash
head "${PREFIX}.map"
```

PLINK MAP columns are chromosome, marker ID, genetic position in centimorgans, and physical position. If column 3 is zero, legacy GONE uses the average recombination-rate setting `cMMb`.

Discuss:

1. Is a species-specific recombination map available?
2. What evidence supports an assumed cM/Mb value?
3. Is a constant rate realistic across the chromosome?
4. How could an incorrect rate distort the inferred time axis?

## Run the official GONE workflow

The instructor supplies a tested copy of the appropriate folder from the [official GONE repository](https://github.com/esrud/GONE).

```bash
cp "${PREFIX}.ped" software/GONE/
cp "${PREFIX}.map" software/GONE/
cd software/GONE
less INPUT_PARAMETERS_FILE
```

Before running, identify:

- number of replicates;
- maximum recombination distance;
- `cMMb` if no genetic map is supplied;
- maximum SNPs sampled per chromosome;
- number of generations reported.

For each, state the classroom value, a possible final-analysis value, and whether changing it affects runtime, resolution, bias, or several of these.

Run using the basename without `.ped` or `.map`:

```bash
bash script_GONE.sh population_chrN
```

The instructor must test the basename and exact behavior with the downloaded version. For the 30-minute exercise, use fewer replicates and compare with a precomputed full run.

## Diagnose before interpreting

- Were close relatives removed?
- Does `Fis` suggest departure from panmixia?
- Could structure or immigration mimic a sharp recent decline?
- How variable are replicate trajectories?
- How much confidence is possible from one chromosome?

Continue to [NeEstimator](03-neestimator.md).

