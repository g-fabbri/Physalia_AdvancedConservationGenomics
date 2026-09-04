# NeEstimator — contemporary effective population size

Estimated practical time: 20 minutes.

The single-sample LD estimator measures non-random allele associations in a population sample and corrects for LD expected from finite sampling. It requires multilocus diploid genotypes from multiple individuals.

## Predict

1. What happens to uncertainty as the number of individuals decreases?
2. Why can rare alleles strongly affect LD estimates?
3. Why should physically linked markers not automatically be treated as independent?
4. Why does GONE need LD by distance while this analysis may use unlinked loci?

## Create a transparent teaching panel

### Step 1 — Define the files

```bash
CHROM=chrN
VCF=data/teaching/population.${CHROM}.vcf.gz
OUTDIR=results/neestimator
mkdir -p "$OUTDIR"
```

**Expected:** nothing is printed. Confirm the input with `ls -lh "$VCF"`.

### Step 2 — Apply variant and frequency filters

**Input:** the population VCF. **Output:** a compressed VCF containing biallelic SNPs above the selected frequency threshold.

```bash
bcftools view \
  -m2 -M2 \
  -v snps \
  -q 0.02:minor \
  "$VCF" \
  -Oz -o "$OUTDIR/population.maf02.vcf.gz"
```

**Expected:** no table is printed because output is written to the named file. Confirm it exists:

```bash
ls -lh "$OUTDIR/population.maf02.vcf.gz"
```

```bash
bcftools index -f "$OUTDIR/population.maf02.vcf.gz"
```

**Expected:** an index file ending `.csi` or `.tbi`. Check with:

```bash
ls "$OUTDIR"/population.maf02.vcf.gz.*
```

### Step 3 — Reduce physical linkage

For a simple classroom comparison, keep at most one SNP per 100 kb:

```bash
bcftools +prune \
  "$OUTDIR/population.maf02.vcf.gz" \
  -Oz -o "$OUTDIR/population.thinned.vcf.gz" \
  -- -n 1 -w 100kb
```

**Expected:** a new compressed VCF. The section after `--` contains options for the `prune` plugin rather than global bcftools options.

| Flag | Effect | Question |
|---|---|---|
| `-m2 -M2` | Exactly two observed alleles | Is this consistent with the estimator? |
| `-v snps` | Retains SNPs | Are indels wanted? |
| `-q 0.02:minor` | Minimum minor-allele frequency | How many copies is this in the sample? |
| `-Oz` | Writes block-gzipped VCF | Why is indexing useful? |
| `+prune -n 1 -w 100kb` | One variant per 100-kb window | Is physical spacing a good proxy for independence? |

This thinning rule is deliberately visible, not universally correct. Better choices can use different chromosomes, a recombination map, or empirical LD decay.

```bash
bcftools index -f "$OUTDIR/population.thinned.vcf.gz"
bcftools index -n "$OUTDIR/population.thinned.vcf.gz"
```

**Expected:** an integer smaller than or equal to the count in the MAF-filtered file. Compare explicitly:

```bash
bcftools index -n "$OUTDIR/population.maf02.vcf.gz"
bcftools index -n "$OUTDIR/population.thinned.vcf.gz"
```

**Check:** if thinning retains almost everything, the window may be too small or SNP density low. If it retains very few loci, precision may be poor.

## Check the converted input

### Step 4 — Validate GENEPOP formatting

The instructor supplies a tested VCF-to-GENEPOP conversion because file-format debugging is not the biological objective.

```bash
head -30 data/teaching/population.thinned.genepop
```

**Expected:** a title, locus names, a `Pop` separator, and one genotype row per individual. The exact allele width depends on the converter.

Check population assignments, unique locus names, missing-data encoding, and allele encoding.

## Run NeEstimator V2.1

### Step 5 — Launch and compare parameter choices

The official package normally uses a Java interface that calls the platform-specific executable:

```bash
java -jar NeEstimator2x1.jar
```

**Expected:** the graphical interface opens. A terminal message about a missing Java runtime or an unreadable JAR means the installation—not the genetic analysis—failed.

In the interface:

1. select the GENEPOP teaching file;
2. choose the single-sample LD method;
3. run with `Pcrit = 0.02`;
4. run again with `Pcrit = 0.05`;
5. retain the jackknife confidence intervals.

Save the two outputs with distinct names such as `LD_Pcrit002.txt` and `LD_Pcrit005.txt`. Confirm that each reports population/sample name, loci used, point estimate, and confidence limits.

The package contains batch-processing support, but control files must be prepared for the exact release used in class. Do not invent terminal syntax; inspect the bundled help and example batch files.

## Interpret

| `Pcrit` | Individuals | Loci | `Ne` | Lower CI | Upper CI |
|---:|---:|---:|---:|---:|---:|
| 0.02 | | | | | |
| 0.05 | | | | | |

Did the estimate change more or less than its uncertainty? Is an infinite upper bound informative? Does precision guarantee that structure, relatives, and linkage were handled correctly?

Continue to [synthesis](04-synthesis.md).
