# Terminal orientation and data QC

Estimated time: 10 minutes.


## Step 1 — Name the input files

**Purpose:** use short, consistent names and avoid repeatedly typing paths.

**Input:** the chromosome identifier selected.

```bash
CHROM=Scaffold_34
SINGLE=data/UrArMa_4573_s34.vcf.gz
POP=data/UrArMa_10i_s34.vcf.gz
MASK=data/UrArMa_callable.bed.gz
```

**Expected:** nothing is printed because these commands assign shell variables.

**Check:** display their values:

```bash
printf 'Chromosome: %s\nSingle VCF: %s\nPopulation VCF: %s\nMask: %s\n' \
  "$CHROM" "$SINGLE" "$POP" "$MASK"
```

Expected shape:

```text
Chromosome: Scaffold_34
Single VCF: data/UrArMa_4573_s34.vcf.gz
Population VCF: data/UrArMa_10i_s34.vcf.gz
Mask: data/UrArMa_callable.bed.gz
```

## Step 2 — Verify that the files exist

**Input:** the three paths above.

```bash
ls -lh "$SINGLE" "$POP" "$MASK"
```

**Expected:** three lines containing file sizes greater than zero.

**Check:** `No such file or directory` means the chromosome label, working directory, or input path is wrong. Run `pwd` and `ls data/teaching` before asking for help.

## Step 3 — Inspect the samples

```bash
bcftools query -l "$SINGLE"
bcftools query -l "$POP"
```

**Expected:** the first command prints one sample ID. The second prints one ID per population sample.

Count them:

```bash
bcftools query -l "$SINGLE" | wc -l
bcftools query -l "$POP" | wc -l
```

**Expected:** `1` for the single-individual VCF and a value greater than `1` for the population VCF. Record the actual population count.

## Step 4 — Inspect chromosomes and record counts

```bash
bcftools query -f '%CHROM\n' "$POP" | sort -u
```

**Expected:** exactly one line matching `$CHROM`.

```bash
bcftools index -n "$SINGLE"
bcftools index -n "$POP"
```

**Expected:** one integer from each command. This is the number of VCF records, not the number of callable bases.


## Step 5 — Inspect callable sequence

```bash
zcat "$MASK" | head
```

Expected format—your coordinates will differ:

```text
Scaffold_34     0       740
Scaffold_34     769     1567
Scaffold_34     1914    2013
Scaffold_34     2102    3531
Scaffold_34     3575    4060
Scaffold_34     4171    6326
Scaffold_34     6403    6506
Scaffold_34     6890    7221
Scaffold_34     7874    8952
Scaffold_34     8989    9807
```

Sum interval lengths:

```bash
zcat "$MASK" | awk '{bp += $3-$2} END {print "Callable bp:", bp}'
```

**Expected:** `Callable bp:` followed by a positive integer. BED intervals are zero-based and half-open, so interval length is column 3 minus column 2.

## Stop and discuss

1. Is the MSMC2 VCF truly single-sample?
2. Is the population sample likely to contain enough individuals for reliable LD estimation?
3. Does the mask represent callable sequence or only variant positions?
4. What happens if BED and VCF chromosome labels differ?

Record sample counts, variant counts, and callable bases in the [answer sheet](../answers/student_answers.md). Continue to [MSMC2](01-msmc2.md).

