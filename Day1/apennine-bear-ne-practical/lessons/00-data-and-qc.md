# Terminal orientation and data QC

Estimated time: 10 minutes.

## Predict

Which dataset should contain one individual, and which should contain several? Why?

Define the inputs, replacing `chrN` with the chromosome selected by the instructor:

```bash
CHROM=chrN
SINGLE=data/teaching/single_bear.${CHROM}.vcf.gz
POP=data/teaching/population.${CHROM}.vcf.gz
MASK=data/teaching/single_bear.${CHROM}.callable.bed.gz
```

Check the files and samples:

```bash
ls -lh "$SINGLE" "$POP" "$MASK"
bcftools query -l "$SINGLE"
bcftools query -l "$POP"
```

Count samples and variants:

```bash
bcftools query -l "$SINGLE" | wc -l
bcftools query -l "$POP" | wc -l
bcftools index -n "$SINGLE"
bcftools index -n "$POP"
```

Inspect chromosome labels:

```bash
bcftools query -f '%CHROM\n' "$POP" | sort -u
```

Summarize the callable mask:

```bash
zcat "$MASK" | awk '{bp += $3-$2} END {print "Callable bp:", bp}'
```

| Element | Meaning |
|---|---|
| `bcftools query -l` | Print sample names |
| `bcftools index -n` | Count indexed VCF records |
| `sort -u` | Return unique chromosome labels |
| `awk` expression | Sum BED interval lengths |

## Discuss

1. Is the MSMC2 VCF truly single-sample?
2. Is the population sample large enough for LD estimation?
3. Does the mask represent callable sequence or only variant positions?
4. What happens if BED and VCF chromosome names differ?

Continue to [MSMC2](01-msmc2.md).

