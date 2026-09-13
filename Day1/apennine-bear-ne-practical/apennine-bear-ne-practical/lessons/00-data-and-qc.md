# Terminal orientation and data QC

Estimated time: 10 minutes.

## Step 1 — Name the input files

**Purpose:** use short, consistent names and avoid repeatedly typing paths.

**Input:** the selected chromosome and the three teaching files.

~~~bash
CHROM=Scaffold_34
APN_SINGLE=data/UrArMa_4573_s34.vcf.gz
SVK_SINGLE=data/UrArMa_U1916_s34.vcf.gz
POP=data/UrArMa_18i_s34.vcf.gz
MASK=data/UrArMa_callable.bed.gz
~~~

**Expected:** nothing is printed because these commands assign shell variables.

**Check:** display their values:

~~~bash
printf 'Chromosome: %s\nApennine VCF: %s\nSlovak VCF: %s\nPopulation VCF: %s\nMask: %s\n' \
  "$CHROM" "$APN_SINGLE" "$SVK_SINGLE" "$POP" "$MASK"
~~~

Expected shape:

~~~text
Chromosome: Scaffold_34
Apennine VCF: data/UrArMa_4573_s34.vcf.gz
Slovak VCF: data/UrArMa_U1916_s34.vcf.gz
Population VCF: data/UrArMa_18i_s34.vcf.gz
Mask: data/UrArMa_callable.bed.gz
~~~

## Step 2 — Verify that the files exist

**Purpose:** confirm that all inputs are present and non-empty.

**Input:** the three paths defined above.

~~~bash
ls -lh "$APN_SINGLE" "$SVK_SINGLE" "$POP" "$MASK"
~~~

**Expected:** four lines containing file sizes greater than zero.

**Check:** “No such file or directory” means the working directory or input path is wrong. Run:

~~~bash
pwd
ls -lh data
~~~

## Step 3 — Inspect the samples

**Purpose:** identify the individuals represented in each VCF.

**Input:** the single-individual and population VCFs.

~~~bash
bcftools query -l "$APN_SINGLE"
bcftools query -l "$SVK_SINGLE"
bcftools query -l "$POP"
~~~

**Expected:** the first two commands print `4573` and `U1916`. The third prints one ID per population sample.

**Check:** count the samples:

~~~bash
bcftools query -l "$APN_SINGLE" | wc -l
bcftools query -l "$SVK_SINGLE" | wc -l
bcftools query -l "$POP" | wc -l
~~~

**Expected:** **1**, **1**, and **18**. Record all observed counts.

Create population sample lists from the naming convention:

~~~bash
bcftools query -l "$POP" | awk '/^U/' > data/slovak.samples
bcftools query -l "$POP" | awk '!/^U/' > data/apennine.samples
~~~

Check the population assignments:

~~~bash
printf 'Apennine bears: '; wc -l < data/apennine.samples
printf 'Slovak bears: '; wc -l < data/slovak.samples
cat data/apennine.samples
cat data/slovak.samples
~~~

**Expected:** 10 Apennine identifiers and 8 Slovak identifiers. `U1916` should occur in the Slovak list, while `4573` should occur in the Apennine list.

## Step 4 — Inspect chromosome labels

**Purpose:** confirm that the two VCFs and callable mask use the same chromosome identifier.

**Input:** both VCFs and the BED mask.

~~~bash
bcftools query -f '%CHROM\n' "$APN_SINGLE" | sort -u
bcftools query -f '%CHROM\n' "$SVK_SINGLE" | sort -u
bcftools query -f '%CHROM\n' "$POP" | sort -u
gzip -cd "$MASK" | cut -f1 | sort -u | head
~~~

**Expected:** both VCF commands include **Scaffold_34**. The mask command displays chromosome or scaffold labels present in the BED file.

**Check:** verify specifically that the mask contains the selected scaffold:

~~~bash
gzip -cd "$MASK" | awk -v chrom="$CHROM" '$1==chrom {found=1; exit} END {if (found) print chrom, "found"; else print chrom, "NOT FOUND"}'
~~~

Expected:

~~~text
Scaffold_34 found
~~~

Do not continue if the result is **NOT FOUND**.

## Step 5 — Count VCF records

**Purpose:** count the variant records in each VCF.

**Input:** the two indexed VCFs.

~~~bash
bcftools index -n "$APN_SINGLE"
bcftools index -n "$SVK_SINGLE"
bcftools index -n "$POP"
~~~

**Expected:** one integer from each command. These are counts of VCF records, not counts of callable base pairs.

**Check:** if an index is missing, create it and repeat:

~~~bash
bcftools index -f "$APN_SINGLE"
bcftools index -f "$SVK_SINGLE"
bcftools index -f "$POP"
bcftools index -n "$APN_SINGLE"
bcftools index -n "$SVK_SINGLE"
bcftools index -n "$POP"
~~~

Inspect several records:

~~~bash
bcftools view -H "$APN_SINGLE" | head
bcftools view -H "$SVK_SINGLE" | head
bcftools view -H "$POP" | head
~~~

## Step 6 — Inspect the callable mask

**Purpose:** inspect callable intervals from the selected scaffold only.

**Input:** the callable BED mask.

~~~bash
gzip -cd "$MASK" | \
  awk -v chrom="$CHROM" '$1==chrom' | head
~~~

**Expected:** records in BED format:

~~~text
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
~~~

Each row describes one callable interval:

~~~text
chromosome    start    end
~~~

**Check:** all displayed rows should begin with **Scaffold_34**, and each end coordinate should be larger than its start coordinate.

## Step 7 — Count BED intervals and callable base pairs

**Purpose:** distinguish the number of BED intervals from the amount of callable sequence.

**Input:** callable intervals belonging to **Scaffold_34**.

Count the intervals:

~~~bash
gzip -cd "$MASK" | \
  awk -v chrom="$CHROM" '$1==chrom' | wc -l
~~~

Sum their lengths:

~~~bash
gzip -cd "$MASK" | \
  awk -v chrom="$CHROM" '$1==chrom {bp += $3-$2} END {print "Callable bp:", bp}'
~~~

**Expected:**

- the first command prints the number of callable intervals on **Scaffold_34**;
- the second prints **Callable bp:** followed by a positive integer.

**Check:** BED intervals are zero-based and half-open. Their lengths are therefore calculated as column 3 minus column 2. The callable-base total should be positive and should not exceed the length of **Scaffold_34**.

## Why are VCF and BED files counted differently?

A VCF normally contains one row per variant record. Counting VCF rows tells us how many recorded variant sites are available.

A BED mask contains genomic intervals. Counting its rows gives the number of intervals, not the number of covered bases. Callable base pairs must be calculated by summing **end − start** across the intervals belonging to the selected scaffold.

## Stop and discuss

1. What is the main difference between the single-individual VCF and the population VCF?
2. How were individuals assigned to the two populations, and what assumption does that naming rule make?
3. Why must callable sequence in the BED mask be counted differently from variant records in the VCF?
4. Does the mask represent callable sequence or only variant positions?
5. What happens if the BED and VCF chromosome labels differ?

Record the sample counts, VCF record counts, BED interval count, and callable-base total in the [answer sheet](../answers/student_answers.md). Continue to the [first analysis](01-msmc2.md).
