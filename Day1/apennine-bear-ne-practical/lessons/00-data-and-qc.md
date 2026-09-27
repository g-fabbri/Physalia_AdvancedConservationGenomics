# Terminal orientation and data QC


## Start your terminal

Prepare your working directory:

~~~bash
pwd
~~~

~~~bash
mkdir -p AdvConGen/Day1
cd AdvConGen/Day1
~~~

~~~bash
cp -r /home/ubuntu/Share/Physalia_AdvancedConservationGenomics/Day1/data .
~~~

~~~bash
ln -s /home/ubuntu/Share/Physalia_AdvancedConservationGenomics/Day1/software .
ln -s /home/ubuntu/Share/Physalia_AdvancedConservationGenomics/Day1/scripts .
~~~

~~~bash
mkdir results
~~~


From the Day 1 course directory, activate the course environment:

~~~bash
conda activate bear-ne-practical
~~~

Run this command whenever you open a new terminal. It makes the command-line programs used in this lesson, including BCFtools, available.

**Check:**

~~~bash
command -v bcftools
~~~

The command should print the path to `bcftools` inside the `bear-ne-practical` environment.


## Step 1 — Name the input files

**Purpose:** use short, consistent names and avoid repeatedly typing paths.

**Input:** the selected chromosome and the three teaching files.

~~~bash
CHROM=Scaffold_25
APN_SINGLE=data/UrArMa_4573_s25.vcf.gz
SVK_SINGLE=data/UrArMa_U1916_s25.vcf.gz
POP=data/UrArMa_18i_s25.vcf.gz
MASK=data/UrArMa_callable_s25.bed.gz
~~~

**Expected:** nothing is printed because these commands assign shell variables.

**Check:** display their values:

~~~bash
printf 'Chromosome: %s\nApennine VCF: %s\nSlovak VCF: %s\nPopulation VCF: %s\nMask: %s\n' \
  "$CHROM" "$APN_SINGLE" "$SVK_SINGLE" "$POP" "$MASK"
~~~

Expected shape:

~~~text
Chromosome: Scaffold_25
Apennine VCF: data/UrArMa_4573_s25.vcf.gz
Slovak VCF: data/UrArMa_U1916_s25.vcf.gz
Population VCF: data/UrArMa_18i_s25.vcf.gz
Mask: data/UrArMa_callable_s25.bed.gz
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
~~~

~~~bash
bcftools query -l "$SVK_SINGLE"
~~~

~~~bash
bcftools query -l "$POP"
~~~

**Expected:** the first two commands print `4573` and `U1916`. The third prints one ID per population sample.

**Check:** count the samples:

~~~bash
bcftools query -l "$APN_SINGLE" | wc -l
~~~

~~~bash
bcftools query -l "$SVK_SINGLE" | wc -l
~~~

~~~bash
bcftools query -l "$POP" | wc -l
~~~

**Expected:** **1**, **1**, and **18**. Record all observed counts.

Create population sample lists from the naming convention:

~~~bash
bcftools query -l "$POP" | awk '/^U/' > data/slovak.samples
~~~

~~~bash
bcftools query -l "$POP" | awk '!/^U/' > data/apennine.samples
~~~

Check the population assignments:

~~~bash
printf 'Apennine bears: '; wc -l < data/apennine.samples
printf 'Slovak bears: '; wc -l < data/slovak.samples
~~~

~~~bash
cat data/apennine.samples
~~~

~~~bash
cat data/slovak.samples
~~~

**Expected:** 10 Apennine identifiers and 8 Slovak identifiers. `U1916` should occur in the Slovak list, while `4573` should occur in the Apennine list.

## Step 4 — Inspect chromosome labels

**Purpose:** confirm that the two VCFs and callable mask use the same chromosome identifier.

**Input:** both VCFs and the BED mask.

~~~bash
bcftools query -f '%CHROM\n' "$APN_SINGLE" | sort -u
~~~

~~~bash
bcftools query -f '%CHROM\n' "$SVK_SINGLE" | sort -u
~~~

~~~bash
bcftools query -f '%CHROM\n' "$POP" | sort -u
~~~

~~~bash
gzip -cd "$MASK" | cut -f1 | sort -u | head
~~~

**Expected:** both VCF commands include **Scaffold_25**. The mask command displays chromosome or scaffold labels present in the BED file.

**Check:** verify specifically that the mask contains the selected scaffold:

~~~bash
gzip -cd "$MASK" | awk -v chrom="$CHROM" '$1==chrom {found=1; exit} END {if (found) print chrom, "found"; else print chrom, "NOT FOUND"}'
~~~

Expected:

~~~text
Scaffold_25 found
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
Scaffold_25     0       2
Scaffold_25     519     588
Scaffold_25     1130    1270
Scaffold_25     2005    2006
Scaffold_25     8267    8269
Scaffold_25     8509    8565
Scaffold_25     8611    8616
Scaffold_25     8885    9006
Scaffold_25     11667   11781
Scaffold_25     11793   11802
~~~

Each row describes one callable interval:

~~~text
chromosome    start    end
~~~

**Check:** all displayed rows should begin with **Scaffold_25**, and each end coordinate should be larger than its start coordinate.

## Step 7 — Count BED intervals and callable base pairs

**Purpose:** distinguish the number of BED intervals from the amount of callable sequence.

**Input:** callable intervals belonging to **Scaffold_25**.

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

- the first command prints the number of callable intervals on **Scaffold_25**;
- the second prints **Callable bp:** followed by a positive integer.

**Check:** BED intervals are zero-based and half-open. Their lengths are therefore calculated as column 3 minus column 2. The callable-base total should be positive and should not exceed the length of **Scaffold_25** (45,511,629 bp).

## Why are VCF and BED files counted differently?

A VCF normally contains one row per variant record. Counting VCF rows tells us how many recorded variant sites are available.

A BED mask contains genomic intervals. Counting its rows gives the number of intervals, not the number of covered bases. Callable base pairs must be calculated by summing **end − start** across the intervals belonging to the selected scaffold.


