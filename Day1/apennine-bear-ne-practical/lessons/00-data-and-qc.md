# Terminal orientation and data QC


## Prepare your working directory

These commands create a personal working directory and connect it to the course files.

### 1 — Check your current location

```bash
pwd
```

`pwd` means **print working directory**. It displays the directory in which you are currently working.

### 2 — Create the Day 1 directory

```bash
mkdir -p AdvConGen/Day1
```

`mkdir` creates a directory. The `-p` option also creates any missing parent directories and prevents an error if they already exist.

### 3 — Enter the new directory

```bash
cd AdvConGen/Day1
```

`cd` means **change directory**. All subsequent commands will run from your personal Day 1 directory.

### 4 — Copy the input data

```bash
cp -r /home/user*/Share/Physalia_AdvancedConservationGenomics/Day1/data .
```

`cp` copies files, while `-r` copies the complete directory and its contents. The final `.` means “copy it into the current directory.”

Each student receives a personal copy of the data that can be modified without affecting the shared course files.

### 5 — Link the software directory

```bash
ln -s /home/ubuntu/Share/Physalia_AdvancedConservationGenomics/Day1/software .
```

`ln -s` creates a **symbolic link**, or shortcut, to the shared software directory. The software is not duplicated.

### 6 — Link the scripts directory

```bash
ln -s /home/ubuntu/Share/Physalia_AdvancedConservationGenomics/Day1/scripts .
```

This creates a shortcut to the shared course scripts. Any updates made by the instructor will therefore be immediately available to everyone.

### 7 — Create the results directory

```bash
mkdir -p results
```

This creates a personal directory where the outputs generated during the exercises will be stored.

### 8 — Check the directory structure

```bash
ls
```

`ls` lists the files and directories in the current location. You should see:

```text
data
results
scripts
software
```

You can obtain a more detailed listing with:

```bash
ls -lh
```

The `-l` option displays details, while `-h` reports file sizes in a human-readable format. Symbolic links are shown with an arrow pointing to their shared location.

**From the Day 1 course directory, activate the course environment:**

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


