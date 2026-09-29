# Part 0 — Meet the four-population dataset

In Day 2, we will work with a VCF containing four bear populations or species:

- **ABB:** Apennine brown bears — focal population;
- **SBB:** Slovak brown bears — focal population;
- **BLB:** black bears — outgroup;
- **POB:** polar bears — outgroup.

ABB and SBB are the populations whose genetic variation we want to compare. Black and polar bears will be used as outgroups to help determine which alleles are ancestral and which are derived.

All samples were called against the Apennine-bear reference genome. Therefore, the VCF `REF` allele is the allele in the reference assembly, but it is not necessarily the ancestral allele.

To keep the practical manageable, we will analyse only **Scaffold 25**.


<img width="1672" height="941" alt="bears" src="https://github.com/user-attachments/assets/fe759e80-2d70-4a5c-9bfe-15cde52b9686" />



## Prepare your working directory

These commands create a personal working directory and connect it to the course files.

**Check your current location**

```bash
pwd
```

**Create the Day 2 directory**

```bash
mkdir -p AdvConGen/Day2
```

**Enter the new directory**

```bash
cd AdvConGen/Day2
```

**Copy the input data**

**Replace USER with your username**

```bash
cp -r /home/USER/Share/Physalia_AdvancedConservationGenomics/Day2/data .
```

Each student receives a personal copy of the data that can be modified without affecting the shared course files.

**Link the software/scripts directories**

**Replace USER with your username**

```bash
ln -s /home/USER/Share/Physalia_AdvancedConservationGenomics/Day2/software .
ln -s /home/USER/Share/Physalia_AdvancedConservationGenomics/Day2/scripts .
ln -s /home/USER/Share/Physalia_AdvancedConservationGenomics/Day2/gerp_input .
```

This creates a shortcut to the shared course scripts. Any updates made by the instructor will therefore be immediately available to everyone.

**Create the results directory**

```bash
mkdir -p results
```

This creates a personal directory where the outputs generated during the exercises will be stored.

**Check the directory structure**

```bash
ls
```


**From your Day 2 directory, run:**

```bash
conda activate bear-load-practical

COURSE_DIR=$(pwd)
export PATH="$COURSE_DIR/software/bin:$PATH"
```
## Step 1 — Define the input files

The prepared VCF and sample lists are stored in the `data/` directory:

```bash
VCF="$COURSE_DIR/data/Bears_4pops_s25.vcf.gz"

ABB_LIST="$COURSE_DIR/data/ABB.samples"
SBB_LIST="$COURSE_DIR/data/SBB.samples"
BLB_LIST="$COURSE_DIR/data/BLB.samples"
POB_LIST="$COURSE_DIR/data/POB.samples"
```

The VCF contains the genotypes, while each sample list contains one sample ID per line.

Check that the files exist:

```bash
ls -lh \
  "$VCF" \
  "$VCF.csi" \
  "$ABB_LIST" \
  "$SBB_LIST" \
  "$BLB_LIST" \
  "$POB_LIST"
```

## Step 2 — Check the samples

Count the total number of samples in the VCF:

```bash
bcftools query -l "$VCF" | wc -l
```

Count the samples assigned to each group:

```bash
for LIST in \
  "$ABB_LIST" \
  "$SBB_LIST" \
  "$BLB_LIST" \
  "$POB_LIST"; do

  printf '%s: ' "$LIST"
  wc -l < "$LIST"
done
```

The counts from the four sample lists should add up to the total number of samples in the VCF.

## Step 3 — Check the variants

Inspect the first three variant records:

```bash
bcftools view -H "$VCF" | head -n 3
```

Count the variants:

```bash
bcftools index -n "$VCF"
```

Check which scaffold is present:

```bash
bcftools query -f '%CHROM\n' "$VCF" | sort -u
```

The final command should print only:

```text
Scaffold_25
```

The scaffold name must match the reference FASTA and the SnpEff database used in the next lesson.

Because this practical uses only one scaffold, the results are illustrative and should not automatically be interpreted as genome-wide estimates of genetic load.

