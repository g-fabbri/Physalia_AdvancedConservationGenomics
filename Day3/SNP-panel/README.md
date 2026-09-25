# Day 3 — SNP Panel Creation and Geographic Assignment

Welcome to the **Day 3 tutorial**!

Today we are embarking on the creation of a SNP panel that will help us manage our study species. In our example, we will be analyzing samples of Hermann’s tortoise (*Testudo hermanni*).

First, we need to describe the population structure of our study species. Then, we want to identify the SNPs that best characterize each subpopulation. This allows us to assign, with associated probability, the geographic origin of a new sample.

Indeed, in the last part of the tutorial, we will apply the newly generated SNP panel to a new sample of confiscated Hermann’s tortoise in order to find out if any of them can be considered suitable for possible future release in the wild in a specific geographic area.

## Objectives of the tutorial

By the end of this tutorial, you will:

* Get familiar with sample data and metadata
* Study population structure and identify possible substructure patterns
* Create a panel of informative SNPs to discern among different geographic areas
* Apply the SNP panel to new samples

## Datasets

We are going to work with two datasets:

1. **Wild population samples** — samples taken from wild populations across the species’ range to obtain a list of informative SNPs and characterize the different geographic areas.
2. **New samples of unknown origin** — samples to which we will apply our new SNP panel.

The first dataset was obtained by applying a ***de novo* ddRAD sequencing protocol**. The second dataset was produced at a subsequent stage by applying a genotyping technique that reads the genotype directly at the positions identified as informative by the SNP panel.

More on that coming up after the break!

---

# 1. Learn about your dataset

You should try to find out as much metadata as possible, especially about sampling locations, as these are crucial for the creation of the SNP panel.

Now, let’s get familiar with handling the dataset.

The file is in **STRUCTURE (`.stru`) format**. It is a tab-separated text file with marker names in the first row, followed by the genotype data distributed on one or two rows per sample.

Sometimes samples are separated into two lines. This is an important detail to notice when reading the file into R.

From the second row, you will see something like this:

```text
Ind1  PopA   0 1  0 0  1 1
Ind2  PopA   0 0  0 1  1 1
Ind3  PopB   1 1  0 1  0 0
Ind4  PopB   1 0  1 1  0 1
```

Here, `PopA`, `PopB`, etc. are also progressive numbers (e.g. `PopA = 1`, `PopB = 2`, `PopC = 3`), so be careful not to confuse them with alleles.

In order to correctly read the dataset into R, we need to know:

* How many samples are present?
* How many loci are present?
* How is missing data coded? (e.g. `NA`, `0`, `-9`)

So, let’s answer these questions:

* **How many samples?**
* **How many loci?**

Take note of these values in the `metadata.xlsx` Excel file that you can download from the GitHub folder:

```text
Day3/tutorial_snpchip/metadata.xlsx
```

> 💡 **Tip:** If you struggle with the coding, you can find the solution in:
>
> `Day3/tutorial_snpchip/solutions`

---

# 2. Study population structure of wild populations

## 2a. PCA

In this section, we will use **PCA** as the first exploratory analysis, as we will leverage it later for the SNP chip creation.

### R code

> 📄 **R script:** `script_pca.R`
>
> Add/link the R script here once your GitHub folder structure is finalized.

---

## 2b. sNMF

Now we will use another approach to delve into the population structure characterizing our dataset.

This is a much faster version of **STRUCTURE/ADMIXTURE**, but the output is similar. We need to understand the best number of ancestral components based on the **cross-entropy value** and then visualize how our samples are described based on that model.

### R code

> 📄 **R script:** `script_snmf.R`
>
> Add/link the R script here once your GitHub folder structure is finalized.

---

## 2c. F<sub>ST</sub>

We end this first part by computing **pairwise F<sub>ST</sub> values** for predefined groups.

This will allow us to evaluate, with numerical values, how differentiated the clusters identified with the previous analyses are.

### R code

> 📄 **R script:** `script_fst_allpop.R`
>
> Add/link the R script here once your GitHub folder structure is finalized.

---

## Questions

1. How many clusters are identified according to the **Mclust** analysis?
2. How many clusters are identified according to the **sNMF** analysis?
3. How do the clusters look according to **F<sub>ST</sub>**?

### Extra

You can have fun trying out many different approaches to study population structure, such as: Multidimensional scaling (MDS), Discriminant analysis of principal components (DAPC), STRUCTURE/ADMIXTURE, Phylogenetic trees, AMOVA

---

# 3. Create the SNP panel

We will use **three different approaches** to identify informative SNPs:

1. PCA loadings
2. F<sub>ST</sub>
3. Random Forest

Each method identifies informative markers according to a different criterion, allowing us to compare their performance and determine which approach is most suitable for the final SNP panel.

---

## 3a. Loadings from PCA

PCA can be used to identify SNPs that contribute most strongly to the genetic differentiation observed among individuals.

For each SNP, we will calculate its **loading** to the principal components that best capture the structure present in our dataset.

SNPs with the highest absolute loadings are considered the most informative because they have the strongest influence on the separation of individuals along the selected principal components.

We will therefore rank SNPs according to their loadings and select the most informative markers.

### R code

> 📄 **R script:** `script_loadings.R`
>
> Add/link the R script here once your GitHub folder structure is finalized.

### Extra

Have a look at the scatter plots including **PC1–PC2** and **PC3–PC4**.

> **Do you think PC2 and PC3 are useful to distinguish possible substructure in the dataset?**

---

## 3b. F<sub>ST</sub> by locus

F<sub>ST</sub> measures genetic differentiation between populations or groups.

It can be computed as a single value for the overall differentiation between two (or more) groups, but also at every marker in the dataset.

SNPs with high F<sub>ST</sub> values show a strong difference in allele frequencies between the groups and can therefore be particularly useful for distinguishing them.

We will calculate F<sub>ST</sub> for each SNP and rank the markers according to their level of genetic differentiation.

SNPs with the highest F<sub>ST</sub> values will be considered the most informative candidates for the SNP panel.

### R code

> 📄 **R script:** `script_fst_node.R`
>
> Add/link the R script here once your GitHub folder structure is finalized.

---

## 3c. Random Forest

**Random Forest** is a machine-learning approach used in classification and regression problems to make accurate predictions.

It can be used to identify SNPs that are particularly useful for discriminating between predefined groups.

The algorithm builds multiple decision trees using a subset of the available features (i.e. SNPs) and evaluates how informative each of them is for predicting group membership.

SNPs with higher importance scores contribute more strongly to the classification.

We will use Random Forest to rank SNPs according to their importance and select the most informative markers.

### R code

> 📄 **R script:** `script_rf.R`
>
> Add/link the R script here once your GitHub folder structure is finalized.

### Extra

If you have extra time, you can look for the SNPs useful to distinguish one of the next levels of differentiation within each subspecies:

* **GRE vs MES**
* **ITP + NCA + CCA vs SCA + SIC + SAR**

See the cladogram in the `Day3/tutorial_snpchip` folder for a summary of dataset subdivisions.

Choose the method you prefer, extract the samples belonging to the groups you want to compare, and get rid of the **invariable sites** before looking for the informative SNPs.

---

# 4. Test the accuracy in geographic assignment with samples of known origin

Now that we have created three different lists of informative SNPs according to different methodologies, we need to choose one of them based on their accuracy in assigning the correct geographic origin.

To do this, we will use **assignment tests**.

The idea is to provide the known origin of the samples, randomly divide each population into a **training set** and a **test set** for several iterations, and each time evaluate how well the test individuals were assigned to their corresponding cluster based on the training individuals from the same cluster.

### R code

> 📄 **R script:** `script_assignPOP.R`
>
> Add/link the R script here once your GitHub folder structure is finalized.

## Question

> **Based on the results from the assignment tests, which SNP panel would you choose?**

### Extra

If you have extra time, you can try providing different proportions of:

* Training individuals
* Loci

It can be particularly relevant to see if you can reduce the number of SNPs selected while keeping a good assignment accuracy, as this could reduce the costs of future genotyping of new samples.

---

# 5. Apply the SNP panel to new samples

This is where all our efforts finally will repay us!

Now you will be split into groups and you’ll work together to find out where your samples of unknown origin come from.

Notice their **assignment scores** to the different clusters and discuss how you would manage each sample based on them.

### R code

> 📄 **R script:** `script_assignPOP.R`
>
> Add/link the R script here once your GitHub folder structure is finalized.
