# Software setup

Complete this setup before class. Students should only activate the prepared environment and follow the lessons.

## 1 — Create one Conda environment

From the course root:

~~~bash
conda env create -f environment.yml
conda activate bear-ne-practical
~~~

This environment supplies BCFtools, tabix/bgzip, PLINK 1.9, MSMC2, R, Python, Git, Make, and a C++ compiler. GONE, GONE2, msmc-tools, and NeEstimator require the additional installation steps below.

## 2 — Install msmc-tools

The helper scripts, including **generate_multihetsep.py** and **multihetsep_bootstrap.py**, are distributed separately from MSMC2:

~~~bash
git clone --depth 1 https://github.com/stschiff/msmc-tools.git \
  software/msmc-tools
~~~

Add them to the active shell:

~~~bash
export PATH="$PWD/software/msmc-tools:$PATH"
~~~

## 3 — Compile GONE2

~~~bash
git clone --depth 1 https://github.com/esrud/GONE2.git software/GONE2
make -C software/GONE2 gone
~~~

## 4 — Install original GONE on Linux

The original GONE programs are platform-specific. Download the official repository and copy its Linux programs into the course driver directory:

~~~bash
git clone --depth 1 https://github.com/esrud/GONE.git \
  software/GONE_official

cp -a software/GONE_official/Linux/PROGRAMMES/. \
  software/GONE/PROGRAMMES/

chmod u+x software/GONE/PROGRAMMES/*
~~~

The optional original GONE lesson uses the supplied **script_GONE.sh** and **INPUT_PARAMETERS_FILE**, not the copies in **GONE_official/Linux**.

## 5 — Install NeEstimator

NeEstimator is not installed through this Conda environment. Download the appropriate release from the [official NeEstimator page](https://www.molecularfisherieslaboratory.com/neestimator-software/) and follow its platform-specific instructions. If it is unavailable on the teaching cluster, provide precomputed output for the interpretation exercise.

## 6 — Verify everything

~~~bash
conda activate bear-ne-practical

bcftools --version | head -n 1
plink --version
msmc2 --help | head -n 1
Rscript --version

generate_multihetsep.py --help | head
multihetsep_bootstrap.py --help | head

test -x software/GONE2/gone2 &&
  echo "GONE2 ready"

test -x software/GONE/PROGRAMMES/MANAGE_CHROMOSOMES2 &&
test -x software/GONE/PROGRAMMES/LD_SNP_REAL3 &&
test -x software/GONE/PROGRAMMES/SUMM_REP_CHROM3 &&
test -x software/GONE/PROGRAMMES/GONEparallel.sh &&
  echo "Original GONE ready"
~~~

Run these checks on the same operating system and compute nodes used for teaching. The Conda environment manages shared dependencies, while the two GONE implementations remain local under **software/**.
