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


## 4 — Install NeEstimator

NeEstimator is not installed through this Conda environment. Download the appropriate release from the [official NeEstimator page](https://www.molecularfisherieslaboratory.com/neestimator-software/) and follow its platform-specific instructions. If it is unavailable on the teaching cluster, provide precomputed output for the interpretation exercise.

## 5 — Verify everything

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

~~~

Run these checks on the same operating system and compute nodes used for teaching. The Conda environment manages shared dependencies, while the GONE2 implementation remain local under **software/**.
