# Software setup

Complete this setup before class. Students should only activate the prepared environment and follow the lessons.

## 1 — Create one Conda environment

From the course root:

~~~bash
conda env create -f environment.yml
conda activate bear-ne-practical
~~~

This environment supplies BCFtools, BEDTools, tabix/bgzip, PLINK 1.9, MSMC2, R, Python, Git, Make, and a C++ compiler. GONE2 and msmc-tools require the additional installation steps below. SMC++ and currentNe2 are optional instructor-tested extensions.

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

## 4 — Optional SMC++ installation

Install the [Bioconda SMC++ package](https://anaconda.org/bioconda/smcpp) in the existing course environment before offering the [optional lesson](../lessons/03-smcpp-optional.md):

~~~bash
conda activate bear-ne-practical
conda install bioconda::smcpp
conda list smcpp
smc++ vcf2smc -h | head
smc++ estimate -h | head
smc++ plot -h | head
~~~

This is an optional install, so it is not included in `environment.yml`. Record the installed version and test all lesson commands on the teaching compute node. The [upstream project](https://github.com/popgenmethods/smcpp) recommends a container for its own latest release; the Conda package may not be that exact build.

## 5 — Optional currentNe2 installation

Compile the [official currentNe2 repository](https://github.com/esrud/currentNe2) on the teaching Linux node:

~~~bash
git clone --depth 1 https://github.com/esrud/currentNe2.git software/currentNe2
make -C software/currentNe2
test -x software/currentNe2/currentne2 && echo 'currentNe2 ready'
~~~

Check the compiler output and `software/currentNe2/currentne2 -h` before using the [optional lesson](../lessons/03b-currentne2-optional.md). The one-scaffold teaching input is useful to compare software output, not to establish a reliable contemporary population size.

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

~~~

Run these checks on the same operating system and compute nodes used for teaching. The Conda environment manages shared dependencies, while GONE2 and optional currentNe2 remain local under **software/**.
