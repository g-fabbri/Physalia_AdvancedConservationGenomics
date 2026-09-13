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

The [SMC++ project](https://github.com/popgenmethods/smcpp) currently recommends a versioned container; its latest release is **not** supplied by this Conda environment. On Jarvis, ask the cluster administrator whether Apptainer/Singularity can run that container, or follow the project's source-build instructions in a separate tested environment. Do not assume `conda install smcpp` provides the current release. Before offering the [optional lesson](../lessons/03-smcpp-optional.md), confirm that `smc++ vcf2smc -h`, `smc++ estimate -h`, and `smc++ plot -h` work on the teaching compute node, and record the exact version and runtime.

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
