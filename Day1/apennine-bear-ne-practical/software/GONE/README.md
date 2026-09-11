# GONE software directory

This directory contains the Apennine brown bear driver and the classroom parameter values. The compiled GONE programs are platform-specific and are not currently included in this repository.

## Instructor setup on Linux

Before class, download the [official GONE repository](https://github.com/esrud/GONE) and copy the contents of its **Linux/PROGRAMMES** directory into:

```text
software/GONE/PROGRAMMES/
```

Required files include:

```text
MANAGE_CHROMOSOMES2
LD_SNP_REAL3
SUMM_REP_CHROM3
GONE
GONEaverage
GONEparallel.sh
```

From the course root, one possible installation is:

~~~bash
git clone --depth 1 https://github.com/esrud/GONE.git software/GONE_official
cp -a software/GONE_official/Linux/PROGRAMMES/. software/GONE/PROGRAMMES/
chmod u+x software/GONE/PROGRAMMES/*
~~~

Confirm the installation:

~~~bash
ls -l software/GONE/PROGRAMMES
test -x software/GONE/PROGRAMMES/MANAGE_CHROMOSOMES2
test -x software/GONE/PROGRAMMES/LD_SNP_REAL3
test -x software/GONE/PROGRAMMES/SUMM_REP_CHROM3
test -x software/GONE/PROGRAMMES/GONEparallel.sh
~~~

Test the complete directory on the same operating system used during the class. Every student has a separate course directory, so analyses can run sequentially inside their own **software/GONE** directory.
