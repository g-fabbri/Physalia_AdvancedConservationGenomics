# Installing GONE2

GONE2 is not bundled because it is under active development and should be compiled on the teaching platform.

From the course root:

~~~bash
git clone --depth 1 https://github.com/esrud/GONE2.git software/GONE2
cd software/GONE2
make gone
cd ../..
~~~

Confirm the installation:

~~~bash
test -x software/GONE2/gone2 && echo "GONE2 is ready"
~~~

Official source and documentation: https://github.com/esrud/GONE2
