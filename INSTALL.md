# Dependencies

* diff
* bash
* perl
* pandoc (to build documentation)

If you don't have pandoc, the documentation will not be built and installed.

To run the test suite, you will additionally need the following:

* tabix
* bcftools

# Installation

## Conda

biodiff is available as a conda package via bioconda so it may be installed as follows:

~~~
$ conda install -c bioconda biodiff
~~~

## Release Distribution

The release distribution is *not* the same as the source code downloads produced by Gitlab.
On the [Tags](https://gitlab.com/LPCDRP/biodiff/tags) page, this section's directions apply to the downloads that are named "biodiff-*version*.tar.gz".

~~~
./configure
make
make check
make install
~~~

## From Git or source code download

You will additionally need autoconf and automake installed.

~~~
autoreconf -i
~~~

Then you can follow the steps above for a release distribution.


# Testing the Installed Version

~~~
make installcheck
~~~
