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

## Release Distribution

The release distribution is *not* the same as the "Source Code" downloads produced by Github.
On the [Releases](https://github.com/valafarlab/biodiff/releases) page, this section's directions apply to the downloads that are named "biodiff-*version*.tar.gz".

~~~
./configure
make
make check
make install
~~~

## From Git or "Source Code" download

You will additionally need autoconf and automake installed.

~~~
autoreconf -i
~~~

Then you can follow the steps above for a release distribution.
