% UDIFF2VCF(1) 0.2
% Afif Elghraoui
% February 2017

# NAME

udiff2vcf - convert unified diff of nucleotide sequences to Variant Call Format

# SYNOPSIS

**udiff2vcf** < *query.diff*

# DESCRIPTION

Genomes can be aligned to find unique exact matches using **diff**(1)
if they are written to files with one nucleotide per line. This program
converts the result of **diff -r --unified=1 --ignore-case** to the
variant call format for downstream bioinformatic analysis. Input is
read from stdin and output is written to stdout.

# EXAMPLE

~~~
diff --recursive --unified=1 --ignore-case refdir qrydir | \
     udiff2vcf > query.vcf
~~~
