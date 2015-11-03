% UDIFF2VCF(1) 0.1
% Afif Elghraoui
% November 2015

# NAME

udiff2vcf - convert unified diff of nucleotide sequences to vcf

# SYNOPSIS

**udiff2vcf** < *query.diff*

# DESCRIPTION

Genomes can be aligned to find unique exact matches using **diff**(1)
if they are written to files with one nucleotide per line. This program
converts the result of **diff --unified=1 --ignore-case** to the
variant calling format for downstream bioinformatic analysis. Input is
read from stdin and output is written to stdout.

# BUGS

* The output may contain fencepost errors.
* Insertions and the deletions at the beginning and end of the sequence will
  have a base that is incorrectly treated as context.

# EXAMPLE

~~~
diff --unified=1 --ignore-case reference.txt query.txt | udiff2vcf > query.vcf
~~~
