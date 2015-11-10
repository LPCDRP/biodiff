% UDIFF2VCF(1) 0.1
% Afif Elghraoui
% November 2015

# NAME

udiff2vcf - convert unified diff of nucleotide sequences to Variant Call Format

# SYNOPSIS

**udiff2vcf** **-r** *ref_id* < *query.diff*

# DESCRIPTION

Genomes can be aligned to find unique exact matches using **diff**(1)
if they are written to files with one nucleotide per line. This program
converts the result of **diff --unified=1 --ignore-case** to the
variant call format for downstream bioinformatic analysis. Input is
read from stdin and output is written to stdout.

# ARGUMENTS

**-r** *ref_id*
:    The reference sequence ID. This will be the first part of the fasta
     header of the reference sequence file. This is important in order for
     VCF parsers to recognize the reference sequence.

# EXAMPLE

~~~
diff --unified=1 --ignore-case reference.txt query.txt | \
     udiff2vcf -r 'gi|448814763|ref|NC_000962.3|' > query.vcf
~~~
