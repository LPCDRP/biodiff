% BIODIFF(1) 0.1
% Afif Elghraoui <aelghraoui@sdsu.edu>
% November 2015

# NAME

biodiff - compare biological sequences

# SYNOPSIS

**biodiff** *reference.fasta* *query.fasta* > *out.vcf*

# DESCRIPTION

**biodiff** is a simple utility originally developed for whole genome
comparisons. It uses **diff**(1), an implementation of Myer's
algorithm to find longest common substrings and determine the minimal difference
between the sequences. Output is in the Variant Call Format.

