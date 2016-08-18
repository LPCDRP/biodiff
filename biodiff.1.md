% BIODIFF(1) 0.1.0
% Afif Elghraoui <aelghraoui@sdsu.edu>
% August 2016

# NAME

biodiff - compare biological sequences

# SYNOPSIS

**biodiff** *reference.fasta* *query.fasta* > *out.vcf*

# DESCRIPTION

**biodiff** determines the exact differences between two biological sequences.
It can operate on DNA and protein sequences, as long as they are in fasta format.
It uses **diff**(1), an implementation of Myer's algorithm to find longest common substrings and determine the minimal difference between the sequences.
Output is in the Variant Call Format.

**biodiff** is especially useful for exact genome comparison, as standard genome comparison tools are often vague regarding the positions of large insertions and deletions.
It can be helpful to first get an accurate picture of the plain insertions and deletions that differentiate two sequences, before trying to decide whether they represent translocations, tandem copy number variation, or anything else.

# EXAMPLES

## Quick evaluation
You might want to quickly see the difference between two revisions of an NCBI reference sequence.

~~~
biodiff NC_123456.1.fasta NC_123456.2.fasta
~~~

The output goes to your terminal.

## Protein Sequence Comparison

It works the same way with amino acid sequences:

~~~
biodiff wild-type.faa mutant.faa
~~~

## Comparison of variants between samples
You will typically want to normalize and left-align variants, especially when comparing variants between different samples, as variants within repetitive sequences can be accurately represented in multiple ways and falsely appear to be different.

~~~
for query in query1 query2
do
	biodiff reference.fasta $query.fasta | bcftools norm --fasta-ref reference.fasta - > $query.vcf
	bgzip $query.vcf && tabix -p vcf $query.vcf.gz
done

bcftools isec query1.vcf.gz query2.vcf.gz
~~~

# BUGS

Please see our issue tracker for known issues.
https://github.com/valafarlab/biodiff/issues

# SEE ALSO

**bcftools**(1)
**dnadiff**(1)
