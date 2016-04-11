% BIODIFF(1) 0.1.0
% Afif Elghraoui <aelghraoui@sdsu.edu>
% April 2016

# NAME

biodiff - compare biological sequences

# SYNOPSIS

**biodiff** *reference.fasta* *query.fasta* > *out.vcf*

# DESCRIPTION

Genome comparison tools are often vague regarding the positions of large insertions and deletions.
It can be helpful to first get an accurate picture of the plain insertions and deletions that differentiate two sequences, before trying to decide whether they represent translations, tandem copy number variation, or anything else.

**biodiff** is a simple utility developed for such whole genome comparisons.
It uses **diff**(1), an implementation of Myer's algorithm to find longest common substrings and determine the minimal difference between the sequences.
Output is in the Variant Call Format.

# EXAMPLES

## Quick evaluation
You might want to quickly see the difference between two revisions of an NCBI reference sequence.

~~~
biodiff NC_123456.1.fasta NC_123456.2.fasta
~~~

The output goes to your terminal.

## Quick comparisons
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

* **biodiff**, despite the name, is currently not biologically-aware and knows nothing of reverse-complemented sequences.
If the genomes you are comparing have inversions with respect to each other, those regions will be shown as having many small scattered variants.

* There still seems to be some whitespace issues with regard to input fasta files, like missing a newline at the end of your fasta file.
If you get an error when doing a comparison, it may be due to such issues.
A quick workaround is to swap the input arguments, thereby using the query as the reference instead.
That may or may not help.

* The VCF header does not properly set the chromosome name.
It currently uses the fasta header as such (rather than assigning it the proper number and adding that annotation to the VCF header).
This doesn't create any issues for using the output with **bcftools**(1) or **vcftools**(1), however.

* **biodiff** hasn't been tested with multi-fasta files and probably doesn't handle them correctly.

# SEE ALSO

**bcftools**(1)
