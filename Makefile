.DELETE_ON_ERROR:

include common.mk

final.vcf: vcf-header out.vcf
	cat $^ > $@

out.vcf: out.diff
	./udiff2vcf -r $(REFERENCE_ID) < $< > $@

out.diff: $(subst .fasta,.txt,$(REFERENCE) $(QUERY))
# * Ignore the return code of diff because it returns 1 if the files
#   are different and that will be seen as an error by Make
# * The VCF requires one position of context for insertions/deletions
	-diff					\
		--ignore-case			\
		--unified=1			\
		--minimal			\
		$^ > $@

%.txt: %.fasta
	tail -n+2 $< | \
	perl -n -e 'chomp; print join("\n",split("",$$_))."\n"' > $@

# http://vcftools.sourceforge.net/perl_module.html
check: final.vcf.gz final.vcf.gz.tbi
	vcf-validator $<
	cat $(REFERENCE) | vcf-consensus $< > out.fasta
	dnadiff $(QUERY) out.fasta

final.vcf.gz: final.vcf
	bgzip $<

final.vcf.gz.tbi: final.vcf.gz
	tabix -p vcf $<

.PHONY: check
