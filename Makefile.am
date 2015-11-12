.DELETE_ON_ERROR:

dist_bin_SCRIPTS = biodiff \
                   udiff2vcf

man_MANS = biodiff.1 \
           udiff2vcf.1

TESTS = $(foreach type,vcf patch, $(subst qry.fasta,$(type).sh,$(wildcard $(srcdir)/test/lambda-phage/*.qry.fasta)))

TESTS_ENVIRONMENT = PATH="$(srcdir):$$PATH"; export PATH;

%.1: %.1.md
	pandoc -s --from=markdown --to=man $< > $@

%.vcf.sh: %.qry.fasta
	echo "#!/bin/sh" > $@
	echo "biodiff test/lambda-phage/reference.fasta $< > $*.vcf" >> $@
	chmod +x $@

%.patch.sh: %.qry.fasta
	echo "#!/bin/sh" > $@
	echo "$(BGZIP) -f $*.vcf" >> $@
	echo "$(TABIX) -p vcf $*.vcf.gz" >> $@
	chmod +x $@

test/lambda-phage/%.vcf.gz: test/lambda-phage/%.vcf
	$(BGZIP) $<
	$(TABIX) -p vcf $@

test/lambda-phage/%.vcf: test/lambda-phage/%.qry.fasta
	biodiff test/lambda-phage/reference.fasta $< > $@

clean-local:
	$(RM) test/lambda-phage/*.vcf* test/lambda-phage/*.sh
