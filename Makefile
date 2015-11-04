.DELETE_ON_ERROR:

include ../common.mk

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
