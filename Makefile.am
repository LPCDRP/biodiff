.DELETE_ON_ERROR:

dist_bin_SCRIPTS = biodiff \
                   udiff2vcf

PANDOC = @PANDOC@

if HAVE_PANDOC
man_MANS = biodiff.1 \
           udiff2vcf.1
endif

EXTRA_DIST = LICENSE $(addsuffix .1.md,$(dist_bin_SCRIPTS))

TESTS = test/integration-wrapper

AM_TESTS_ENVIRONMENT = PATH="$(srcdir):$(srcdir)/test/exec:$$PATH"; DATADIR="$(srcdir)/test/lambda-phage"; export PATH DATADIR;

%.1: %.1.md
	$(PANDOC) -s --from=markdown --to=man $< > $@

clean-local:
	$(RM) -r integration.*

# http://gnu-automake.7480.n7.nabble.com/DELETE-ON-ERROR-test-suite-log-and-PRECIOUS-td20856.html
.PRECIOUS: $(TEST_SUITE_LOG)
