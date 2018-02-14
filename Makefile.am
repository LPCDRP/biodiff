.DELETE_ON_ERROR:

dist_bin_SCRIPTS = biodiff \
                   udiff2vcf

PANDOC = @PANDOC@

if HAVE_PANDOC
man_MANS = biodiff.1 \
           udiff2vcf.1
endif

EXTRA_DIST = INSTALL INSTALL.md LICENSE $(addsuffix .1.md,$(dist_bin_SCRIPTS)) \
test/integration-wrapper test/exec/integration \
$(srcdir)/test/lambda-phage/*

TESTS = test/integration-wrapper

AM_TESTS_ENVIRONMENT = PATH="$(srcdir):$(srcdir)/test/exec:$$PATH"; DATADIR="$(srcdir)/test/lambda-phage"; export PATH DATADIR;

installcheck-local:
	DATADIR="$(srcdir)/test/lambda-phage" \
	PATH="$(srcdir)/test/exec:$$PATH" $(srcdir)/test/integration-wrapper

%.1: %.1.md
	sed '/^\[!/d' $< `# Remove conda badge before converting to man` \
	| $(PANDOC) -s --from=markdown --to=man > $@

# Convert markdown documentation to true plain text
# (no matter how much you try, markdown source won't really be as readable as formatted plain text)
INSTALL: INSTALL.md
	$(PANDOC) -t plain $< -o $@

clean-local:
	$(RM) -r integration.*
	$(RM) $(man_MANS)
	$(RM) INSTALL

# http://gnu-automake.7480.n7.nabble.com/DELETE-ON-ERROR-test-suite-log-and-PRECIOUS-td20856.html
.PRECIOUS: $(TEST_SUITE_LOG)
