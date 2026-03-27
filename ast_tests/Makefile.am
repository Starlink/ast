TESTS = ast_test
check_PROGRAMS = ast_test
ast_test_SOURCES = ast_test.c

#ast_test_LDADD = `ast_link`
# Expand ast_link to avoid libast_pass2, which causes problems for Solaris
ast_test_LDADD = @LIBPAL@ @LIBCMINPACK@ \
    $(top_builddir)/libast.la \
    $(top_builddir)/libast_grf_3.2.la \
    $(top_builddir)/libast_grf_5.6.la \
    $(top_builddir)/libast_grf_2.0.la \
    $(top_builddir)/libast_grf3d.la \
    $(top_builddir)/libast_err.la \
    -lm
