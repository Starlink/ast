# Performance benchmarks.  These are not tests and are not run by the build;
# they are built by "make check" so that distcheck keeps them compiling.  Run
# them by hand from the build tree, e.g. "./bench_lock".

AM_CPPFLAGS = -DHAVE_CONFIG_H \
    -I$(top_builddir) -I$(top_srcdir) \
    -I$(top_builddir)/src -I$(top_srcdir)/src

check_PROGRAMS = bench_lock

bench_lock_SOURCES = bench_lock.c
bench_lock_LDADD = $(AST_PLUGIN_KEEP) @LIBPAL@ @LIBCMINPACK@ \
    $(top_builddir)/libast.la \
    $(top_builddir)/libast_grf_3.2.la \
    $(top_builddir)/libast_grf_5.6.la \
    $(top_builddir)/libast_grf_2.0.la \
    $(top_builddir)/libast_grf3d.la \
    $(top_builddir)/libast_err.la \
    -lm

EXTRA_DIST = CMakeLists.txt README.md
