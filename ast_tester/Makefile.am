# ---------------------------------------------------------------------------
# Standalone Fortran tests
#
# These build and run the Fortran test suite without requiring a Starlink
# installation.  The star_compat library provides minimal stubs for EMS
# (error reporting), PSX (memory allocation), and CHR (string handling).
# ---------------------------------------------------------------------------
if !NOFORTRAN

# Compatibility stubs replacing Starlink EMS/PSX/CHR for test-only use.
noinst_LTLIBRARIES = libstar_compat.la
libstar_compat_la_SOURCES = \
    star_compat/err_standalone.c \
    star_compat/psx_standalone.c \
    star_compat/chr_standalone.c \
    star_compat/cnf_standalone.c \
    star_compat/num_standalone.c
libstar_compat_la_CFLAGS = $(AM_CFLAGS) -I$(top_builddir)/src

# Fortran compiler include path for test include files:
#   - .. for AST_PAR (generated at build time)
#   - star_compat/ for SAE_PAR, SAI_ERR, ERR_PAR, PRM_PAR
#   - $(builddir) for AST_ERR and any other generated headers
# For .f (Fortran 77 fixed-form) sources, automake uses F77/FFLAGS.
# Per-target flags for .f files use the _FFLAGS suffix.
FTEST_FFLAGS = \
    -I$(top_builddir) \
    -I$(srcdir)/star_compat \
    -I$(builddir) \
    -fno-range-check \
    -fallow-argument-mismatch

# Libraries for all Fortran tests.
FTEST_LDADD = libstar_compat.la @LIBPAL@ @LIBCMINPACK@ \
    $(top_builddir)/libast.la \
    $(top_builddir)/libast_grf_3.2.la \
    $(top_builddir)/libast_grf_5.6.la \
    $(top_builddir)/libast_grf_2.0.la \
    $(top_builddir)/libast_grf3d.la \
    $(top_builddir)/libast_err.la \
    -lm

FORTRAN_TESTS = \
    testchannel \
    testchebymap \
    testcmpframe \
    testcmpmap \
    testfitschan \
    testfitstable \
    testflux \
    testframeset \
    testkeymap \
    testlutmap \
    testmapping \
    testmoc \
    testmocchan \
    testnormmap \
    testpolymap \
    testrate \
    testratemap \
    testrebin \
    testrebinseq \
    testregions \
    testskyframe \
    testspecflux \
    testspecframe \
    testsplinemap \
    teststc \
    teststcschan \
    testswitchmap \
    testtable \
    testtime \
    testtrangrid \
    testunitnormmap \
    testxmlchan \
    testyamlchan \
    testzoommap

TESTS = $(FORTRAN_TESTS)

check_PROGRAMS = $(FORTRAN_TESTS)

# Source and flag declarations for each Fortran test.
testchannel_SOURCES = testchannel.f
testchannel_FFLAGS  = $(FTEST_FFLAGS)
testchannel_LDADD   = $(FTEST_LDADD)

testchebymap_SOURCES = testchebymap.f
testchebymap_FFLAGS  = $(FTEST_FFLAGS)
testchebymap_LDADD   = $(FTEST_LDADD)

testcmpframe_SOURCES = testcmpframe.f
testcmpframe_FFLAGS  = $(FTEST_FFLAGS)
testcmpframe_LDADD   = $(FTEST_LDADD)

testcmpmap_SOURCES = testcmpmap.f
testcmpmap_FFLAGS  = $(FTEST_FFLAGS)
testcmpmap_LDADD   = $(FTEST_LDADD)

testfitschan_SOURCES = testfitschan.f
testfitschan_FFLAGS  = $(FTEST_FFLAGS)
testfitschan_LDADD   = $(FTEST_LDADD)

testfitstable_SOURCES = testfitstable.f
testfitstable_FFLAGS  = $(FTEST_FFLAGS)
testfitstable_LDADD   = $(FTEST_LDADD)

testflux_SOURCES = testflux.f
testflux_FFLAGS  = $(FTEST_FFLAGS)
testflux_LDADD   = $(FTEST_LDADD)

testframeset_SOURCES = testframeset.f
testframeset_FFLAGS  = $(FTEST_FFLAGS)
testframeset_LDADD   = $(FTEST_LDADD)

testkeymap_SOURCES = testkeymap.f
testkeymap_FFLAGS  = $(FTEST_FFLAGS)
testkeymap_LDADD   = $(FTEST_LDADD)

testlutmap_SOURCES = testlutmap.f
testlutmap_FFLAGS  = $(FTEST_FFLAGS)
testlutmap_LDADD   = $(FTEST_LDADD)

testmapping_SOURCES = testmapping.f
testmapping_FFLAGS  = $(FTEST_FFLAGS)
testmapping_LDADD   = $(FTEST_LDADD)

testmoc_SOURCES = testmoc.f
testmoc_FFLAGS  = $(FTEST_FFLAGS)
testmoc_LDADD   = $(FTEST_LDADD)

testmocchan_SOURCES = testmocchan.f
testmocchan_FFLAGS  = $(FTEST_FFLAGS)
testmocchan_LDADD   = $(FTEST_LDADD)

testnormmap_SOURCES = testnormmap.f
testnormmap_FFLAGS  = $(FTEST_FFLAGS)
testnormmap_LDADD   = $(FTEST_LDADD)

testpolymap_SOURCES = testpolymap.f
testpolymap_FFLAGS  = $(FTEST_FFLAGS)
testpolymap_LDADD   = $(FTEST_LDADD)

testrate_SOURCES = testrate.f
testrate_FFLAGS  = $(FTEST_FFLAGS)
testrate_LDADD   = $(FTEST_LDADD)

testratemap_SOURCES = testratemap.f
testratemap_FFLAGS  = $(FTEST_FFLAGS)
testratemap_LDADD   = $(FTEST_LDADD)

testrebin_SOURCES = testrebin.f
testrebin_FFLAGS  = $(FTEST_FFLAGS)
testrebin_LDADD   = $(FTEST_LDADD)

testrebinseq_SOURCES = testrebinseq.f
testrebinseq_FFLAGS  = $(FTEST_FFLAGS)
testrebinseq_LDADD   = $(FTEST_LDADD)

testregions_SOURCES = testregions.f
testregions_FFLAGS  = $(FTEST_FFLAGS)
testregions_LDADD   = $(FTEST_LDADD)

testskyframe_SOURCES = testskyframe.f
testskyframe_FFLAGS  = $(FTEST_FFLAGS)
testskyframe_LDADD   = $(FTEST_LDADD)

testspecflux_SOURCES = testspecflux.f
testspecflux_FFLAGS  = $(FTEST_FFLAGS)
testspecflux_LDADD   = $(FTEST_LDADD)

testspecframe_SOURCES = testspecframe.f
testspecframe_FFLAGS  = $(FTEST_FFLAGS)
testspecframe_LDADD   = $(FTEST_LDADD)

testsplinemap_SOURCES = testsplinemap.f
testsplinemap_FFLAGS  = $(FTEST_FFLAGS)
testsplinemap_LDADD   = $(FTEST_LDADD)

teststc_SOURCES = teststc.f
teststc_FFLAGS  = $(FTEST_FFLAGS)
teststc_LDADD   = $(FTEST_LDADD)

teststcschan_SOURCES = teststcschan.f
teststcschan_FFLAGS  = $(FTEST_FFLAGS)
teststcschan_LDADD   = $(FTEST_LDADD)

testswitchmap_SOURCES = testswitchmap.f
testswitchmap_FFLAGS  = $(FTEST_FFLAGS)
testswitchmap_LDADD   = $(FTEST_LDADD)

testtable_SOURCES = testtable.f
testtable_FFLAGS  = $(FTEST_FFLAGS)
testtable_LDADD   = $(FTEST_LDADD)

testtime_SOURCES = testtime.f
testtime_FFLAGS  = $(FTEST_FFLAGS)
testtime_LDADD   = $(FTEST_LDADD)

testtrangrid_SOURCES = testtrangrid.f
testtrangrid_FFLAGS  = $(FTEST_FFLAGS)
testtrangrid_LDADD   = $(FTEST_LDADD)

testunitnormmap_SOURCES = testunitnormmap.f
testunitnormmap_FFLAGS  = $(FTEST_FFLAGS)
testunitnormmap_LDADD   = $(FTEST_LDADD)

testxmlchan_SOURCES = testxmlchan.f
testxmlchan_FFLAGS  = $(FTEST_FFLAGS)
testxmlchan_LDADD   = $(FTEST_LDADD)

testyamlchan_SOURCES = testyamlchan.f
testyamlchan_FFLAGS  = $(FTEST_FFLAGS)
testyamlchan_LDADD   = $(FTEST_LDADD)

testzoommap_SOURCES = testzoommap.f
testzoommap_FFLAGS  = $(FTEST_FFLAGS)
testzoommap_LDADD   = $(FTEST_LDADD)

endif  # NOFORTRAN
