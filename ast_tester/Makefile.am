# ---------------------------------------------------------------------------
# AST test suite
#
# C tests run unconditionally.  Fortran tests are gated on !NOFORTRAN.
#
# In standalone mode (no Starlink), star_compat/ provides C header stubs
# (sae_par.h, sai_err.h) and a Fortran stub library (libstar_compat.la)
# for the subset of the Starlink API used by the tests.
#
# In Starlink mode (starconf detected), the real Starlink headers and
# libraries are used.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# C tests
# ---------------------------------------------------------------------------

AM_CPPFLAGS = -DHAVE_CONFIG_H -I$(top_builddir) -I$(top_srcdir)/src \
    $(STAR_CPPFLAGS)

if !STARLINK_BUILD
AM_CPPFLAGS += -I$(srcdir)/star_compat
endif

CTEST_LDADD = @LIBPAL@ @LIBCMINPACK@ \
    $(top_builddir)/libast.la \
    $(top_builddir)/libast_grf_3.2.la \
    $(top_builddir)/libast_grf_5.6.la \
    $(top_builddir)/libast_grf_2.0.la \
    $(top_builddir)/libast_grf3d.la \
    $(top_builddir)/libast_err.la \
    -lm

C_TESTS = \
    ast_test \
    testaxis \
    testconvert \
    testerror \
    testframe \
    testobject \
    testresimp \
    testsplinemap_c \
    testunitnorm

TESTS = $(C_TESTS)
check_PROGRAMS = $(C_TESTS)

ast_test_SOURCES       = ast_test.c
ast_test_LDADD         = $(CTEST_LDADD)

testaxis_SOURCES       = testaxis.c
testaxis_LDADD         = $(CTEST_LDADD)

testconvert_SOURCES    = testconvert.c
testconvert_LDADD      = $(CTEST_LDADD)

testerror_SOURCES      = testerror.c
testerror_LDADD        = $(CTEST_LDADD)

testframe_SOURCES      = testframe.c
testframe_LDADD        = $(CTEST_LDADD)

testobject_SOURCES     = testobject.c
testobject_LDADD       = $(CTEST_LDADD)

testresimp_SOURCES     = testresimp.c
testresimp_LDADD       = $(CTEST_LDADD)

testsplinemap_c_SOURCES = testsplinemap_c.c
testsplinemap_c_LDADD  = $(CTEST_LDADD)

testunitnorm_SOURCES   = testunitnorm.c
testunitnorm_LDADD     = $(CTEST_LDADD)

# ---------------------------------------------------------------------------
# Fortran tests
# ---------------------------------------------------------------------------
if !NOFORTRAN

if STARLINK_BUILD  # Starlink mode: use real Starlink libraries

# Fortran compiler flags for test include files:
#   - $(top_builddir) for AST_PAR (generated at build time)
#   - $(builddir) for AST_ERR and any other generated headers
#   SAE_PAR, ERR_PAR, PRM_PAR, etc. come from $STARLINK/include via STAR_FCFLAGS.
#   -fallow-argument-mismatch is already in STAR_FCFLAGS; only add -fno-range-check
#   and the local build directory.
AM_FCFLAGS = \
    -I$(top_builddir) \
    -I$(builddir) \
    -fno-range-check

# Cursed: Starlink's patched automake forcibly inserts $(STAR_FCFLAGS) at the
# beginning of the flag list for $(FC), and likewise $(STAR_LDFLAGS) and
# $(FCLD). That's fine for the original use case of running the `ast_tester`
# script against an installed AST, but counter to the goals of an in-source
# `make check`. The most reasonable workaround seems to be to override the
# default pattern rule for compiling Fortran sources.
#
# For linking this does not appear to be a problem since locally built
# libraries are passed in by relative path and are handled by libtool.
.f.o:
	$(AM_V_FC)$(FC) $(AM_FCFLAGS) $(FCFLAGS) $(STAR_FCFLAGS) -c -o $@ $<

# Libraries for all Fortran tests (real Starlink libs).
# We link the locally-built AST .la files directly (not via ast_link, which
# targets an installed AST).  The Starlink support libraries mirror what the
# ast_tester shell script obtains from:
#   `err_link` `psx_link` `prm_link` `chr_link`
# These flags must be hard-coded here because $(shell ...) in LDADD causes
# automake to add the expansion to _DEPENDENCIES, making make try to build
# the -l flags as file targets.
FTEST_LDADD = @LIBPAL@ @LIBCMINPACK@ \
    $(top_builddir)/libast.la \
    $(top_builddir)/libast_grf_3.2.la \
    $(top_builddir)/libast_grf_5.6.la \
    $(top_builddir)/libast_grf_2.0.la \
    $(top_builddir)/libast_grf3d.la \
    $(top_builddir)/libast_err.la \
    -lerr_standalone -lchr -lemsf -lems -lcnf -lstarmem -lpthread -lstarutil \
    -lpsx -lprm -lprm_a -lm
else
# Compatibility stubs replacing Starlink EMS/PSX/CHR/CNF/NUM for test-only use.
noinst_LTLIBRARIES = libstar_compat.la
libstar_compat_la_SOURCES = \
    star_compat/err_standalone.c \
    star_compat/psx_standalone.c \
    star_compat/chr_standalone.c \
    star_compat/cnf_standalone.c \
    star_compat/num_standalone.c
libstar_compat_la_CFLAGS = $(AM_CFLAGS) -I$(top_builddir)/src

# Fortran compiler include path for test include files:
#   - $(top_builddir) for AST_PAR (generated at build time)
#   - star_compat/ for SAE_PAR, SAI_ERR, ERR_PAR, PRM_PAR
#   - $(builddir) for AST_ERR and any other generated headers
# Starlink's automake compiles .f files via FC (not F77), so AM_FCFLAGS is
# used (it appears in FCCOMPILE = $(FC) $(STAR_FCFLAGS) $(AM_FCFLAGS) $(FCFLAGS)).
AM_FCFLAGS = \
    -I$(top_builddir) \
    -I$(srcdir)/star_compat \
    -I$(builddir) \
    -fno-range-check \
    -fallow-argument-mismatch

# Automake actually defaults to the f77 compiler interface for bare .f files
# so we actually need to set AM_FFLAGS ; but keep AM_FCFLAGS as well in case
# the use of gfortran was forced
AM_FFLAGS = $(AM_FCFLAGS)

# Libraries for all Fortran tests (standalone stubs).
FTEST_LDADD = libstar_compat.la @LIBPAL@ @LIBCMINPACK@ \
    $(top_builddir)/libast.la \
    $(top_builddir)/libast_grf_3.2.la \
    $(top_builddir)/libast_grf_5.6.la \
    $(top_builddir)/libast_grf_2.0.la \
    $(top_builddir)/libast_grf3d.la \
    $(top_builddir)/libast_err.la \
    -lm
endif  # STARLINK_BUILD

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

TESTS += $(FORTRAN_TESTS)
check_PROGRAMS += $(FORTRAN_TESTS)

# testhuge is always compiled but skipped unless TEST_HUGE=1.
# The wrapper script exits 77 (automake SKIP) when TEST_HUGE is unset.
check_PROGRAMS += testhuge
testhuge_SOURCES = testhuge.f
testhuge_LDADD   = $(FTEST_LDADD)
TESTS += testhuge.sh
EXTRA_DIST = testhuge.sh

# Source and LDADD declarations for each Fortran test.
testchannel_SOURCES    = testchannel.f
testchannel_LDADD      = $(FTEST_LDADD)

testchebymap_SOURCES   = testchebymap.f
testchebymap_LDADD     = $(FTEST_LDADD)

testcmpframe_SOURCES   = testcmpframe.f
testcmpframe_LDADD     = $(FTEST_LDADD)

testcmpmap_SOURCES     = testcmpmap.f
testcmpmap_LDADD       = $(FTEST_LDADD)

testfitschan_SOURCES   = testfitschan.f
testfitschan_LDADD     = $(FTEST_LDADD)

testfitstable_SOURCES  = testfitstable.f
testfitstable_LDADD    = $(FTEST_LDADD)

testflux_SOURCES       = testflux.f
testflux_LDADD         = $(FTEST_LDADD)

testframeset_SOURCES   = testframeset.f
testframeset_LDADD     = $(FTEST_LDADD)

testkeymap_SOURCES     = testkeymap.f
testkeymap_LDADD       = $(FTEST_LDADD)

testlutmap_SOURCES     = testlutmap.f
testlutmap_LDADD       = $(FTEST_LDADD)

testmapping_SOURCES    = testmapping.f
testmapping_LDADD      = $(FTEST_LDADD)

testmoc_SOURCES        = testmoc.f
testmoc_LDADD          = $(FTEST_LDADD)

testmocchan_SOURCES    = testmocchan.f
testmocchan_LDADD      = $(FTEST_LDADD)

testnormmap_SOURCES    = testnormmap.f
testnormmap_LDADD      = $(FTEST_LDADD)

testpolymap_SOURCES    = testpolymap.f
testpolymap_LDADD      = $(FTEST_LDADD)

testrate_SOURCES       = testrate.f
testrate_LDADD         = $(FTEST_LDADD)

testratemap_SOURCES    = testratemap.f
testratemap_LDADD      = $(FTEST_LDADD)

testrebin_SOURCES      = testrebin.f
testrebin_LDADD        = $(FTEST_LDADD)

testrebinseq_SOURCES   = testrebinseq.f
testrebinseq_LDADD     = $(FTEST_LDADD)

testregions_SOURCES    = testregions.f
testregions_LDADD      = $(FTEST_LDADD)

testskyframe_SOURCES   = testskyframe.f
testskyframe_LDADD     = $(FTEST_LDADD)

testspecflux_SOURCES   = testspecflux.f
testspecflux_LDADD     = $(FTEST_LDADD)

testspecframe_SOURCES  = testspecframe.f
testspecframe_LDADD    = $(FTEST_LDADD)

testsplinemap_SOURCES  = testsplinemap.f
testsplinemap_LDADD    = $(FTEST_LDADD)

teststc_SOURCES        = teststc.f
teststc_LDADD          = $(FTEST_LDADD)

teststcschan_SOURCES   = teststcschan.f
teststcschan_LDADD     = $(FTEST_LDADD)

testswitchmap_SOURCES  = testswitchmap.f
testswitchmap_LDADD    = $(FTEST_LDADD)

testtable_SOURCES      = testtable.f
testtable_LDADD        = $(FTEST_LDADD)

testtime_SOURCES       = testtime.f
testtime_LDADD         = $(FTEST_LDADD)

testtrangrid_SOURCES   = testtrangrid.f
testtrangrid_LDADD     = $(FTEST_LDADD)

testunitnormmap_SOURCES = testunitnormmap.f
testunitnormmap_LDADD  = $(FTEST_LDADD)

testxmlchan_SOURCES    = testxmlchan.f
testxmlchan_LDADD      = $(FTEST_LDADD)

testyamlchan_SOURCES   = testyamlchan.f
testyamlchan_LDADD     = $(FTEST_LDADD)

testzoommap_SOURCES    = testzoommap.f
testzoommap_LDADD      = $(FTEST_LDADD)

endif  # NOFORTRAN
