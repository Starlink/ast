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

AM_CPPFLAGS = -DHAVE_CONFIG_H \
    -I$(top_builddir) -I$(top_builddir)/src -I$(top_srcdir)/src \
    $(STAR_CPPFLAGS)

if !STARLINK_BUILD
AM_CPPFLAGS += -I$(srcdir)/star_compat
endif

AM_TESTS_ENVIRONMENT = srcdir='$(srcdir)'; export srcdir;

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
    testchannel_c \
    testchebymap_c \
    testcmpframe_c \
    testcmpmap_c \
    testconvert \
    testerror \
    testfitschan_c \
    testfitstable_c \
    testflux_c \
    testframe \
    testframeset_c \
    testkeymap_c \
    testlutmap_c \
    testmapping_c \
    testmoc_c \
    testmocchan_c \
    testnormmap_c \
    testobject \
    testpolymap_c \
    testrate_c \
    testratemap_c \
    testrebin_c \
    testrebinseq_c \
    testregions_c \
    testresimp \
    testskyframe_c \
    testspecflux_c \
    testspecframe_c \
    testsplinemap_c \
    teststc_c \
    teststcschan_c \
    testswitchmap_c \
    testtable_c \
    testtime_c \
    testtrangrid_c \
    testunitnorm \
    testunitnormmap_c \
    testxmlchan_c \
    testyamlchan \
    testzoommap_c

# testhuge (C version) is compiled unconditionally; testhuge_f is
# compiled when !NOFORTRAN.  Both are driven by testhuge.sh.
check_PROGRAMS = $(C_TESTS) testhuge_c
TESTS = $(C_TESTS)

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

testyamlchan_SOURCES   = testyamlchan.c
testyamlchan_LDADD     = $(CTEST_LDADD)

# Fortran-converted C tests.  Every one shares its source basename with a
# .f counterpart, so each needs a per-target _CFLAGS to force automake to
# use long object names (foo-foo.o) and avoid object-file collisions.

testchannel_c_SOURCES  = testchannel.c
testchannel_c_CFLAGS   = $(AM_CFLAGS)
testchannel_c_LDADD    = $(CTEST_LDADD)

testchebymap_c_SOURCES = testchebymap.c
testchebymap_c_CFLAGS  = $(AM_CFLAGS)
testchebymap_c_LDADD   = $(CTEST_LDADD)

testcmpframe_c_SOURCES = testcmpframe.c
testcmpframe_c_CFLAGS  = $(AM_CFLAGS)
testcmpframe_c_LDADD   = $(CTEST_LDADD)

testcmpmap_c_SOURCES   = testcmpmap.c
testcmpmap_c_CFLAGS    = $(AM_CFLAGS)
testcmpmap_c_LDADD     = $(CTEST_LDADD)

testfitschan_c_SOURCES = testfitschan.c
testfitschan_c_CFLAGS  = $(AM_CFLAGS)
testfitschan_c_LDADD   = $(CTEST_LDADD)

testfitstable_c_SOURCES = testfitstable.c
testfitstable_c_CFLAGS  = $(AM_CFLAGS)
testfitstable_c_LDADD   = $(CTEST_LDADD)

testflux_c_SOURCES     = testflux.c
testflux_c_CFLAGS      = $(AM_CFLAGS)
testflux_c_LDADD       = $(CTEST_LDADD)

testframeset_c_SOURCES = testframeset.c
testframeset_c_CFLAGS  = $(AM_CFLAGS)
testframeset_c_LDADD   = $(CTEST_LDADD)

testkeymap_c_SOURCES   = testkeymap.c
testkeymap_c_CFLAGS    = $(AM_CFLAGS)
testkeymap_c_LDADD     = $(CTEST_LDADD)

testlutmap_c_SOURCES   = testlutmap.c
testlutmap_c_CFLAGS    = $(AM_CFLAGS)
testlutmap_c_LDADD     = $(CTEST_LDADD)

testmapping_c_SOURCES  = testmapping.c
testmapping_c_CFLAGS   = $(AM_CFLAGS)
testmapping_c_LDADD    = $(CTEST_LDADD)

testmoc_c_SOURCES      = testmoc.c
testmoc_c_CFLAGS       = $(AM_CFLAGS)
testmoc_c_LDADD        = $(CTEST_LDADD)

testmocchan_c_SOURCES  = testmocchan.c
testmocchan_c_CFLAGS   = $(AM_CFLAGS)
testmocchan_c_LDADD    = $(CTEST_LDADD)

testnormmap_c_SOURCES  = testnormmap.c
testnormmap_c_CFLAGS   = $(AM_CFLAGS)
testnormmap_c_LDADD    = $(CTEST_LDADD)

testpolymap_c_SOURCES  = testpolymap.c
testpolymap_c_CFLAGS   = $(AM_CFLAGS)
testpolymap_c_LDADD    = $(CTEST_LDADD)

testrate_c_SOURCES     = testrate.c
testrate_c_CFLAGS      = $(AM_CFLAGS)
testrate_c_LDADD       = $(CTEST_LDADD)

testratemap_c_SOURCES  = testratemap.c
testratemap_c_CFLAGS   = $(AM_CFLAGS)
testratemap_c_LDADD    = $(CTEST_LDADD)

testrebin_c_SOURCES    = testrebin.c
testrebin_c_CFLAGS     = $(AM_CFLAGS)
testrebin_c_LDADD      = $(CTEST_LDADD)

testrebinseq_c_SOURCES = testrebinseq.c
testrebinseq_c_CFLAGS  = $(AM_CFLAGS)
testrebinseq_c_LDADD   = $(CTEST_LDADD)

testregions_c_SOURCES  = testregions.c
testregions_c_CFLAGS   = $(AM_CFLAGS)
testregions_c_LDADD    = $(CTEST_LDADD)

testskyframe_c_SOURCES = testskyframe.c
testskyframe_c_CFLAGS  = $(AM_CFLAGS)
testskyframe_c_LDADD   = $(CTEST_LDADD)

testspecflux_c_SOURCES = testspecflux.c
testspecflux_c_CFLAGS  = $(AM_CFLAGS)
testspecflux_c_LDADD   = $(CTEST_LDADD)

testspecframe_c_SOURCES = testspecframe.c
testspecframe_c_CFLAGS  = $(AM_CFLAGS)
testspecframe_c_LDADD   = $(CTEST_LDADD)

teststc_c_SOURCES      = teststc.c
teststc_c_CFLAGS       = $(AM_CFLAGS)
teststc_c_LDADD        = $(CTEST_LDADD)

teststcschan_c_SOURCES = teststcschan.c
teststcschan_c_CFLAGS  = $(AM_CFLAGS)
teststcschan_c_LDADD   = $(CTEST_LDADD)

testswitchmap_c_SOURCES = testswitchmap.c
testswitchmap_c_CFLAGS  = $(AM_CFLAGS)
testswitchmap_c_LDADD   = $(CTEST_LDADD)

testtable_c_SOURCES    = testtable.c
testtable_c_CFLAGS     = $(AM_CFLAGS)
testtable_c_LDADD      = $(CTEST_LDADD)

testtime_c_SOURCES     = testtime.c
testtime_c_CFLAGS      = $(AM_CFLAGS)
testtime_c_LDADD       = $(CTEST_LDADD)

testtrangrid_c_SOURCES = testtrangrid.c
testtrangrid_c_CFLAGS  = $(AM_CFLAGS)
testtrangrid_c_LDADD   = $(CTEST_LDADD)

testunitnormmap_c_SOURCES = testunitnormmap.c
testunitnormmap_c_CFLAGS  = $(AM_CFLAGS)
testunitnormmap_c_LDADD   = $(CTEST_LDADD)

testhuge_c_SOURCES     = testhuge.c
testhuge_c_CFLAGS      = $(AM_CFLAGS)
testhuge_c_LDADD       = $(CTEST_LDADD)

testxmlchan_c_SOURCES  = testxmlchan.c
testxmlchan_c_CFLAGS   = $(AM_CFLAGS)
testxmlchan_c_LDADD    = $(CTEST_LDADD)

testzoommap_c_SOURCES  = testzoommap.c
testzoommap_c_CFLAGS   = $(AM_CFLAGS)
testzoommap_c_LDADD    = $(CTEST_LDADD)


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
    $(AST_FFLAGS_EXTRA)

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
    $(AST_FFLAGS_EXTRA)

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
    testchannel_f \
    testchebymap_f \
    testcmpframe_f \
    testcmpmap_f \
    testfitschan_f \
    testfitstable_f \
    testflux_f \
    testframeset_f \
    testkeymap_f \
    testlutmap_f \
    testmapping_f \
    testmoc_f \
    testmocchan_f \
    testnormmap_f \
    testpolymap_f \
    testrate_f \
    testratemap_f \
    testrebin_f \
    testrebinseq_f \
    testregions_f \
    testskyframe_f \
    testspecflux_f \
    testspecframe_f \
    testsplinemap_f \
    teststc_f \
    teststcschan_f \
    testswitchmap_f \
    testtable_f \
    testtime_f \
    testtrangrid_f \
    testunitnormmap_f \
    testxmlchan_f \
    testzoommap_f

TESTS += $(FORTRAN_TESTS)
check_PROGRAMS += $(FORTRAN_TESTS)

check_PROGRAMS += testhuge_f
testhuge_f_SOURCES = testhuge.f
testhuge_f_LDADD   = $(FTEST_LDADD)
TESTS += testhuge.sh
EXTRA_DIST = testhuge.sh test_wcsconverter.sh test_simplify.sh numdiff.py

# Include all test fixture files also in EXTRA_DIST
EXTRA_DIST += \
    2dspline_c.dat \
    2dspline.dat \
    alt.header \
    image_wcs.asdf \
    lsst_wcs.txt \
    sip.head \
    sparse.ast \
    splittest1.ast \
    stcchan-test1-doc3.ast \
    stcchan-test1-doc3-props.ast \
    tanSipWcs.txt \
    testmocchan_com \
    teststcom \
    teststc_eg1 teststc_eg2 teststc_eg3 teststc_eg4 teststc_eg5 \
    teststc_eg6 teststc_eg7 teststc_eg8 teststc_eg9 teststc_eg10 \
    testxmlchanom

# Data-driven tests: WCS encoding conversion and Mapping simplification.
check_PROGRAMS += wcsconverter simplify
wcsconverter_SOURCES = wcsconverter.f
wcsconverter_LDADD   = $(FTEST_LDADD)
simplify_SOURCES     = simplify.f
simplify_LDADD       = $(FTEST_LDADD)
TESTS += test_wcsconverter.sh test_simplify.sh

# Source and LDADD declarations for each Fortran test.
testchannel_f_SOURCES    = testchannel.f
testchannel_f_LDADD      = $(FTEST_LDADD)

testchebymap_f_SOURCES   = testchebymap.f
testchebymap_f_LDADD     = $(FTEST_LDADD)

testcmpframe_f_SOURCES   = testcmpframe.f
testcmpframe_f_LDADD     = $(FTEST_LDADD)

testcmpmap_f_SOURCES     = testcmpmap.f
testcmpmap_f_LDADD       = $(FTEST_LDADD)

testfitschan_f_SOURCES   = testfitschan.f
testfitschan_f_LDADD     = $(FTEST_LDADD)

testfitstable_f_SOURCES  = testfitstable.f
testfitstable_f_LDADD    = $(FTEST_LDADD)

testflux_f_SOURCES       = testflux.f
testflux_f_LDADD         = $(FTEST_LDADD)

testframeset_f_SOURCES   = testframeset.f
testframeset_f_LDADD     = $(FTEST_LDADD)

testkeymap_f_SOURCES     = testkeymap.f
testkeymap_f_LDADD       = $(FTEST_LDADD)

testlutmap_f_SOURCES     = testlutmap.f
testlutmap_f_LDADD       = $(FTEST_LDADD)

testmapping_f_SOURCES    = testmapping.f
testmapping_f_LDADD      = $(FTEST_LDADD)

testmoc_f_SOURCES        = testmoc.f
testmoc_f_LDADD          = $(FTEST_LDADD)

testmocchan_f_SOURCES    = testmocchan.f
testmocchan_f_LDADD      = $(FTEST_LDADD)

testnormmap_f_SOURCES    = testnormmap.f
testnormmap_f_LDADD      = $(FTEST_LDADD)

testpolymap_f_SOURCES    = testpolymap.f
testpolymap_f_LDADD      = $(FTEST_LDADD)

testrate_f_SOURCES       = testrate.f
testrate_f_LDADD         = $(FTEST_LDADD)

testratemap_f_SOURCES    = testratemap.f
testratemap_f_LDADD      = $(FTEST_LDADD)

testrebin_f_SOURCES      = testrebin.f
testrebin_f_LDADD        = $(FTEST_LDADD)

testrebinseq_f_SOURCES   = testrebinseq.f
testrebinseq_f_LDADD     = $(FTEST_LDADD)

testregions_f_SOURCES    = testregions.f
testregions_f_LDADD      = $(FTEST_LDADD)

testskyframe_f_SOURCES   = testskyframe.f
testskyframe_f_LDADD     = $(FTEST_LDADD)

testspecflux_f_SOURCES   = testspecflux.f
testspecflux_f_LDADD     = $(FTEST_LDADD)

testspecframe_f_SOURCES  = testspecframe.f
testspecframe_f_LDADD    = $(FTEST_LDADD)

testsplinemap_f_SOURCES  = testsplinemap.f
testsplinemap_f_LDADD    = $(FTEST_LDADD)

teststc_f_SOURCES        = teststc.f
teststc_f_LDADD          = $(FTEST_LDADD)

teststcschan_f_SOURCES   = teststcschan.f
teststcschan_f_LDADD     = $(FTEST_LDADD)

testswitchmap_f_SOURCES  = testswitchmap.f
testswitchmap_f_LDADD    = $(FTEST_LDADD)

testtable_f_SOURCES      = testtable.f
testtable_f_LDADD        = $(FTEST_LDADD)

testtime_f_SOURCES       = testtime.f
testtime_f_LDADD         = $(FTEST_LDADD)

testtrangrid_f_SOURCES   = testtrangrid.f
testtrangrid_f_LDADD     = $(FTEST_LDADD)

testunitnormmap_f_SOURCES = testunitnormmap.f
testunitnormmap_f_LDADD  = $(FTEST_LDADD)

testxmlchan_f_SOURCES    = testxmlchan.f
testxmlchan_f_LDADD      = $(FTEST_LDADD)

testzoommap_f_SOURCES    = testzoommap.f
testzoommap_f_LDADD      = $(FTEST_LDADD)

endif  # NOFORTRAN
