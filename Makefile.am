## Process this file with automake to produce Makefile.in

# First declare various groups of files.  These were initially extracted
# from the grp.make file, as constructed by the SDT newdev command
GRP_C_ROUTINES = \
    src/yamlchan.c \
    src/axis.c \
    src/box.c \
    src/channel.c \
    src/chebymap.c \
    src/circle.c \
    src/cmpframe.c \
    src/cmpmap.c \
    src/cmpregion.c \
    src/dsbspecframe.c \
    src/dssmap.c \
    src/ellipse.c \
    src/error.c \
    src/fitschan.c \
    src/fitstable.c \
    src/fluxframe.c \
    src/frame.c \
    src/frameset.c \
    src/globals.c \
    src/grismmap.c \
    src/interval.c \
    src/intramap.c \
    src/keymap.c \
    src/loader.c \
    src/lutmap.c \
    src/mapping.c \
    src/mathmap.c \
    src/matrixmap.c \
    src/memory.c \
    src/moc.c \
    src/mocchan.c \
    src/normmap.c \
    src/nullregion.c \
    src/object.c \
    src/pcdmap.c \
    src/permmap.c \
    src/plot.c \
    src/plot3d.c \
    src/pointlist.c \
    src/pointset.c \
    src/polygon.c \
    src/polymap.c \
    src/prism.c \
    src/ratemap.c \
    src/region.c \
    src/selectormap.c \
    src/shiftmap.c \
    src/skyaxis.c \
    src/skyframe.c \
    src/slamap.c \
    src/specfluxframe.c \
    src/specframe.c \
    src/specmap.c \
    src/sphmap.c \
    src/stc.c \
    src/stccatalogentrylocation.c \
    src/stcobsdatalocation.c \
    src/stcresourceprofile.c \
    src/stcschan.c \
    src/stcsearchlocation.c \
    src/switchmap.c \
    src/table.c \
    src/timeframe.c \
    src/timemap.c \
    src/tranmap.c \
    src/unit.c \
    src/unitmap.c \
    src/unitnormmap.c \
    src/wcsmap.c \
    src/winmap.c \
    src/xml.c \
    src/xmlchan.c \
    src/xphmap.c \
    src/zoommap.c


# The C source files required for the Fortran interface
if !NOFORTRAN
F_C_ROUTINES = \
    src/fyamlchan.c \
    src/c2f77.c \
    src/fbox.c \
    src/fchannel.c \
    src/fchebymap.c \
    src/fcircle.c \
    src/fcmpframe.c \
    src/fcmpmap.c \
    src/fcmpregion.c \
    src/fdsbspecframe.c \
    src/fdssmap.c \
    src/fellipse.c \
    src/ferror.c \
    src/ffitschan.c \
    src/ffitstable.c \
    src/ffluxframe.c \
    src/fframe.c \
    src/fframeset.c \
    src/fgrismmap.c \
    src/finterval.c \
    src/fintramap.c \
    src/fkeymap.c \
    src/flutmap.c \
    src/fmapping.c \
    src/fmathmap.c \
    src/fmatrixmap.c \
    src/fmoc.c \
    src/fmocchan.c \
    src/fnormmap.c \
    src/fnullregion.c \
    src/fobject.c \
    src/fpcdmap.c \
    src/fpermmap.c \
    src/fplot.c \
    src/fplot3d.c \
    src/fpointlist.c \
    src/fpolygon.c \
    src/fpolymap.c \
    src/fprism.c \
    src/fratemap.c \
    src/fregion.c \
    src/fselectormap.c \
    src/fshiftmap.c \
    src/fskyframe.c \
    src/fslamap.c \
    src/fspecfluxframe.c \
    src/fspecframe.c \
    src/fspecmap.c \
    src/fsphmap.c \
    src/fstc.c \
    src/fstccatalogentrylocation.c \
    src/fstcobsdatalocation.c \
    src/fstcresourceprofile.c \
    src/fstcschan.c \
    src/fstcsearchlocation.c \
    src/fswitchmap.c \
    src/ftable.c \
    src/ftimeframe.c \
    src/ftimemap.c \
    src/ftranmap.c \
    src/funitmap.c \
    src/funitnormmap.c \
    src/fwcsmap.c \
    src/fwinmap.c \
    src/fxmlchan.c \
    src/fzoommap.c
else
F_C_ROUTINES =
endif

# Header files which contribute to the "ast.h" file, organised to correspond
# with the class hierarchy.
AST_H_FILES = \
    src/xml.h \
    wcslib/wcstrig.h \
    wcslib/proj.h \
    src/memory.h \
    src/error.h \
    src/globals.h \
    src/unit.h \
    ast_err.h \
    src/version.h \
    src/object.h \
       src/keymap.h \
          src/table.h \
             src/fitstable.h \
       src/pointset.h \
       src/axis.h \
          src/skyaxis.h \
       src/mapping.h \
          src/cmpmap.h \
          src/dssmap.h \
          src/grismmap.h \
          src/intramap.h \
          src/lutmap.h \
          src/mathmap.h \
          src/matrixmap.h \
          src/pcdmap.h \
          src/permmap.h \
          src/polymap.h \
          src/chebymap.h \
          src/ratemap.h \
          src/normmap.h \
          src/shiftmap.h \
          src/slamap.h \
          src/specmap.h \
          src/sphmap.h \
          src/timemap.h \
          src/selectormap.h \
          src/switchmap.h \
          src/tranmap.h \
          src/unitmap.h \
          src/unitnormmap.h \
          src/wcsmap.h \
          src/winmap.h \
          src/xphmap.h \
          src/zoommap.h \
          src/frame.h \
             src/cmpframe.h \
                 src/specfluxframe.h \
             src/fluxframe.h \
             src/frameset.h \
                 src/plot.h \
                     src/plot3d.h \
             src/skyframe.h \
             src/specframe.h \
                 src/dsbspecframe.h \
             src/region.h \
                src/box.h \
                src/circle.h \
                src/cmpregion.h \
                src/ellipse.h \
                src/interval.h \
                src/moc.h \
                src/nullregion.h \
                src/pointlist.h \
                src/polygon.h \
                src/prism.h \
                src/stc.h \
                   src/stcresourceprofile.h \
                   src/stcsearchlocation.h \
                   src/stccatalogentrylocation.h \
                   src/stcobsdatalocation.h \
             src/timeframe.h \
       src/channel.h \
          src/yamlchan.h \
          src/fitschan.h \
          src/mocchan.h \
          src/stcschan.h \
          src/xmlchan.h

# All the (C) include files required to build the library.
GRP_C_INCLUDE_FILES = \
    $(AST_H_FILES) \
    pal2ast.h \
    erfa2ast.h \
    src/ems.h \
    src/err.h \
    src/Ers.h \
    src/f77.h \
    src/grf.h \
    src/grf3d.h \
    src/pg3d.h \
    src/loader.h \
    src/skyaxis.h \
    src/stc.h \
    src/stcresourceprofile.h \
    src/stcsearchlocation.h \
    src/stccatalogentrylocation.h \
    src/stcobsdatalocation.h \
    wcslib/wcsmath.h \
    wcslib/wcstrig.h \
    src/xmlchan.h

if !NOFORTRAN
F_C_INCLUDE_FILES = \
    src/c2f77.h

# The following list should include AST_PAR, but that must not be
# distributed, and so it is listed separately in
# nodist_libast_la_SOURCES below.
GRP_F_INCLUDE_FILES = \
    GRF_PAR \
    AST_ERR

else
F_C_INCLUDE_FILES =
GRP_F_INCLUDE_FILES =
endif

# If we have no Fortran we are not building f77.h and we probably
# do not want a Fortran runtime. This requires that PGPLOT is disabled.
# We replace it with the stub GRF interface. An alternative would
# be to have the PGPLOT wrappers use a preprocessor symbol to build
# the pgplot to always error if used.
if !NOFORTRAN
GRF_PGPLOT_SOURCES = \
    src/grf_pgplot.c
GRF3D_PGPLOT_SOURCES = \
    src/grf3d_pgplot.c
else
GRF_PGPLOT_SOURCES = \
    src/grf_5.6.c
GRF3D_PGPLOT_SOURCES = \
    src/grf3d.c
endif


## Following declaration isn't used
## LATEX_DOCUMENTATION_FILES = \
##     sun210.tex \
##     sun211.tex

DOCUMENTATION_PRODUCTS = $(PAPER_DOCUMENTATION) $(HYPER_DOCUMENTATION)
PAPER_DOCUMENTATION = sun210.tex sun211.tex sun210.pdf sun211.pdf
HYPER_DOCUMENTATION = sun210.htx_tar sun211.htx_tar

PDF_FIGURES = \
    cmpframe.pdf \
    complex.pdf \
    frames.pdf \
    frameset.pdf \
    fronta.pdf \
    fronta_bw.pdf \
    frontb.pdf \
    frontb_bw.pdf \
    frontc.pdf \
    frontc_bw.pdf \
    fsalign.pdf \
    fsconvert.pdf \
    fsexample.pdf \
    fsmerge.pdf \
    fsremap.pdf \
    gridplot.pdf \
    gridplot_bw.pdf \
    mapping.pdf \
    overgrid.pdf \
    overgrid_bw.pdf \
    parallel.pdf \
    series.pdf \
    simpexamp.pdf

WCSLIB_FILES = \
    wcslib/proj.c \
    wcslib/tpn.c \
    wcslib/proj.h \
    wcslib/wcstrig.c \
    wcslib/wcsmath.h \
    wcslib/wcstrig.h

STAR_PAL_FILES = \
    pal/pal.h \
    pal/palAddet.c \
    pal/palAmpqk.c \
    pal/palCaldj.c \
    pal/palDat.c \
    pal/palDe2h.c \
    pal/palDeuler.c \
    pal/palDh2e.c \
    pal/palDjcal.c \
    pal/palDmat.c \
    pal/palDrange.c \
    pal/palDs2tp.c \
    pal/palDtp2s.c \
    pal/palDtps2c.c \
    pal/palDtt.c \
    pal/palEcmat.c \
    pal/palEqgal.c \
    pal/palEtrms.c \
    pal/palEvp.c \
    pal/palFk45z.c \
    pal/palFk524.c \
    pal/palFk54z.c \
    pal/palGaleq.c \
    pal/palGalsup.c \
    pal/palMappa.c \
    pal/palMapqkz.c \
    pal/palOne2One.c \
    pal/palPrebn.c \
    pal/palPrec.c \
    pal/palPrenut.c \
    pal/palPvobs.c \
    pal/palRvgalc.c \
    pal/palRvlg.c \
    pal/palRvlsrd.c \
    pal/palRvlsrk.c \
    pal/palSubet.c \
    pal/palSupgal.c \
    pal/pal1sofa.h \
    pal/palmac.h

ERFA_FILES = \
    erfa/00READ.ME \
    erfa/erfa.h \
    erfa/erfam.h \
    erfa/a2af.c \
    erfa/a2tf.c \
    erfa/af2a.c \
    erfa/anp.c \
    erfa/anpm.c \
    erfa/bi00.c \
    erfa/bp00.c \
    erfa/bp06.c \
    erfa/bpn2xy.c \
    erfa/c2i00a.c \
    erfa/c2i00b.c \
    erfa/c2i06a.c \
    erfa/c2ibpn.c \
    erfa/c2ixy.c \
    erfa/c2ixys.c \
    erfa/c2s.c \
    erfa/c2t00a.c \
    erfa/c2t00b.c \
    erfa/c2t06a.c \
    erfa/c2tcio.c \
    erfa/c2teqx.c \
    erfa/c2tpe.c \
    erfa/c2txy.c \
    erfa/cal2jd.c \
    erfa/cp.c \
    erfa/cpv.c \
    erfa/cr.c \
    erfa/d2dtf.c \
    erfa/d2tf.c \
    erfa/dat.c \
    erfa/dtdb.c \
    erfa/dtf2d.c \
    erfa/ee00.c \
    erfa/ee00a.c \
    erfa/ee00b.c \
    erfa/ee06a.c \
    erfa/eect00.c \
    erfa/eform.c \
    erfa/eo06a.c \
    erfa/eors.c \
    erfa/epb.c \
    erfa/epb2jd.c \
    erfa/epj.c \
    erfa/epj2jd.c \
    erfa/epv00.c \
    erfa/eqeq94.c \
    erfa/era00.c \
    erfa/fad03.c \
    erfa/fae03.c \
    erfa/faf03.c \
    erfa/faju03.c \
    erfa/fal03.c \
    erfa/falp03.c \
    erfa/fama03.c \
    erfa/fame03.c \
    erfa/fane03.c \
    erfa/faom03.c \
    erfa/fapa03.c \
    erfa/fasa03.c \
    erfa/faur03.c \
    erfa/fave03.c \
    erfa/fk52h.c \
    erfa/fk5hip.c \
    erfa/fk5hz.c \
    erfa/fw2m.c \
    erfa/fw2xy.c \
    erfa/gc2gd.c \
    erfa/gc2gde.c \
    erfa/gd2gc.c \
    erfa/gd2gce.c \
    erfa/gmst00.c \
    erfa/gmst06.c \
    erfa/gmst82.c \
    erfa/gst00a.c \
    erfa/gst00b.c \
    erfa/gst06.c \
    erfa/gst06a.c \
    erfa/gst94.c \
    erfa/h2fk5.c \
    erfa/hfk5z.c \
    erfa/ir.c \
    erfa/jd2cal.c \
    erfa/jdcalf.c \
    erfa/num00a.c \
    erfa/num00b.c \
    erfa/num06a.c \
    erfa/numat.c \
    erfa/nut00a.c \
    erfa/nut00b.c \
    erfa/nut06a.c \
    erfa/nut80.c \
    erfa/nutm80.c \
    erfa/obl06.c \
    erfa/obl80.c \
    erfa/p06e.c \
    erfa/p2pv.c \
    erfa/p2s.c \
    erfa/pap.c \
    erfa/pas.c \
    erfa/pb06.c \
    erfa/pdp.c \
    erfa/pfw06.c \
    erfa/plan94.c \
    erfa/pm.c \
    erfa/pmat00.c \
    erfa/pmat06.c \
    erfa/pmat76.c \
    erfa/pmp.c \
    erfa/pn.c \
    erfa/pn00.c \
    erfa/pn00a.c \
    erfa/pn00b.c \
    erfa/pn06.c \
    erfa/pn06a.c \
    erfa/pnm00a.c \
    erfa/pnm00b.c \
    erfa/pnm06a.c \
    erfa/pnm80.c \
    erfa/pom00.c \
    erfa/ppp.c \
    erfa/ppsp.c \
    erfa/pr00.c \
    erfa/prec76.c \
    erfa/pv2p.c \
    erfa/pv2s.c \
    erfa/pvdpv.c \
    erfa/pvm.c \
    erfa/pvmpv.c \
    erfa/pvppv.c \
    erfa/pvstar.c \
    erfa/pvu.c \
    erfa/pvup.c \
    erfa/pvxpv.c \
    erfa/pxp.c \
    erfa/refco.c \
    erfa/rm2v.c \
    erfa/rv2m.c \
    erfa/rx.c \
    erfa/rxp.c \
    erfa/rxpv.c \
    erfa/rxr.c \
    erfa/ry.c \
    erfa/rz.c \
    erfa/s00.c \
    erfa/s00a.c \
    erfa/s00b.c \
    erfa/s06.c \
    erfa/s06a.c \
    erfa/s2c.c \
    erfa/s2p.c \
    erfa/s2pv.c \
    erfa/s2xpv.c \
    erfa/sepp.c \
    erfa/seps.c \
    erfa/sp00.c \
    erfa/starpm.c \
    erfa/starpv.c \
    erfa/sxp.c \
    erfa/sxpv.c \
    erfa/taitt.c \
    erfa/taiut1.c \
    erfa/taiutc.c \
    erfa/tcbtdb.c \
    erfa/tcgtt.c \
    erfa/tdbtcb.c \
    erfa/tdbtt.c \
    erfa/tf2a.c \
    erfa/tf2d.c \
    erfa/tr.c \
    erfa/trxp.c \
    erfa/trxpv.c \
    erfa/tttai.c \
    erfa/tttcg.c \
    erfa/tttdb.c \
    erfa/ttut1.c \
    erfa/ut1tai.c \
    erfa/ut1tt.c \
    erfa/ut1utc.c \
    erfa/utctai.c \
    erfa/utcut1.c \
    erfa/xy06.c \
    erfa/xys00a.c \
    erfa/xys00b.c \
    erfa/xys06a.c \
    erfa/zp.c \
    erfa/zpv.c \
    erfa/zr.c

PAL_FILES = \
    palwrap.c \
    pal.h \
    erfa.h \
    erfam.h

CMINPACK_FILES = \
    cminpack/cminpack.h \
    cminpack/cminpackP.h \
    cminpack/lmder1.c \
    cminpack/lmder.c \
    cminpack/dpmpar.c \
    cminpack/enorm.c \
    cminpack/qrfac.c \
    cminpack/lmpar.c \
    cminpack/qrsolv.c

bin_SCRIPTS = ast_link
dist_bin_SCRIPTS = ast_link_adam
noinst_SCRIPTS = ast_cpp
dist_noinst_SCRIPTS = $(srcdir)/makeh
# Scripts are not distributed by default (since they might be derived objects)
# Add these to the distribution below.  In fact, it would be useful
# and straightforward  to make ast_link{,_adam} derived, since they
# could then have installation directories painlessly edited in to
# them.  This might be a requirement for scripts which supported
# linking against shared libraries.

# Headers required by library users.  Both of the following lines
# indicate headers which are installed.
cincludedir = $(includedir)/star
if !NOFORTRAN
cinclude_HEADERS = GRF_PAR src/grf.h src/grf3d.h
else
cinclude_HEADERS = src/grf.h src/grf3d.h
endif

# Following are generated, so should not be distributed.
if !NOFORTRAN
include_MESSAGES = AST_ERR ast_err.h
nodist_cinclude_HEADERS = ast.h AST_PAR
else
include_MESSAGES = ast_err.h
nodist_cinclude_HEADERS = ast.h
endif

# Also install headers in main include directory for the benefit of
# legacy code that does '#include "ast.h"' instead of '#include
# "star/ast.h"'. Note, currently _MESSAGES files must be installed in
# the main include directory.
if TOPINCLUDE
if !NOFORTRAN
include_HEADERS = GRF_PAR src/grf.h src/grf3d.h
nodist_include_HEADERS = ast.h AST_PAR
else
include_HEADERS = src/grf.h src/grf3d.h
nodist_include_HEADERS = ast.h
endif
endif

if EXTERNAL_PAL
PAL_LIB =
else
PAL_LIB = libast_pal.la
endif

if EXTERNAL_CMINPACK
CMINPACK_LIB =
else
CMINPACK_LIB = libast_cminpack.la
endif

lib_LTLIBRARIES = \
    $(PAL_LIB) \
    $(CMINPACK_LIB) \
    libast.la \
    libast_err.la \
    libast_ems.la \
    libast_drama.la \
    libast_grf3d.la \
    libast_grf_2.0.la \
    libast_grf_3.2.la \
    libast_grf_5.6.la \
    libast_pgplot.la \
    libast_pgplot3d.la

stardocs_DATA = @STAR_LATEX_DOCUMENTATION@
dist_starnews_DATA = ast.news
dist_pkgdata_DATA = COPYING COPYING.LESSER COPYING.LIB

if !NOYAML
yaml = -DYAML
else
yaml =
endif

if !NOTHREADS
threadsafe = -DTHREAD_SAFE
else
threadsafe =
endif

# Make all library code position independent by default. This is handy for
# creating shareable libraries from the static ones (Java JNI libraries).
# Note we do not simply set the AM_CFLAGS variable, as this would also
# apply to programs compiled without using libtool, possibly causing the
# compilation to fail.
if !NOPIC
libast_la_CFLAGS = $(AM_CFLAGS) -prefer-pic $(threadsafe) $(yaml)
else
libast_la_CFLAGS = $(AM_CFLAGS) $(threadsafe) $(yaml)
endif

if !NOPIC
libast_err_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_ems_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_drama_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_grf3d_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_grf_2_0_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_grf_3_2_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_grf_5_6_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_pgplot_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_pgplot3d_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_pal_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
libast_cminpack_la_CFLAGS = $(AM_CFLAGS) -prefer-pic
endif

#  The module containing the main AST classes
libast_la_SOURCES = \
    $(GRP_C_ROUTINES) \
    $(F_C_ROUTINES) \
    $(GRP_C_INCLUDE_FILES) \
    $(F_C_INCLUDE_FILES) \
    $(GRP_F_INCLUDE_FILES) \
    $(WCSLIB_FILES) \
    ast_err.h
libast_la_LDFLAGS = -version-info @version_info@

# Ensure libast links against libraries containing functions used within
# libast. If AST is configured --with-external-pal, then the internal
# libast_pal library will be empty, and we link to an external PAL
# library instead. Do the same for cminpack
if EXTERNAL_PAL
if EXTERNAL_CMINPACK
libast_la_LIBADD =  -lpal -lcminpack
else
libast_la_LIBADD = -lpal libast_cminpack.la
endif
else
if EXTERNAL_CMINPACK
libast_la_LIBADD =  libast_pal.la -lcminpack
else
libast_la_LIBADD = libast_pal.la libast_cminpack.la
endif
endif



# AST_PAR is really part of GRP_F_INCLUDE_FILES, but it must not be
# distributed, so list it separately.
if !NOFORTRAN
nodist_libast_la_SOURCES = \
    ast.h \
    AST_PAR
else
nodist_libast_la_SOURCES = ast.h
endif

#  The default error reporting module.
libast_err_la_SOURCES = src/err_null.c
libast_err_la_LDFLAGS = -version-info @version_info@

#  The error reporting module that uses EMS to deliver errors.
libast_ems_la_SOURCES = src/err_ems.c
libast_ems_la_LDFLAGS = -version-info @version_info@

#  The error reporting module that uses DRAMA Ers to deliver errors.
libast_drama_la_SOURCES = src/err_drama.c
libast_drama_la_LDFLAGS = -version-info @version_info@

#  The module containing null implementations of the 3D graphics routines
#  required by AST
libast_grf3d_la_SOURCES = src/grf3d.c
libast_grf3d_la_LDFLAGS = -version-info @version_info@

#  The module containing null implementations of the graphics routines
#  required by AST V2.0
libast_grf_2_0_la_SOURCES = src/grf_2.0.c
libast_grf_2_0_la_LDFLAGS = -version-info @version_info@

#  The module containing null implementations of the graphics routines
#  added by AST V3.2 and not present in V2.0
libast_grf_3_2_la_SOURCES = src/grf_3.2.c
libast_grf_3_2_la_LDFLAGS = -version-info @version_info@

#  The module containing null implementations of the graphics routines
#  added by AST V5.6 and not present in V3.2
libast_grf_5_6_la_SOURCES = src/grf_5.6.c
libast_grf_5_6_la_LDFLAGS = -version-info @version_info@

#  The graphics module that uses PGPLOT for 2D graphical output.
libast_pgplot_la_SOURCES = $(GRF_PGPLOT_SOURCES)
libast_pgplot_la_LDFLAGS = -version-info @version_info@

#  The graphics module that uses PGPLOT for 3D graphical output.
libast_pgplot3d_la_SOURCES = $(GRF3D_PGPLOT_SOURCES)
libast_pgplot3d_la_LDFLAGS = -version-info @version_info@

#  Positional astronomy libraries.
libast_pal_la_SOURCES = $(PAL_FILES)
libast_pal_la_LDFLAGS = -version-info @version_info@

#  Cminpack
libast_cminpack_la_SOURCES = $(CMINPACK_FILES)
libast_cminpack_la_LDFLAGS = -version-info @version_info@

# The following files are built by the targets in this makefile.
MAINTAINERCLEANFILES = src/version.h builddocs addversion \
	ast.h $(DOCUMENTATION_PRODUCTS)
if !NOFORTRAN
CLEANFILES = AST_PAR ast.h
else
CLEANFILES = ast.h
endif

# Special cases start here

# The AST_PAR include file is produced by compiling the astbad.c
# program and editing its output into the ast_par.source file (while
# also changing the "E" exponent character to a "D").  The astbad.c
# and ast_par.source must be distributed (the generation of the
# AST__BAD value must be done on the build host) but not installed.
if !NOFORTRAN
noinst_PROGRAMS = astbad$(EXEEXT)
astbad_SOURCES = src/astbad.c src/pointset.h
AST_PAR: $(srcdir)/ast_par.source astbad$(EXEEXT)
	ASTBAD=`./astbad$(EXEEXT) AST__BAD` && \
	ASTNAN=`./astbad$(EXEEXT) AST__NAN` && \
	ASTNANF=`./astbad$(EXEEXT) AST__NANF` && \
	sed -e 's/<AST__BAD>/'`echo $$ASTBAD | tr 'E' 'D'`'/' \
	    -e 's/<AST__NAN>/'`echo $$ASTNAN | tr 'E' 'D'`'/' \
	    -e 's/<AST__NANR>/'`echo $$ASTNANF | tr 'E' 'D'`'/' \
	    $(srcdir)/ast_par.source >$@
endif

# ast_link is generated from ast_link.in; ast_link_adam does not
# need configuration, and so is not mentioned in AC_CONFIG_FILES within
# configure.ac, and so is not distributed automatically.
#
# makeh is required in order to build ast.h after distribution (see below).
EXTRA_DIST = $(srcdir)/ast_par.source sun210_figures sun211_figures pal erfa cminpack

# ast.h depends on the error-facility files.  ast.h _could_ be made
# before distribution, but since it's generated from other distributed
# files, it's more robust to distribute makeh and make ast.h on the
# build host.
ast.h: $(AST_H_FILES) ast_err.h $(srcdir)/makeh config.h
	@PERL@ $(srcdir)/makeh -s $(srcdir) $(AST_H_FILES) >$@

# AST_ERR and ast_err.h are `generated source files', and so must be
# generated before any other compilations are done.  Note that these
# files are generated on the distribution host, and so this
# declaration has no effect post-distribution.
#
# AST_PAR is also a generated source file, but it should _not_ be
# included in the list of BUILT_SOURCES, otherwise `make' tries to make
# it before it makes the `astbad' program it depends on.  Instead,
# just rely on the dependencies expressed in the main body above to
# have AST_PAR built before it is needed.
#
# version.h is included by object.h, and thus indirectly by most modules.
# It's most straightforward to build it at the beginning.
if !NOFORTRAN
BUILT_SOURCES = AST_ERR ast_err.h src/version.h
else
BUILT_SOURCES = ast_err.h src/version.h
endif

#  Form a second link to the main object library (static and shared). This
#  is used when a second pass through the library is needed during linking
#  to resolve references made within other AST libraries (e.g. the grf
#  modules, etc), to functions defined within libast (e.g. memory management
#  and error reporting functions). Do not forget to change the contents of
#  the libast_pass2.la file to refer to libast_pass2.* rather than libast.*.
#  Use target install-exec-hook rather than install-exec-local so that the
#  soft links and files are created *after* the main library has been
#  installed.
install-exec-hook: libast.la
	$(mkdir_p) $(DESTDIR)$(libdir)
	cd $(DESTDIR)$(libdir); \
	for f in `ls libast.*`; do \
	   ff=`echo $$f | sed -e 's/libast/libast_pass2/'`; \
	   if test -f "$$ff"; then rm "$$ff"; fi; \
	   $(LN_S) $$f $$ff; \
	   $(MANIFEST) && echo "MANIFEST:$(DESTDIR)$(libdir)/$$ff" || :; \
	done; \
	if test -f "libast.la"; then \
	   if test -f "libast_pass2.la"; then rm "libast_pass2.la"; fi; \
	   sed -e 's/libast\./libast_pass2\./g' libast.la > libast_pass2.la; \
	fi

# Make pre-distribution files.  These are files which are required for
# building the distribution, but are not themselves distributed.
# The source files here should be mentioned in STAR_PREDIST_SOURCES in
# configure.ac
@PREDIST@predist_subs = sed \
@PREDIST@    -e 's,@PACKAGE_VERSION\@,$(PACKAGE_VERSION),' \
@PREDIST@    -e 's,@PACKAGE_VERSION_MAJOR\@,$(PACKAGE_VERSION_MAJOR),' \
@PREDIST@    -e 's,@PACKAGE_VERSION_MINOR\@,$(PACKAGE_VERSION_MINOR),' \
@PREDIST@    -e 's,@PACKAGE_VERSION_RELEASE\@,$(PACKAGE_VERSION_RELEASE),' \
@PREDIST@    -e 's,@PERL\@,$(PERL),' \
@PREDIST@    -e 's,@STARLINK\@,$(STARLINK),'

@PREDIST@src/version.h: src/version.h.in configure.ac
@PREDIST@	rm -f $@; $(predist_subs) src/version.h.in >$@
@PREDIST@builddocs: builddocs.in configure.ac
@PREDIST@	rm -f $@; $(predist_subs) builddocs.in >$@; chmod +x $@
@PREDIST@addversion: addversion.in configure.ac
@PREDIST@	rm -f $@; $(predist_subs) addversion.in >$@; chmod +x $@

# Documentation
@PREDIST@$(PAPER_DOCUMENTATION): sun211_figures builddocs addversion
@PREDIST@	./builddocs

# The contents of the sun211_figures directory is identical to that
# sun210_figures
sun211_figures: sun210_figures
	$(LN_S) sun210_figures sun211_figures

# Installation check

TESTS = ast_test
check_PROGRAMS = ast_test
ast_test_SOURCES = ast_test.c

test: install
	cd ast_tester && STARLINK=@STARLINK@ PGPLOT_DIR=@STARLINK@/bin FFLAGS="@STAR_FFLAGS_TEST@" ./ast_tester -nd

#ast_test_LDADD = `ast_link`
# Expand ast_link to avoid libast_pass2, which causes problems for Solaris
ast_test_LDADD = @LIBPAL@ @LIBCMINPACK@ libast.la libast_grf_3.2.la libast_grf_5.6.la libast_grf_2.0.la libast_grf3d.la libast_err.la -lm

#  Need to include latex support files in the distribution tar ball so
#  that the docs can be built from the tex source files. Requires environment
#  variable STARLATEXSUPPORT to be deined. Is there a better way to do this?
dist-hook:
	cp -p $(STARLATEXSUPPORT)/starlink.cls $(distdir)/
	cp -p $(STARLATEXSUPPORT)/starabbrev.sty $(distdir)/
	cp -p $(STARLATEXSUPPORT)/starstyle.sty $(distdir)/
	cp -p $(STARLATEXSUPPORT)/sst.sty $(distdir)/
