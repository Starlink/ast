dnl ============================================================
dnl Starlink fallback stubs (used when starlink macros not present)
dnl ============================================================
dnl
dnl This file is ONLY used in standalone (no-Starlink) builds.
dnl It is copied into m4/ by bootstrap.local before autoreconf.
dnl
dnl In Starlink builds, starconf.m4 from the Starlink buildsupport
dnl aclocal directory is used directly (the Starlink aclocal finds it
dnl automatically via its default search path), and this file is not
dnl present in m4/.
dnl ============================================================

m4_ifndef([STAR_DEFAULTS],
[AC_DEFUN([STAR_DEFAULTS], [
dnl Provide the minimal substitution variables needed for Makefile.am
dnl and the generated files (component.xml.in, etc.) when building
dnl without the real Starlink build infrastructure (starconf).

      AS_VAR_SET_IF([stardocsdir], [],
        [AS_VAR_SET([stardocsdir], [${datadir}/stardocs])])
      AC_SUBST([stardocsdir])

      AS_VAR_SET_IF([starnewsdir], [],
        [AS_VAR_SET([starnewsdir], [${datadir}/starnews])])
      AC_SUBST([starnewsdir])

      dnl Fortran include flags: empty for standalone builds
      AC_SUBST([STAR_FFLAGS], [""])

      dnl MANIFEST is a Starlink install tracker; 'false' makes it a no-op
      AC_SUBST([MANIFEST], [false])

      dnl STARLINK installation root: empty for standalone builds
      AC_SUBST([STARLINK], [""])

      dnl component.xml substitution variables
      AC_SUBST([STAR_DEPENDENCIES_ATTRIBUTES], [""])
      AC_SUBST([STAR_DEPENDENCIES_CHILDREN], [""])
      AC_SUBST([STAR_DOCUMENTATION], [""])

      dnl Predist status initialization, for STAR_PREDIST_SOURCES
      _star_predist_status=unknown
      PREDIST='#'  # safe default
      AC_SUBST(PREDIST)
])])

m4_ifndef([STAR_DECLARE_DEPENDENCIES],
[AC_DEFUN([STAR_DECLARE_DEPENDENCIES], [dnl
dnl no-op fallback
])])

m4_ifndef([STAR_LATEX_DOCUMENTATION],
[AC_DEFUN([STAR_LATEX_DOCUMENTATION], [dnl
dnl No documentation built in standalone mode.
  AC_SUBST([STAR_LATEX_DOCUMENTATION], [""])
])])

m4_ifndef([STAR_PREDIST_SOURCES],
[AC_DEFUN([STAR_PREDIST_SOURCES], [dnl
m4_ifval([$1], [], [AC_FATAL([$0: called with no stamp file])])dnl
_star_tmp='$1'
for marker in $_star_tmp
do
    if test -f "$srcdir/$marker"; then
        _star_predist_marker_present=:
        AC_MSG_NOTICE([found predist marker file $marker])
    else
        _star_predist_marker_present=false
    fi
    case $_star_predist_status in
        unknown)
            if $_star_predist_marker_present; then
                # we do want to build sourceset files
                _star_predist_status=predist
                PREDIST=
                AC_MSG_NOTICE([in pre-distribution state])
            else
                _star_predist_status=postdist
                PREDIST='#'
                AC_MSG_NOTICE([in post-distribution state])
            fi
            ;;
        predist)
            if $_star_predist_marker_present; then
                : OK
            else
                AC_MSG_WARN([Building predist, but marker file $marker is not present])
            fi
            ;;
        postdist)
            if $_star_predist_marker_present; then
                AC_MSG_WARN([In postdistribution state, but predist marker file $marker is present])
            fi
            ;;
        *)
            AC_MSG_ERROR([impossible predist status $_star_predist_status])
            ;;
    esac
done
])])

m4_ifndef([STAR_CHECK_PROGS],
[AC_DEFUN([STAR_CHECK_PROGS], [dnl
AC_CHECK_PROGS([$1], [$2], [:])
])])

m4_ifndef([STAR_CNF_F2C_COMPATIBLE],
[AC_DEFUN([STAR_CNF_F2C_COMPATIBLE], [dnl
dnl Determine whether the Fortran compiler returns REAL function results
dnl as double (f2c-compatible) or float.  gfortran is f2c-compatible.
dnl Sets REAL_FUNCTION_TYPE to 'double' (f2c) or 'float'.
  if test "$ac_cv_fc_compiler_gnu" = yes; then
    AC_SUBST([REAL_FUNCTION_TYPE], [double])
  else
    AC_SUBST([REAL_FUNCTION_TYPE], [float])
  fi
])])

m4_ifndef([STAR_CNF_TRAIL_TYPE],
[AC_DEFUN([STAR_CNF_TRAIL_TYPE], [dnl
dnl Determine the C type used for hidden Fortran character string length
dnl arguments (the "TRAIL" type).  gfortran >= 8.1 uses size_t; older
dnl versions use int.  This substitutes TRAIL_TYPE into src/f77.h.
  AC_CACHE_CHECK([type used for Fortran string lengths],
    [star_cv_cnf_trail_type],
    [star_cv_cnf_trail_type=int
     if test "$ac_cv_fc_compiler_gnu" = yes; then
       star_gfortran_major=$("$FC" -dumpversion 2>/dev/null | sed 's/\..*//')
       if test -n "$star_gfortran_major" && test "$star_gfortran_major" -ge 8 2>/dev/null; then
         star_cv_cnf_trail_type=size_t
       fi
     fi])
  TRAIL_TYPE=$star_cv_cnf_trail_type
  AC_SUBST([TRAIL_TYPE])
])])

m4_ifndef([STAR_CNF_BLANK_COMMON],
[AC_DEFUN([STAR_CNF_BLANK_COMMON], [dnl
dnl Determine the symbol for the Fortran blank common block.
dnl gfortran >= 4.x uses __BLNK__; older/other compilers use _BLNK__.
  if test "$ac_cv_fc_compiler_gnu" = yes; then
    AC_SUBST([BLANK_COMMON_SYMBOL], [__BLNK__])
  else
    AC_SUBST([BLANK_COMMON_SYMBOL], [_BLNK__])
  fi
])])

m4_ifndef([STAR_MESSGEN],
[AC_DEFUN([STAR_MESSGEN], [dnl
dnl no-op fallback
])])
