# ===========================================================================
#  AST YAML backend detection
# ===========================================================================
#
# SYNOPSIS
#
#   AST_CHECK_YAML
#
# DESCRIPTION
#
#   Detect a YAML backend for AST's YamlChan support--either libyaml or
#   libfyaml--and configure the build to use it.  The backend is chosen
#   with:
#
#       --with-yaml[=BACKEND]
#
#   where BACKEND is one of:
#
#       auto (default)  use libyaml if available, otherwise libfyaml;
#                       building without YAML support is not an error
#       yaml            require libyaml
#       fyaml           require libfyaml
#       yes             same as auto
#       no              disable YAML support
#
#   A library installed outside the compiler's default search paths can be
#   located with the precious variables YAML_CFLAGS / YAML_LIBS (libyaml) or
#   FYAML_CFLAGS / FYAML_LIBS (libfyaml), for example:
#
#       ./configure YAML_CFLAGS=-I/opt/foo/include \
#                   YAML_LIBS='-L/opt/foo/lib -lyaml'
#
#   Setting either variable of a backend (and none of the other's) implies a
#   request for that backend, so a failure to use it is reported as an error
#   rather than silently ignored.
#
#   On success the chosen backend's flags are applied to CPPFLAGS/LIBS so
#   the whole package is built and linked against it, and these outputs are
#   produced:
#
#       Makefile : @AST_YAML_DEFINE@ the backend defines for the compiler
#                                     command line -- -DYAML or -DFYAML (which
#                                     yaml_backend.h switches on) together with
#                                     -DYAML_BACKEND="libyaml"/"libfyaml"; empty
#                                     when no backend was selected
#                  @YAML@ / @FYAML@   "1" or "0"
#                  @YAML_CFLAGS@ @YAML_LIBS@ @FYAML_CFLAGS@ @FYAML_LIBS@
#       automake : NOYAML            true when no backend was selected
#
AC_DEFUN([AST_CHECK_YAML],
[
  AC_ARG_VAR([YAML_CFLAGS],  [C compiler flags for libyaml (e.g. -I/opt/foo/include)])
  AC_ARG_VAR([YAML_LIBS],    [linker flags for libyaml (e.g. -L/opt/foo/lib -lyaml)])
  AC_ARG_VAR([FYAML_CFLAGS], [C compiler flags for libfyaml (e.g. -I/opt/foo/include)])
  AC_ARG_VAR([FYAML_LIBS],   [linker flags for libfyaml (e.g. -L/opt/foo/lib -lfyaml)])

  AC_ARG_WITH([yaml],
    [AS_HELP_STRING([--with-yaml@<:@=BACKEND@:>@],
       [YAML backend: auto (default), yaml or fyaml; --without-yaml disables])],
    [], [with_yaml=auto])

  AS_CASE([$with_yaml],
    [yes|auto], [ast_yaml_backend=auto],
    [yaml],     [ast_yaml_backend=yaml],
    [fyaml],    [ast_yaml_backend=fyaml],
    [no],       [ast_yaml_backend=no],
    [AC_MSG_ERROR([unknown --with-yaml backend '$with_yaml' (use auto, yaml, fyaml or no)])])

  dnl  Setting a backend's location variables (and not the other's) is an
  dnl  implicit request for that backend when none was named explicitly.
  AS_IF([test "x$ast_yaml_backend" = xauto],
    [AS_IF([test -n "$YAML_CFLAGS$YAML_LIBS" && test -z "$FYAML_CFLAGS$FYAML_LIBS"],
       [ast_yaml_backend=yaml],
     [test -n "$FYAML_CFLAGS$FYAML_LIBS" && test -z "$YAML_CFLAGS$YAML_LIBS"],
       [ast_yaml_backend=fyaml])])

  use_yaml=0
  use_fyaml=0

  AS_IF([test "x$ast_yaml_backend" != xno],
  [
    dnl  Try libyaml unless fyaml was specifically requested.
    AS_IF([test "x$ast_yaml_backend" = xauto || test "x$ast_yaml_backend" = xyaml],
      [_AST_YAML_TRY_BACKEND([libyaml], [yaml.h], [yaml_parser_initialize],
                             [-lyaml], [YAML_CFLAGS], [YAML_LIBS], [use_yaml])])
    AS_IF([test "x$use_yaml" = x0 && test "x$ast_yaml_backend" = xyaml],
      [AC_MSG_ERROR([libyaml requested (--with-yaml=yaml) but not usable])])

    dnl  Fall back to (or specifically require) libfyaml.
    AS_IF([test "x$use_yaml" = x0],
      [AS_IF([test "x$ast_yaml_backend" = xauto || test "x$ast_yaml_backend" = xfyaml],
         [_AST_YAML_TRY_BACKEND([libfyaml], [libfyaml.h], [fy_parser_create],
                                [-lfyaml], [FYAML_CFLAGS], [FYAML_LIBS], [use_fyaml])])])
    AS_IF([test "x$use_fyaml" = x0 && test "x$ast_yaml_backend" = xfyaml],
      [AC_MSG_ERROR([libfyaml requested (--with-yaml=fyaml) but not usable])])
  ])

  dnl  yamlchan.c uses YAML_BACKEND before it includes any header, so YAML/FYAML
  dnl  and YAML_BACKEND must all be command-line defines to be visible there.
  AST_YAML_DEFINE=
  AS_IF([test "x$use_yaml" = x1],
    [AST_YAML_DEFINE='-DYAML -DYAML_BACKEND=\"libyaml\"'],
   [test "x$use_fyaml" = x1],
    [AST_YAML_DEFINE='-DFYAML -DYAML_BACKEND=\"libfyaml\"'])

  AC_SUBST([AST_YAML_DEFINE])
  AC_SUBST([YAML],  [$use_yaml])
  AC_SUBST([FYAML], [$use_fyaml])
  AM_CONDITIONAL([NOYAML], [test "x$use_yaml" = x0 && test "x$use_fyaml" = x0])
])

# _AST_YAML_TRY_BACKEND(PRETTY, HEADER, FUNC, DEFAULT_LIBS, CFLAGS_VAR, LIBS_VAR, RESULT_VAR)
# ------------------------------------------------------------------------------------------
# Probe one YAML backend with its CFLAGS/LIBS applied.  Both the header and the
# link must succeed.  On success set RESULT_VAR=1 and leave the flags applied
# to CPPFLAGS/LIBS (so the rest of configure and the whole build use them); on
# failure restore CPPFLAGS/LIBS and clear the backend's location variables.
AC_DEFUN([_AST_YAML_TRY_BACKEND],
[
  AS_IF([test -z "$$6"], [$6=$4])

  ast_yaml_save_CPPFLAGS=$CPPFLAGS
  ast_yaml_save_LIBS=$LIBS
  CPPFLAGS="$CPPFLAGS $$5"
  LIBS="$$6 $LIBS"

  $7=1
  AC_CHECK_HEADER([$2], [], [$7=0])
  AS_IF([test "x$$7" = x1],
    [AC_MSG_CHECKING([for $3 in $1])
     AC_LINK_IFELSE([AC_LANG_CALL([], [$3])],
       [AC_MSG_RESULT([yes])],
       [AC_MSG_RESULT([no]); $7=0])])

  AS_IF([test "x$$7" = x0],
    [CPPFLAGS=$ast_yaml_save_CPPFLAGS
     LIBS=$ast_yaml_save_LIBS
     $5=
     $6=])
])
