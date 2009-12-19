#
# check for lp_solve
#

AC_DEFUN([AC_LPSOLVE],
  [
  AC_SUBST(LPSOLVE_CPPFLAGS)
  AC_SUBST(LPSOLVE_LDFLAGS)
  AC_SUBST(LPSOLVE_LIBS)
    
  AC_ARG_WITH(lp_solve,
    AC_HELP_STRING([--with-lp_solve=DIR],
                   [path to lp_solve library]),
    [
    if test "x$withval" = xno; then
      lpsolve=no
    else
      lpsolve=yes
      if test "x$withval" != xyes; then
        lpsolve_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $lpsolve_dir in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute directory name for --with-lp_solve : $lpsolve_dir]);;
        esac
      fi
    fi
    ]
  )

  ac_cv_have_lpsolve=no

  if test "x$lpsolve" != xno; then
    if test -n "$lpsolve_dir"; then
      if test -d "$lpsolve_dir"; then :; else
        lpsolve_dir=
      fi
    else
      lpsolve_dir=
      for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local/src
      do
	for s in lp_solve lp_solve_5.1 lp_solve_5.0 lp_solve_4.0 lp_solve_3.2
        do
          if test -d "$d/$s"; then
            lpsolve_dir="$d/$s"
            break
          fi
        done
	test -n "$lpsolve_dir" && break
      done
    fi

    if test -n "$lpsolve_dir"; then
      LPSOLVE_CPPFLAGS="-I$lpsolve_dir"
      LPSOLVE_LDFLAGS="-L$lpsolve_dir"
    fi

    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS=$CPPFLAGS
    ac_save_LDFLAGS=$LDFLAGS
    ac_save_LIBS=$LIBS
    CPPFLAGS="$CPPFLAGS $LPSOLVE_CPPFLAGS"
    LDFLAGS="$LDFLAGS $LPSOLVE_LDFLAGS"

    AC_CHECK_HEADER([lpkit.h])
    AC_CHECK_LIB([lpsolve51], [make_lp])
    AC_CHECK_LIB([lpk5], [make_lp])
    AC_CHECK_LIB([lpk], [make_lp])

    CPPFLAGS=$ac_save_CPPFLAGS
    LDFLAGS=$ac_save_LDFLAGS
    LIBS=$ac_save_LIBS
    AC_LANG_RESTORE

  fi

  if test "$ac_cv_header_lpkit_h" = yes && test "$ac_cv_lib_lpsolve51_make_lp" = yes; then
    ac_cv_have_lpsolve=yes
    LPSOLVE_LIBS="-llpsolve51"
  elif test "$ac_cv_header_lpkit_h" = yes && test "$ac_cv_lib_lpk5_make_lp" = yes; then
    ac_cv_have_lpsolve=yes
    LPSOLVE_LIBS="-llpk5"
  elif test "$ac_cv_header_lpkit_h" = yes && test "$ac_cv_lib_lpk_make_lp" = yes; then
    ac_cv_have_lpsolve=yes
    LPSOLVE_LIBS="-llpk"
  else
    ac_cv_have_lpsolve=no
    LPSOLVE_CPPFLAGS=
    LPSOLVE_LDFLAGS=
    LPSOLVE_LIBS=
  fi
  ]
)
