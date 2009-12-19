#
# check for IETL Library
#

AC_DEFUN([AC_IETL],
  [
  AC_SUBST(IETL_CPPFLAGS)
    
  AC_ARG_WITH(ietl,
    AC_HELP_STRING([--with-ietl=DIR],
                   [path to IETL library]),
    [
    if test "x$withval" = xno; then
      ietl=no
    else
      ietl=yes
      if test "x$withval" != xyes; then
        ietl_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $ietl_dir in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute directory name for --with-ietl : $ietl_dir]);;
        esac
      fi
    fi
    ]
  )

  found=no

  if test "x$ietl" != xno; then
    AC_MSG_CHECKING([for IETL library])
    if test -n "$ietl_dir"; then
      if test -f "$ietl_dir/ietl/lanczos.h"; then
        found=yes
      else
        ietl_dir=no
      fi
    else
      ietl_dir=no
      for d in $prefix $prefix/src $HOME $HOME/src /usr/local /usr/local/src
      do
        for s in include ietl-2.0 ietl
        do
          if test -f "$d/$s/ietl/lanczos.h"; then
            found=yes
            ietl_dir="$d/$s"
            break
          fi
        done
        if test "$found" = yes; then
          break
        fi
      done
    fi
    AC_MSG_RESULT([$ietl_dir])
  fi

  if test "$found" = yes; then
    ac_cv_have_ietl=yes
    ac_cv_ietl_dir="$ietl_dir"
    IETL_CPPFLAGS="-I$ietl_dir"
  else
    ac_cv_have_ietl=no
    ac_cv_ietl_dir=
    IETL_CPPFLAGS=
  fi

  if test "$ac_cv_have_ietl" = yes; then
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $IETL_CPPFLAGS"
    AC_CHECK_HEADER([ietl/lanczos.h])
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_LANG_RESTORE
  fi
  ]
)
