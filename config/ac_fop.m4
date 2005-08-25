#
# check for Apache FOP
#

AC_DEFUN([AC_FOP],
  [
  AC_ARG_WITH(fop,
    AC_HELP_STRING([--with-fop=DIR],
                   [path to Apache FOP]),
    [
    if test "x$withval" = xno; then
      fop=no
    else
      fop=yes
      if test "x$withval" != xyes; then
        fop_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )

  found=no

  if test "x$fop" != xno; then
    AC_MSG_CHECKING([for Apache FOP])
    if test -n "$fop_dir"; then
      if test -f "$fop_dir/fop.sh"; then
        found=yes
      else
        fop_dir=no
      fi
    else
      fop_dir=no
      for d in $prefix $prefix/src $HOME $HOME/src /usr/local /usr/local/src
      do
        for s in fop-0.20.5
        do
          if test -f "$d/$s/fop.sh"; then
            found=yes
            fop_dir="$d/$s"
            break
          fi
        done
        if test "$found" = yes; then
          break
        fi
      done
    fi
    AC_MSG_RESULT([$fop_dir])
  fi

  if test "$found" = yes; then
    ac_cv_have_fop=yes
    ac_cv_fop_dir="$fop_dir"
  else
    ac_cv_have_fop=no
    ac_cv_fop_dir=
  fi
  ]
)
