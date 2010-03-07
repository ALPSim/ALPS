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
        # Be sure to have absolute paths.
        case $fop_dir in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute directory name for --with-fop : $fop_dir]);;
        esac
      fi
    fi
    ]
  )

  found=no

  if test "x$fop" != xno; then
    AC_MSG_CHECKING([for Apache FOP])
    if test -n "$fop_dir"; then
      for s in fop fop.sh
      do
      	if test -f "$fop_dir/$s"; then
          found=yes
	  fop_bin="$s"
        else
          fop_dir=no
	  fop_bin=no
        fi
      done
    else
      fop_dir=no
      fop_bin=no
      for d in $prefix $prefix/src $HOME $HOME/src /usr/local /usr/local/src /usr $BOOST_DIR/tools/boostbook
      do
        for s in fop fop-0.20.5 bin
        do
	  for k in fop fop.sh
	  do
            if test -f "$d/$s/$k"; then
              found=yes
              fop_dir="$d/$s"
	      fop_bin="$k"
              break
            fi
	  done
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
    ac_cv_fop_bin="$fop_bin"
    ac_cv_fop=$fop_dir/$fop_bin"
  else
    ac_cv_have_fop=no
    ac_cv_fop_dir=
    ac_cv_fop_bin=
    ac_cv_fop=
  fi
  ]
)
