#
# check for DocBook
#

AC_DEFUN([AC_DOCBOOK_DTD],
  [
  AC_ARG_WITH(docbook-dtd,
    AC_HELP_STRING([--with-docbook-dtd=DIR],
                   [path to DocBook DTD]),
    [
    if test "x$withval" = xno; then
      docbook_dtd=no
    else
      docbook_dtd=yes
      if test "x$withval" != xyes; then
        docbook_dtd_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $docbook_dtd_dir in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute directory name for --with-docbook-dtd : $docbook_dtd_dir]);;
        esac
      fi
    fi
    ]
  )

  found=no

  if test "x$docbook_dtd" != xno; then
    AC_MSG_CHECKING([for DocBook DTD])
    if test -n "$docbook_dtd_dir"; then
      if test -f "$docbook_dtd_dir/docbookx.dtd"; then
        found=yes
      else
        docbook_dtd_dir=no
      fi
    else
      docbook_dtd_dir=no
      for d in $prefix $prefix/src $HOME $HOME/src /usr/local /usr/local/src $BOOST_DIR/tools/boostbook
      do
        for s in docbook-dtd docbook-dtd-4.4 docbook-dtd-4.2
        do
          if test -f "$d/$s/docbookx.dtd"; then
            found=yes
            docbook_dtd_dir="$d/$s"
            break
          fi
        done
        if test "$found" = yes; then
          break
        fi
      done
    fi
    AC_MSG_RESULT([$docbook_dtd_dir])
  fi

  if test "$found" = yes; then
    ac_cv_have_docbook_dtd=yes
    ac_cv_docbook_dtd_dir="$docbook_dtd_dir"
  else
    ac_cv_have_docbook_dtd=no
    ac_cv_docbook_dtd_dir=
  fi
  ]
)

AC_DEFUN([AC_DOCBOOK_XSL],
  [
  AC_ARG_WITH(docbook-xsl,
    AC_HELP_STRING([--with-docbook-xsl=DIR],
                   [path to DocBook XSL]),
    [
    if test "x$withval" = xno; then
      docbook_xsl=no
    else
      docbook_xsl=yes
      if test "x$withval" != xyes; then
        docbook_xsl_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $docbook_xsl_dir in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute directory name for --with-docbook-xsl : $docbook_xsl_dir]);;
        esac
      fi
    fi
    ]
  )

  found=no

  if test "x$docbook_xsl" != xno; then
    AC_MSG_CHECKING([for DocBook XSL])
    if test -n "$docbook_xsl_dir"; then
      if test -f "$docbook_xsl_dir/fo/docbook.xsl"; then
        found=yes
      else
        docbook_xsl_dir=no
      fi
    else
      docbook_xsl_dir=no
      for d in $prefix $prefix/src $HOME $HOME/src /usr/local /usr/local/src $BOOST_DIR/tools/boostbook
      do
        for s in docbook-xsl docbook-xsl-1.69.1 docbook-xsl-1.68.1
        do
          if test -f "$d/$s/fo/docbook.xsl"; then
            found=yes
            docbook_xsl_dir="$d/$s"
            break
          fi
        done
        if test "$found" = yes; then
          break
        fi
      done
    fi
    AC_MSG_RESULT([$docbook_xsl_dir])
  fi

  if test "$found" = yes; then
    ac_cv_have_docbook_xsl=yes
    ac_cv_docbook_xsl_dir="$docbook_xsl_dir"
  else
    ac_cv_have_docbook_xsl=no
    ac_cv_docbook_xsl_dir=
  fi
  ]
)
