dnl
dnl libxslt
dnl

AC_DEFUN([AC_LIBXSLT],
  [
  AC_ARG_WITH(libxslt,
    AC_HELP_STRING([--with-libxslt=DIR],[specify the path to libxslt library]),
    [
    if test "x$withval" = "xno"; then
      ac_cv_have_libxslt=no
    else
      if test "x$withval" != "xyes"; then
        libxslt_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $libxslt_dir in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute directory name for --with-libxslt : $libxslt_dir]);;
        esac
      fi
    fi
    ]
  )
  
  if test "$ac_cv_have_libxslt" != no; then
    test -z "$ac_cv_have_libxml" && AC_LIBXML
    test "$ac_cv_have_libxml" = no && ac_cv_have_libxslt=no
  fi

  if test -n "$ac_cv_have_libxslt"; then
    AC_MSG_CHECKING([for libxslt library])
    AC_MSG_RESULT([$ac_cv_have_libxslt])
  else
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS="$CPPFLAGS"
    ac_save_LDFLAGS="$LDFLAGS"
    ac_save_LIBS="$LIBS"

    cppflags=
    ldflags=
    libs=

    found=

    cppflags_try=
    if test -n "$libxslt_dir" && \
       test -f "$libxslt_dir/include/libxslt/xslt.h"; then
      AC_MSG_CHECKING([for libxslt include dir])
      AC_MSG_RESULT([$libxslt_dir/include/libxslt])
      cppflags_try="-I$libxslt_dir/include"
    fi
    CPPFLAGS="$ac_save_CPPFLAGS $ac_cv_libxml_cppflags $cppflags_try"
    AC_CHECK_HEADER([libxslt/xslt.h], [cppflags="$cppflags_try"], [found=no])
  
    case "$host_os" in
      darwin* )
        # for Mac OS X
        if test -z "$found"; then
          CPPFLAGS="$ac_save_CPPFLAGS $ac_cv_libxml_cppflags $cppflags"
          LDFLAGS="$ac_save_LDFLAGS $ac_cv_libxml_ldflags -framework libxslt"
          LIBS="$ac_cv_libxml_libs $ac_save_LIBS"
          AC_MSG_CHECKING([for xsltApplyStylesheet in -framework libxslt])
          AC_TRY_LINK([extern "C" char xsltApplyStylesheet();],
            [xsltApplyStylesheet();],
            [AC_MSG_RESULT(yes); ldflags="-framework libxslt"; found=yes],
            [AC_MSG_RESULT(no)])
        fi
        break ;;
      * )
        break ;;
    esac

    # check for libxslt.a
    if test -z "$found"; then
      CPPFLAGS="$ac_save_CPPFLAGS $ac_cv_libxml_cppflags $cppflags"
      ldflags_try=
      if test -n "$libxslt_dir" && test -d "$libxslt_dir/lib"; then
        AC_MSG_CHECKING([for libxslt library dir])
        AC_MSG_RESULT([$libxslt_dir/lib])
        ldflags_try="-L$libxslt_dir/lib"
      fi
      LDFLAGS="$ac_save_LDFLAGS $ac_cv_libxml_ldflags $ldflags_try"
      LIBS="$ac_cv_libxml_libs $ac_save_LIBS"
      AC_CHECK_LIB([xslt], [xsltApplyStylesheet],
        [ldflags="$ldflags_try"; libs="-lxslt"; found=yes])
    fi
  
    if test "$found" = yes; then
      ac_cv_libxslt_cppflags="$ac_cv_libxml_cppflags $cppflags"
      ac_cv_libxslt_ldflags="$ac_cv_libxml_ldflags $ldflags"
      ac_cv_libxslt_libs="$libs $ac_cv_libxml_libs"
      ac_cv_have_libxslt=yes
      AC_DEFINE(HAVE_LIBXSLT, [], [Define if you have a libxslt library.])
    else
      ac_cv_have_libxslt=no
    fi

    CPPFLAGS=$ac_save_CPPFLAGS
    LDFLAGS=$ac_save_LDFLAGS
    LIBS=$ac_save_LIBS
    AC_LANG_RESTORE
  fi
  ]
)


dnl
dnl Xalan C++ library
dnl

AC_DEFUN([AC_LIBXALAN],
  [
  AC_ARG_WITH(libxalan,
    AC_HELP_STRING([--with-libxalan=DIR],
      [specify the path to Xalan C++ library]),
    [
    if test "x$withval" = "xno"; then
      ac_cv_have_libxalan=no
    else
      if test "x$withval" != "xyes"; then
        libxalan_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $libxalan_dir in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute directory name for --with-libxalan : $libxalan_dir]);;
        esac
      fi
    fi
    ]
  )

  if test "$ac_cv_have_libxalan" != no; then
    test -z "$ac_cv_have_libxerces" && AC_XERCES
    test "$ac_cv_have_xerces" = no && ac_cv_have_libxalan=no
  fi

  if test -n "$ac_cv_have_libxalan"; then
    AC_MSG_CHECKING([for Xalan C++ library])
    AC_MSG_RESULT([$ac_cv_have_libxalan])
  else
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS="$CPPFLAGS"
    ac_save_LDFLAGS="$LDFLAGS"
    ac_save_LIBS="$LIBS"

    cppflags=
    ldflags=
    libs=

    found=

    cppflags_try=
    if test -n "$libxalan_dir" && \
       test -f "$libxalan_dir/include/xalanc/XSLT/XSLTInit.hpp"; then
      AC_MSG_CHECKING([for Xalan C++ include dir])
      AC_MSG_RESULT([$libxalan_dir/include/xalanc])
      cppflags_try="-I$libxalan_dir/include"
    fi
    CPPFLAGS="$ac_save_CPPFLAGS $ac_cv_xerces_cppflags $cppflags_try"
    AC_CHECK_HEADER([xalanc/XSLT/XSLTInit.hpp],
      [cppflags="$cppflags_try"], [found=no])

    if test -z "$found"; then
      CPPFLAGS="$ac_save_CPPFLAGS $ac_cv_xerces_cppflags $cppflags"
      ldflags_try=
      if test -n "$xalan_dir" && test -d "$xalan_dir/lib"; then
        AC_MSG_CHECKING([for Xalan C++ library dir])
        AC_MSG_RESULT([$xalan_dir/lib])
        ldflags_try="-L$xalan_dir/lib"
      fi
      LDFLAGS="$ac_save_LDFLAGS $ldflags_try"
      LIBS="-lxalan-c $ac_save_LIBS"
      AC_MSG_CHECKING([for XMLPlatoformUtils::Initialize() in -lxalan-c])
      AC_TRY_LINK(
        [#include <xalanc/XalanTransformer/XalanTransformer.hpp>],
        [XALAN_USING_XALAN(XalanTransformer); XalanTransformer::initialize();],
        [AC_MSG_RESULT(yes); ldflags="$ldflags_try"; libs="-lxalan-c";
         found=yes],
        [AC_MSG_RESULT(no)])
    fi

    if test "$found" = yes; then
      ac_cv_libxalan_cppflags="$ac_cv_libxerces_cppflags $cppflags"
      ac_cv_libxalan_ldflags="$ac_cv_libxerces_ldflags $ldflags"
      ac_cv_libxalan_libs="$libs $ac_cv_libxerces_libs"
      ac_cv_have_libxalan=yes
      AC_DEFINE(HAVE_LIBXALAN, [], [Define if you have an Xalan C++ library.])
    else
      ac_cv_have_libxalan=no
    fi

    CPPFLAGS=$ac_save_CPPFLAGS
    LDFLAGS=$ac_save_LDFLAGS
    LIBS=$ac_save_LIBS
    AC_LANG_RESTORE
  fi
  ]
)


dnl
dnl xsltproc 
dnl

AC_DEFUN([AC_XSLTPROC],
  [
  AC_ARG_WITH(xsltproc,
    AC_HELP_STRING([--with-xsltproc=PATH],
      [specify the path to xsltproc]),
    [
    if test "x$withval" = "xno"; then
      ac_cv_have_xsltproc=no
    else
      if test "x$withval" != "xyes"; then
        xsltproc_path=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $xsltproc_path in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute path name for --with-xsltproc : $xsltproc_path]);;
        esac
      fi
    fi
    ]
  )

  if test -n "$ac_cv_have_xsltproc"; then
    AC_MSG_CHECKING([for xsltproc])
    AC_MSG_RESULT([$ac_cv_have_xsltproc])
  else
    if test -z "$xsltproc_path"; then
      AC_PATH_PROG([ac_cv_prog_xsltproc], [xsltproc])
    else
      if test -x "$xsltproc_path"; then
        AC_MSG_CHECKING([for $xsltproc_path])
        AC_MSG_RESULT([yes])
        ac_cv_prog_xsltproc="$xsltproc_path"
      else
        AC_PATH_PROG([ac_cv_prog_xsltproc], [$xsltproc_path])
      fi
    fi
    if test -n "$ac_cv_prog_xsltproc"; then
      ac_cv_have_xsltproc=yes
      AC_DEFINE_UNQUOTED([XSLTPROC], ["$ac_cv_prog_xsltproc"])
    else
      ac_cv_have_xsltproc=no
    fi
  fi
  ]
)


dnl
dnl xalan 
dnl

AC_DEFUN([AC_XALAN],
  [
  AC_ARG_WITH(xalan,
    AC_HELP_STRING([--with-xalan=PATH],
      [specify the path to xalan]),
    [
    if test "x$withval" = "xno"; then
      ac_cv_have_xalan=no
    else
      if test "x$withval" != "xyes"; then
        xalan_path=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $xalan_path in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute path name for --with-xalan : $xalan_path]);;
        esac
      fi
    fi
    ]
  )

  if test -n "$ac_cv_have_xalan"; then
    AC_MSG_CHECKING([for xalan])
    AC_MSG_RESULT([$ac_cv_have_xalan])
  else
    if test -z "$xalan_path"; then
      AC_PATH_PROG([ac_cv_prog_xalan], [xalan])
    else
      if test -x "$xalan_path"; then
        AC_MSG_CHECKING([for $xalan_path])
        AC_MSG_RESULT([yes])
        ac_cv_prog_xalan="$xalan_path"
      else
        AC_PATH_PROG([ac_cv_prog_xalan], [$xalan_path])
      fi
    fi
    if test -n "$ac_cv_prog_xalan"; then
      ac_cv_have_xalan=yes
      AC_DEFINE_UNQUOTED([XALAN], ["$ac_cv_prog_xalan"])
    else
      ac_cv_have_xalan=no
    fi
  fi
  ]
)


dnl
dnl XSLT processor
dnl

AC_DEFUN([AC_XSLTP],
  [
  AC_ARG_WITH(xslt-proc,
    AC_HELP_STRING([--with-xslt-proc=TYPE],
      [specify XSLT processor (TYPE = libxslt, libxalan, xsltproc, xalan)]),
    [
    case "x$withval" in
      xlibxslt* )
        ac_cv_xsltp=libxslt
        ;;
      xlibxalan* )
        ac_cv_xsltp=libxalan
        ;;
      xxsltproc* )
        ac_cv_xsltp=xsltproc
        ;;
      xxalan* )
        ac_cv_xsltp=xalan
        ;;
      *)
        AC_MSG_ERROR([unknown XSLT processor type])
        ;;
    esac
    ]
  )

  dnl libxslt
  if test -z "$ac_cv_xsltp" || test "$ac_cv_xsltp" = libxslt; then
    AC_LIBXSLT
    if test "$ac_cv_have_libxslt" = yes; then
      ac_cv_have_xsltp=yes
      ac_cv_xsltp=libxslt
      ac_cv_xsltp_cppflags="$ac_cv_libxslt_cppflags"
      ac_cv_xsltp_ldflags="$ac_cv_libxslt_ldflags"
      ac_cv_xsltp_libs="$ac_cv_libxslt_libs"
    else
      test -n "$ac_cv_xsltp" && AC_MSG_ERROR([XSLT processor not found.])
    fi
  fi

  dnl libxalan
  if test -z "$ac_cv_xsltp" || test "$ac_cv_xsltp" = libxalan; then
    AC_LIBXALAN
    if test "$ac_cv_have_libxalan" = yes; then
      ac_cv_have_xsltp=yes
      ac_cv_xsltp=libxalan
      ac_cv_xsltp_cppflags="$ac_cv_libxalan_cppflags"
      ac_cv_xsltp_ldflags="$ac_cv_libxalan_ldflags"
      ac_cv_xsltp_libs="$ac_cv_libxalan_libs"
    else
      test -n "$ac_cv_xsltp" && AC_MSG_ERROR([XSLT processor not found.])
    fi
  fi

  dnl xsltproc
  if test -z "$ac_cv_xsltp" || test "$ac_cv_xsltp" = xsltproc; then
    AC_XSLTPROC
    if test "$ac_cv_have_xsltproc" = yes; then
      ac_cv_have_xsltp=yes
      ac_cv_xsltp=xsltproc
    else
      test -n "$ac_cv_xsltp" && AC_MSG_ERROR([XSLT processor not found.])
    fi
  fi

  dnl xalan
  if test -z "$ac_cv_xsltp" || test "$ac_cv_xsltp" = xalan; then
    AC_XALAN
    if test "$ac_cv_have_xalan" = yes; then
      ac_cv_have_xsltp=yes
      ac_cv_xsltp=xalan
    else
      test -n "$ac_cv_xsltp" && AC_MSG_ERROR([XSLT processor not found.])
    fi
  fi
  ]
)
