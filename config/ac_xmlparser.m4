AC_DEFUN([AC_XMLPARSER],
  [
  AC_SUBST(XML_CPPFLAGS)
  AC_SUBST(XML_LDFLAGS)
  AC_SUBST(XML_LIBS)

  AC_ARG_WITH(expat,
    AC_HELP_STRING([--with-expat=DIR],[use Expat XML parser by James Clark]),
    [
    if test "x$withval" = "xno"; then
      expat=no
    else
      if test -n "$xml_parser"; then
        AC_MSG_ERROR([two (or more) XML parsers are specified.])
      fi
      expat=yes
      xml_parser="expat"
      if test "x$withval" != "xyes"; then
        expat_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )

  AC_ARG_WITH(xerces,
    AC_HELP_STRING([--with-xerces=DIR],[use Xerces C++ XML parser by Apache Software Foundation]),
    [
    if test "x$withval" = "xno"; then
      xerces=no
    else
      if test -n "$xml_parser"; then
        AC_MSG_ERROR([two (or more) XML parsers are specified.])
      fi
      xerces=yes
      xml_parser="xerces"
      if test "x$withval" != "xyes"; then
        xerces_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )

  if test "$xerces" != "no"; then
    if test -z "$xml_parser" || test "$xml_parser" = "xerces"; then
      AC_MSG_CHECKING([for Xerces C++ root directory])
      if test "x$xerces_dir" = "x"; then
        for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local/src
        do
          if test -f "$d/include/xercesc/parsers/SAXParser.hpp"; then
            xerces_dir="$d"
            break
          fi
          if test -f "$d/xercesc/include/xercesc/parsers/SAXParser.hpp"; then
            xerces_dir="$d/xercesc"
            break
          fi
        done
        if test -n "$xerces_dir"; then
          AC_MSG_RESULT([$xerces_dir])
        else
          AC_MSG_RESULT()
        fi
      else
        AC_MSG_RESULT([$xerces_dir])
        if test -f "$xerces_dir/include/xercesc/parsers/SAXParser.hpp"; then :; else
          AC_MSG_ERROR([$xerces_dir/include/xercesc/parsers/SAXParser.hpp not found])
        fi
      fi

      if test -n "$xerces_dir"; then
        XML_CPPFLAGS="-I$xerces_dir/include"
        XML_LDFLAGS="-L$xerces_dir/lib"
      fi

      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      ac_save_CPPFLAGS=$CPPFLAGS
      ac_save_LDFLAGS=$LDFLAGS
      ac_save_LIBS=$LIBS
      CPPFLAGS="$XML_CPPFLAGS $CPPFLAGS"
      LDFLAGS="$XML_LDFLAGS $LDFLAGS"

      AC_CHECK_HEADER([xercesc/parsers/SAXParser.hpp],,
                      [
                      if test "$xerces" = yes; then
                        AC_MSG_ERROR([check for Xerces C++ XML parser failed])
                      else
                        xerces=no
                      fi
                      ])

      if test "$xerces" != no; then
        AC_MSG_CHECKING([for Xerces C++ version])
        xerces_version_maj=`grep "^#define XERCES_VERSION_MAJOR" $xerces_dir/include/xercesc/util/XercesVersion.hpp | awk "{print $ 3}" | sed 's/\"//g'`
	xerces_version_min=`grep "^#define XERCES_VERSION_MINOR" $xerces_dir/include/xercesc/util/XercesVersion.hpp | awk "{print $ 3}" | sed 's/\"//g'`
	xerces_version_rev=`grep "^#define XERCES_VERSION_REVISION" $xerces_dir/include/xercesc/util/XercesVersion.hpp | awk "{print $ 3}" | sed 's/\"//g'`
	if test -n "$xerces_version_maj"; then
 	  xerces_version="$xerces_version_maj.$xerces_version_min.$xerces_version_rev"
          AC_MSG_RESULT([$xerces_version])
	else
          AC_MSG_RESULT([unknown])
        fi

        found=no

        if test "$found" = no; then
          XML_LIBS="-lxerces-c"
          LIBS="$XML_LIBS $ac_save_LIBS"
          AC_MSG_CHECKING([for SAXParser in $XML_LIBS])
          AC_TRY_LINK([#include <xercesc/parsers/SAXParser.hpp>],
                      [XERCES_CPP_NAMESPACE_QUALIFIER SAXParser();],
                      [AC_MSG_RESULT(yes); found=yes],
                      AC_MSG_RESULT(no))
        fi

        if test "$found" = no; then
          XML_LIBS="-lxerces-c$xerces_version"
          LIBS="$XML_LIBS $ac_save_LIBS"
          AC_MSG_CHECKING([for SAXParser in $XML_LIBS])
          AC_TRY_LINK([#include <xercesc/parsers/SAXParser.hpp>],
                      [XERCES_CPP_NAMESPACE_QUALIFIER SAXParser();],
                      [AC_MSG_RESULT(yes); found=yes],
                      AC_MSG_RESULT(no))
        fi

        if test "$found" = no; then
          if test "$xerces" = yes; then
            AC_MSG_ERROR([check for Xerces C++ XML library failed])
          fi
          xerces=no
        else
          xml_parser="xerces"
        fi
      fi

      CPPFLAGS=$ac_save_CPPFLAGS
      LDFLAGS=$ac_save_LDFLAGS
      LIBS=$ac_save_LIBS
      AC_LANG_RESTORE
    fi
  fi

  if test "$expat" != "no"; then
    if test -z "$xml_parser" || test "$xml_parser" = "expat"; then
      AC_MSG_CHECKING([for expat root directory])
      if test "x$expat_dir" = "x"; then
        for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local/src
        do
          if test -f "$d/include/expat.h"; then
            expat_dir="$d"
            break
          fi
          if test -f "$d/expat/include/expat.h"; then
            expat_dir="$d/expat"
            break
          fi
        done
        if test -n "$expat_dir"; then
          AC_MSG_RESULT([$expat_dir])
        else
          AC_MSG_RESULT()
        fi
      else
        AC_MSG_RESULT([$expat_dir])
        if test -f "$expat_dir/include/expat.h"; then :; else
          AC_MSG_ERROR([$expat_dir/include/expat.h not found])
        fi
      fi

      if test -n "$expat_dir"; then
        XML_CPPFLAGS="-I$expat_dir/include"
        XML_LDFLAGS="-L$expat_dir/lib"
      fi

      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      ac_save_CPPFLAGS=$CPPFLAGS
      ac_save_LDFLAGS=$LDFLAGS
      ac_save_LIBS=$LIBS
      CPPFLAGS="$XML_CPPFLAGS $CPPFLAGS"
      LDFLAGS="$XML_LDFLAGS $LDFLAGS"

      AC_CHECK_HEADER([expat.h],,
                      [
                      if test "$expat" = yes; then
                        AC_MSG_ERROR([check for expat XML parser failed])
                      else
                        expat=no
                      fi
                      ])

      if test "$expat" != no; then
        found=no

        if test "$found" = no; then
          XML_LIBS="-lexpat"
          LIBS="$XML_LIBS $LIBS"
          AC_MSG_CHECKING([for XML_ParserCreate() in $XML_LIBS])
          AC_TRY_LINK([#include <expat.h>],
                      [XML_Parser p = XML_ParserCreate(NULL);],
                      [AC_MSG_RESULT(yes); found=yes],
                      AC_MSG_RESULT(no))
        fi

        if test "$found" = no; then
          if test "$expat" = yes; then
            AC_MSG_ERROR([check for expat XML library failed])
          fi
          expat=no
        else
          xml_parser="expat"
        fi
      fi

      CPPFLAGS=$ac_save_CPPFLAGS
      LDFLAGS=$ac_save_LDFLAGS
      LIBS=$ac_save_LIBS
      AC_LANG_RESTORE
    fi
  fi
  
  if test -n "$xml_parser"; then
    ac_cv_xml_parser=$xml_parser
    ac_cv_xml_cppflags="$XML_CPPFLAGS"
    ac_cv_xml_ldflags="$XML_LDFLAGS"
    ac_cv_xml_libs="$XML_LIBS"
    if test "$xml_parser" = "expat"; then
      AC_DEFINE(HAVE_EXPAT_PARSER)
    fi
    if test "$xml_parser" = "xerces"; then
      AC_DEFINE(HAVE_XERCES_PARSER)
    fi
    AC_MSG_NOTICE([using XML parser : $xml_parser])
  else
    ac_cv_xml_parser=native
    ac_cv_xml_cppflags=
    ac_cv_xml_ldflags=
    ac_cv_xml_libs=
    XML_CPPFLAGS=
    XML_LDFLAGS=
    XML_LIBS=
    AC_MSG_NOTICE([XML parser not found])
    AC_MSG_NOTICE([using native XML parser])
  fi
  ]
)
