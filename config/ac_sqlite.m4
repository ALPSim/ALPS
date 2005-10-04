#
# check for SQLite Library
#

AC_DEFUN([AC_SQLITE],
  [
  AC_ARG_WITH(sqlite-incdir,
    AC_HELP_STRING([--with-sqlite-incdir=DIR],[SQLite include directory]),
    [sqlite_incdir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`])
  AC_ARG_WITH(sqlite-libdir,
    AC_HELP_STRING([--with-sqlite-libdir=DIR],[SQLite library directory]),
    [sqlite_libdir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`])

  if test -z "$sqlite_incdir"; then
    for d in /usr/local; do
      if test -f "$d/include/sqlite3.h"; then
        sqlite_incdir="$d/include"
        if test -z "$sqlite_libdir"; then
          sqlite_libdir="$d/lib"
        fi
      fi
    done
  fi

  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_LIBS="$LIBS"

  if test -n "$sqlite_incdir" && test -d "$sqlite_incdir"; then
    CPPFLAGS="$CPPFLAGS -I$sqlite_incdir"
    ac_cv_sqlite_cppflags="-I$sqlite_incdir"
  fi
  if test -n "$sqlite_libdir" && test -d "$sqlite_libdir"; then
    LDFLAGS="$LDFLAGS -L$sqlite_libdir"
    ac_cv_sqlite_ldflags="-L$sqlite_libdir"
  fi

  AC_CHECK_HEADER([sqlite3.h])
  AC_CHECK_LIB(sqlite3, sqlite3_libversion_number)

  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
  LIBS="$ac_save_LIBS"
  AC_LANG_RESTORE

  if test "$ac_cv_header_sqlite3_h" = yes && test "$ac_cv_lib_sqlite3_sqlite3_libversion_number" = yes; then
    ac_cv_have_sqlite=yes
    ac_cv_sqlite_libs="-lsqlite3"
  else
    ac_cv_have_sqlite=no
    ac_cv_sqlite_cppflags=
    ac_cv_sqlite_ldflags=
    ac_cv_sqlite_libs=
  fi

  AC_SUBST(SQLITE_CPPFLAGS)
  AC_SUBST(SQLITE_LDFLAGS)
  AC_SUBST(SQLITE_LIBS)
  test -n $ac_cv_sqlite_cppflags && SQLITE_CPPFLAGS="$ac_cv_sqlite_cppflags"
  test -n $ac_cv_sqlite_ldflags && SQLITE_LDFLAGS="$ac_cv_sqlite_ldflags"
  test -n $ac_cv_sqlite_libs && SQLITE_LIBS="$ac_cv_sqlite_libs"
  ]
)
