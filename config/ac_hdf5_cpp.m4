AC_DEFUN([AC_HDF5_CPP],
  [
  AC_SUBST(HDF5_CPPFLAGS)
  AC_SUBST(HDF5_LDFLAGS)
  AC_SUBST(HDF5_LIBS)

  AC_ARG_WITH(hdf5,
    AC_HELP_STRING([--with-hdf5=DIR],[use HDF5 library]),
    [
    if test "x$withval" = "xno"; then
      hdf5=no
    else
      hdf5=yes
      if test "x$withval" != "xyes"; then
        hdf5_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
        # Be sure to have absolute paths.
        case $hdf5_dir in
          [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
          *)  AC_MSG_ERROR([expected an absolute directory name for --with-hdf5 : $hdf5_dir]);;
        esac
      fi
    fi
    ]
  )

  # default is 'without HDF5'
  test -z "$hdf5" && hdf5=no

  if test "$hdf5" != "no"; then
    AC_MSG_CHECKING([for HDF5 c++ interface])
    if test "x$hdf5_dir" = "x"; then
      for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local/src /sw; do
        if test -f "$d/include/H5Cpp.h"; then
          hdf5_dir="$d"
          break
        fi
        if test -f "$d/hdf5/include/H5Cpp.h"; then
          hdf5_dir="$d/hdf5"
          break
        fi
      done
      if test -n "$hdf5_dir"; then
        AC_MSG_RESULT([$hdf5_dir])
      else
        AC_MSG_RESULT([yes])
      fi
    else
      AC_MSG_RESULT([$hdf5_dir])
      if test -f "$hdf5_dir/include/H5Cpp.h"; then :; else
        AC_MSG_ERROR([$hdf5_dir/include/H5Cpp.h not found])
      fi
    fi

    if test -n "$hdf5_dir"; then
      HDF5_CPPFLAGS="-I$hdf5_dir/include"
      HDF5_LDFLAGS="-L$hdf5_dir/lib"
    fi

    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS=$CPPFLAGS
    ac_save_LDFLAGS=$LDFLAGS
    ac_save_LIBS=$LIBS
    CPPFLAGS="$HDF5_CPPFLAGS $CPPFLAGS"
    LDFLAGS="$HDF5_LDFLAGS $LDFLAGS"

    AC_CHECK_HEADER([H5Cpp.h],,
                    [
                    if test "$hdf5" = yes; then
                      AC_MSG_ERROR([check for HDF5 library failed])
                    else
                      hdf5=no
                    fi
                    ])
                     

    if test "$hdf5" != no; then
      found=no

      if test "$found" = no; then
        HDF5_LIBS="-lhdf5_cpp -lhdf5"
        LIBS="$HDF5_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for H5File() in $HDF5_LIBS])
        AC_TRY_LINK([#include <H5Cpp.h>],
                    [#ifndef H5_NO_NAMESPACE
                     using namespace H5;
                     #endif
                     H5File file("tmp.h5",H5F_ACC_TRUNC);],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi

      if test "$found" = no; then
        HDF5_LIBS="-lhdf5_cpp -lhdf5 -lz"
        LIBS="$HDF5_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for H5File() in $HDF5_LIBS])
        AC_TRY_LINK([#include <H5Cpp.h>],
                    [#ifndef H5_NO_NAMESPACE
                     using namespace H5;
                     #endif
                     H5File file("tmp.h5",H5F_ACC_TRUNC);],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi

      if test "$found" = no; then
        if test "$hdf5" = yes; then
          AC_MSG_ERROR([check for HDF5 library failed])
        fi
        hdf5=no
      fi

      CPPFLAGS=$ac_save_CPPFLAGS
      LDFLAGS=$ac_save_LDFLAGS
      LIBS=$ac_save_LIBS
      AC_LANG_RESTORE
    fi
  fi

  if test "$hdf5" != no; then
    ac_cv_have_hdf5_cpp=yes
    ac_cv_hdf5_cppflags="$HDF5_CPPFLAGS"
    ac_cv_hdf5_ldflags="$HDF5_LDFLAGS"
    ac_cv_hdf5_libs="$HDF5_LIBS"
    AC_DEFINE(HAVE_HDF5_CPP)
    AC_MSG_NOTICE([enabling HDF5 c++ support])
  else
    ac_cv_have_hdf5_cpp=no
    ac_cv_hdf5_cppflags=
    ac_cv_hdf5_ldflags=
    ac_cv_hdf5_libs=
    HDF5_CPPFLAGS=
    HDF5_LDFLAGS=
    HDF5_LIBS=
    AC_MSG_NOTICE([disabling HDF5 c++ support])
  fi
  ]
)
