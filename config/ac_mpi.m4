AC_DEFUN([AC_MPI],
  [
  AC_SUBST(MPI_CPPFLAGS)
  AC_SUBST(MPI_LDFLAGS)
  AC_SUBST(MPI_LIBS)
  AC_SUBST(OBJ_MPI)
  AC_SUBST(LIB_MPI)
  AC_SUBST(EXAMPLE_MPI)
  
  AC_ARG_WITH(mpi,
    AC_HELP_STRING([--with-mpi=DIR],[use MPI library]),
    [
    if test "x$withval" = "xno"; then
      mpi=no
    else
      mpi=yes
      if test "x$withval" != "xyes"; then
        mpi_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
)
  
  AC_ARG_WITH(mpi-incdir,
    AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory]),
    [mpi_incdir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`])
  AC_ARG_WITH(mpi-libdir,
    AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory]),
    [mpi_libdir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`])
  AC_ARG_WITH(mpi-libs,
    AC_HELP_STRING([--with-mpi-libs=LIBS],[MPI libraries]),
    [mpi_libs="$withval"])
  
  if test "$mpi" != no; then
    AC_MSG_CHECKING([for MPI root directory])
    ll_root="/usr/lpp/ppe.poe"
    if test "x$MP_PREFIX" != x; then
      ll_root="$MP_PREFIX/ppe.poe"
    fi
    if test -z "$mpi_dir"; then
      for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local /usr/local/src $ll_root; do
        if test -f "$d/include/mpi.h" && test -d "$d/lib"; then
          mpi_dir="$d"
          break
        fi
        if test -f "$d/mpi/include/mpi.h" && test -d "$d/mpi/lib"; then
          mpi_dir="$d/mpi"
          break
        fi
        if test -f "$d/mpich/include/mpi.h" && test -d "$d/mpich/lib"; then
          mpi_dir="$d/mpich"
          break
        fi
        if test -f "$d/lammpi/include/mpi.h" && test -d "$d/lammpi/lib"; then
          mpi_dir="$d/lammpi"
          break
        fi
      done
      if test -n "$mpi_dir"; then
        AC_MSG_RESULT([$mpi_dir])
      else
        AC_MSG_RESULT([yes])
      fi
    else
      AC_MSG_RESULT([$mpi_dir])
      if test -f "$mpi_dir/include/mpi.h"; then :; else
        AC_MSG_ERROR([$mpi_dir/include/mpih.h not found])
      fi
      if test -d "$mpi_dir/lib"; then :; else
        AC_MSG_ERROR([$mpi_dir/lib not found])
      fi
    fi

    if test -n "$mpi_incdir"; then
      AC_MSG_CHECKING([for MPI header directory])
      AC_MSG_RESULT([$mpi_incdir])
      if test -f "$mpi_incdir/mpi.h"; then :; else
        AC_MSG_ERROR([$mpi_incdir/mpi.h not found])
      fi
    else
      if test -n "$mpi_dir"; then
        mpi_incdir="$mpi_dir/include"
        AC_MSG_CHECKING([for MPI header directory])
        AC_MSG_RESULT([$mpi_incdir])
      fi
    fi

    if test -n "$mpi_libdir"; then
      AC_MSG_CHECKING([for MPI library directory])
      AC_MSG_RESULT([$mpi_libdir])
      if test -d "$mpi_libdir"; then :; else
        AC_MSG_ERROR([$mpi_libdir not found])
      fi
    else
      if test -n "$mpi_dir"; then
        mpi_libdir="$mpi_dir/lib"
        AC_MSG_CHECKING([for MPI library directory])
        AC_MSG_RESULT([$mpi_libdir])
      fi
    fi

    if test -n "$mpi_incdir"; then
      if test -d "$mpi_incdir/mpi2c++"; then
        MPI_CPPFLAGS="-I$mpi_incdir -I$mpi_incdir/mpi2c++" # for LAM MPI
      else
        MPI_CPPFLAGS="-I$mpi_incdir"
      fi
    fi
    if test -n "$mpi_libdir"; then
      MPI_LDFLAGS="-L$mpi_libdir"
    fi
    
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS=$CPPFLAGS
    ac_save_LDFLAGS=$LDFLAGS
    ac_save_LIBS=$LIBS
    CPPFLAGS="$MPI_CPPFLAGS $CPPFLAGS"
    LDFLAGS="$MPI_LDFLAGS $LDFLAGS"

    AC_CHECK_HEADER([mpi.h],,
                    [
                    if test "$mpi" = yes; then
                      AC_MSG_ERROR([check for MPI library failed])
                    else
                      mpi=no
                    fi
                    ])

    if test "$mpi" != no; then
      found=no

      if test -n "$mpi_libs"; then
        MPI_LIBS="$mpi_libs"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    [AC_MSG_RESULT(no);
                     AC_MSG_ERROR([check for MPI library failed])])
      fi

      if test "$found" = no; then
        if test "$ac_cv_compiler" = ibm32 || test "$ac_cv_compiler" = ibm64; then
	  # for IBM LoadLeveler
          MPI_LIBS="-binitfini:poe_remote_main -lmpi -lvtd"
          LIBS="$MPI_LIBS $ac_save_LIBS"
          AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
          AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                      [AC_MSG_RESULT(yes); found=yes],
                      AC_MSG_RESULT(no))
        fi
      fi
      if test "$found" = no; then
        MPI_LIBS=
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize()])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi
      if test "$found" = no; then
        MPI_LIBS="-lmpi"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi
      if test "$found" = no; then
        MPI_LIBS="-lmpi -lmpi++"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi
      if test "$found" = no; then
        MPI_LIBS="-lmpich"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi
      if test "$found" = no; then
        MPI_LIBS="-llammpi++ -llammpio -lpmpi -llamf77mpi -lmpi -llam -lnsl -lutil"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi
      if test "$found" = no; then
        MPI_LIBS="-llammpi++ -llammpio -lpmpi -llamf77mpi -lmpi -llam"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi
      if test "$found" = no; then
        # for lammpi 7.1 in Fedora Core 3
        MPI_LIBS="-llammpi++ -llammpio -llamf77mpi -lmpi -llam -ldl"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi
      if test "$found" = no; then
        # for lammpi 7.0.2 on MacOSX
        MPI_LIBS="-llammpi++ -llammpio -lmpi -llam"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi
      if test "$found" = no; then
        MPI_LIBS="-lmpi++ -lmpi -ltstdio -ltrillium -largs -lt"
        LIBS="$MPI_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for MPI_Finalize() in $MPI_LIBS])
        AC_TRY_LINK([#include <mpi.h>],[MPI_Finalize();],
                    [AC_MSG_RESULT(yes); found=yes],
                    AC_MSG_RESULT(no))
      fi

      if test "$found" = no; then
        if test "$mpi" = yes; then
          AC_MSG_ERROR([check for MPI library failed])
        fi
        mpi=no
      fi

      CPPFLAGS=$ac_save_CPPFLAGS
      LDFLAGS=$ac_save_LDFLAGS
      LIBS=$ac_save_LIBS
      AC_LANG_RESTORE
    fi
  fi
  
  if test "$mpi" != no; then
    ac_cv_have_mpi=yes
    ac_cv_mpi_cppflags="$MPI_CPPFLAGS"
    ac_cv_mpi_ldflags="$MPI_LDFLAGS"
    ac_cv_mpi_libs="$MPI_LIBS"
    AC_DEFINE(HAVE_MPI)
    OBJ_MPI='$(OBJ_MPI)'
    LIB_MPI='$(LIB_MPI)'
    EXAMPLE_MPI='$(EXAMPLE_MPI)'
    AC_MSG_NOTICE([enabling MPI support])
  else
    ac_cv_have_mpi=no
    ac_cv_mpi_cppflags=
    ac_cv_mpi_ldflags=
    ac_cv_mpi_libs=
    MPI_CPPFLAGS=
    MPI_LDFLAGS=
    MPI_LIBS=
    AC_MSG_NOTICE([disabling MPI support])
  fi
  ]
)
