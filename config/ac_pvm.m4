
AC_DEFUN([AC_PVM],
  [
  AC_SUBST(PVM_CPPFLAGS)
  AC_SUBST(PVM_LDFLAGS)
  AC_SUBST(PVM_LIBS)
  AC_SUBST(OBJ_PVM)
  AC_SUBST(LIB_PVM)
  AC_SUBST(EXAMPLE_PVM)
  
  AC_ARG_WITH(pvm,
    AC_HELP_STRING([--with-pvm=DIR],[use PVM library]),
    [
    if test "x$withval" = "xno"; then
      pvm=no
    else
      pvm=yes
      if test "x$withval" != "xyes"; then
        pvm_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )

  AC_ARG_WITH(pvm-incdir,
    AC_HELP_STRING([--with-pvm-incdir=DIR],[PVM include directory]),
    [pvm_incdir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`])
  AC_ARG_WITH(pvm-libdir,
    AC_HELP_STRING([--with-pvm-libdir=DIR],[PVM library directory]),
    [pvm_libdir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`])
  AC_ARG_WITH(pvm-libs,
    AC_HELP_STRING([--with-pvm-libs=LIBS],[PVM libraries]),
    [pvm_libs="$withval"])
  
  if test "$pvm" != no; then
    AC_MSG_CHECKING([for PVM root directory])
    if test -z "$pvm_dir"; then
      for d in $HOME $HOME/src $prefix $prefix/src /usr/share /usr/array/PVM /usr/local /usr/local /usr/local/src; do
        if test -f "$d/include/pvm3.h" && test -d "$d/lib"; then
          pvm_dir="$d"
          break
        fi
        if test -f "$d/pvm/include/pvm3.h" && test -d "$d/pvm/lib"; then
          pvm_dir="$d/pvm"
          break
        fi
        if test -f "$d/pvm3/include/pvm3.h" && test -d "$d/pvm3/lib"; then
          pvm_dir="$d/pvm3"
          break
        fi
      done
      if test -n "$pvm_dir"; then
        AC_MSG_RESULT([$pvm_dir])
      else
        AC_MSG_RESULT([yes])
      fi
    else
      AC_MSG_RESULT([$pvm_dir])
      if test -f "$pvm_dir/include/pvm3.h"; then :; else
        AC_MSG_ERROR([$pvm_dir/include/pvmh3.h not found])
      fi
    fi

    AC_MSG_CHECKING([for PVM arch])
    if test -n "$PVM_ARCH"; then
      pvm_arch="$PVM_ARCH"
    else
      if test -x "$pvm_dir/lib/pvmgetarch"; then
        pvm_arch=`$pvm_dir/lib/pvmgetarch`
      fi
    fi
    if test -n "$pvm_arch"; then
      AC_MSG_RESULT([$pvm_arch])
    else
      AC_MSG_RESULT([unknown])
    fi
    
    if test -n "$pvm_incdir"; then
      AC_MSG_CHECKING([for PVM header directory])
      AC_MSG_RESULT([$pvm_incdir])
      if test -f "$pvm_incdir/pvm3.h"; then :; else
        AC_MSG_ERROR([$pvm_incdir/pvm3.h not found])
      fi
    else
      if test -n "$pvm_dir"; then
        pvm_incdir="$pvm_dir/include"
        AC_MSG_CHECKING([for PVM header directory])
        AC_MSG_RESULT([$pvm_incdir])
      fi
    fi

    if test -n "$pvm_libdir"; then
      AC_MSG_CHECKING([for PVM library directory])
      AC_MSG_RESULT([$pvm_libdir])
      if test -d "$pvm_libdir"; then :; else
        AC_MSG_ERROR([$pvm_libdir not found])
      fi
    else
      if test -n "$pvm_dir"; then
        if test -n "$pvm_arch"; then
          if test -d "$pvm_dir/lib/$pvm_arch"; then 
            pvm_libdir="$pvm_dir/lib/$pvm_arch"
          fi
        else
          if test -d "$pvm_dir/lib"; then 
            pvm_libdir="$pvm_dir/lib"
          fi
        fi
        if test -n "$pvm_libdir"; then
          AC_MSG_CHECKING([for PVM library directory])
          AC_MSG_RESULT([$pvm_libdir])
        fi
      fi
    fi

    if test -n "$pvm_incdir"; then
      PVM_CPPFLAGS="-I$pvm_incdir"
    fi
    if test -n "$pvm_libdir"; then
      PVM_LDFLAGS="-L$pvm_libdir"
    fi
    
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    ac_save_CPPFLAGS=$CPPFLAGS
    ac_save_LDFLAGS=$LDFLAGS
    ac_save_LIBS=$LIBS
    CPPFLAGS="$PVM_CPPFLAGS $CPPFLAGS"
    LDFLAGS="$PVM_LDFLAGS $LDFLAGS"

    AC_CHECK_HEADER([pvm3.h],,
                    [
                    if test "$pvm" = yes; then
                      AC_MSG_ERROR([check for PVM library failed])
                    else
                      pvm=no
                    fi
                    ])

    if test "$pvm" != no; then
      found=no

      if test -n "$pvm_libs"; then
        PVM_LIBS="$pvm_libs"
        LIBS="$PVM_LIBS $ac_save_LIBS"
        AC_MSG_CHECKING([for pvm_barrier() in $PVM_LIBS])
        AC_TRY_LINK([#include <pvm3.h>],[pvm_barrier(0,0);],
                    [AC_MSG_RESULT(yes); found=yes],
                    [AC_MSG_RESULT(no);
                     AC_MSG_ERROR([check for PVM library failed])])
      fi

      if test "$found" = no; then
        AC_CHECK_LIB(pvm3, pvm_spawn,
        [
        LIBS="-lpvm3 $ac_save_LIBS"
        AC_CHECK_LIB(pvm3, pvm_barrier,
        [PVM_LIBS="-lpvm3"; found=yes],
        [
        AC_CHECK_LIB(gpvm3, pvm_barrier,
        [PVM_LIBS="-lgpvm3 -lpvm3"; found=yes])
        ])
        ])
      fi

      if test "$found" = no; then
        if test "$pvm" = yes; then
          AC_MSG_ERROR([check for PVM library failed])
        fi
        pvm=no
      fi

      CPPFLAGS=$ac_save_CPPFLAGS
      LDFLAGS=$ac_save_LDFLAGS
      LIBS=$ac_save_LIBS
      AC_LANG_RESTORE
    fi
  fi
  
  if test "$pvm" != no; then
    ac_cv_have_pvm=yes
    ac_cv_pvm_cppflags="$PVM_CPPFLAGS"
    ac_cv_pvm_ldflags="$PVM_LDFLAGS"
    ac_cv_pvm_libs="$PVM_LIBS"
    AC_DEFINE(HAVE_PVM)
    OBJ_PVM='$(OBJ_PVM)'
    LIB_PVM='$(LIB_PVM)'
    EXAMPLE_PVM='$(EXAMPLE_PVM)'
    AC_MSG_NOTICE([enabling PVM support])
  else
    ac_cv_have_pvm=no
    ac_cv_pvm_cppflags=
    ac_cv_pvm_ldflags=
    ac_cv_pvm_libs=
    PVM_CPPFLAGS=
    PVM_LDFLAGS=
    PVM_LIBS=
    AC_MSG_NOTICE([disabling PVM support])
  fi
  ]
)
