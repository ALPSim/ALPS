AC_DEFUN([AC_PTHREAD],
  [
  AC_SUBST(OBJ_PTHREAD)
  AC_SUBST(LIB_PTHREAD)
  AC_SUBST(EXAMPLE_PTHREAD)

  # check for configure option
  pthread=no
  AC_ARG_WITH(pthread,
    AC_HELP_STRING([--with-pthread=DIR],[use pthread library]),
    [
    if test "x$withval" != "xno"; then
      pthread=yes
    fi
    ]
  )

  if test "$pthread" = yes; then
    AC_CHECK_LIB(pthread, pthread_create, ,
      [AC_CHECK_FUNC(pthread_create, ,
        [AC_CHECK_LIB(pthread, __pthread_create)] # for Tru64
      )]
    )
    if test "$ac_cv_lib_pthread_pthread_create" = yes; then
      ac_cv_have_pthread=yes
    else
      if test "$ac_cv_func_pthread_create" = yes; then
        ac_cv_have_pthread=yes
      else
        if test "$ac_cv_lib_pthread___pthread_create" = yes; then
          ac_cv_have_pthread=yes
        fi
      fi
    fi
  fi

  if test "$ac_cv_have_pthread" = yes; then 
    OBJ_PTHREAD='$(OBJ_PTHREAD)'
    LIB_PTHREAD='$(LIB_PTHREAD)'
    EXAMPLE_PTHREAD='$(EXAMPLE_PTHREAD)'
    AC_MSG_NOTICE([enabling pthread support])
  else
    ac_cv_have_pthread=no
    AC_MSG_NOTICE([disabling pthread support])
  fi
  ]
)
