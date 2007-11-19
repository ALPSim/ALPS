dnl check for ALPS Library (full version)
AC_DEFUN([AC_ALPS],
  [
  AC_MSG_CHECKING([for ALPS Libraries])
  if test -f "$ac_cv_alps_includedir/alps/config.h"; then
    AC_MSG_RESULT([OK])
  else
    AC_MSG_RESULT([not found])
    AC_MSG_ERROR([ALPS Libraries not found.])
  fi
  AC_MSG_CHECKING([for ALPS version])
  AC_MSG_RESULT([$ac_cv_alps_version])
  if test "$ac_cv_alps_type" = "full"; then :; else
    AC_MSG_NOTICE([detected ALPS Light Libraries])
    AC_MSG_ERROR([The full version of ALPS Libraries is required.])
  fi
  CXXFLAGS="$ac_cv_alps_cxxflags"
  CPPFLAGS="$ac_cv_alps_cppflags"
  LDFLAGS="$ac_cv_alps_ldflags"
  LIBS=
  AC_SUBST([LDADD], "$ac_cv_alps_libs")
  if test "$ac_cv_alps_have_mpi" = yes; then
    AC_SUBST([LDFLAGS_MPI], "$ac_cv_alps_ldflags_mpi")
    AC_SUBST([LIBS_MPI], "$ac_cv_alps_libs_mpi")
  fi

  dnl check for MPI library
  AC_MSG_CHECKING([for MPI library])
  AC_MSG_RESULT([$ac_cv_alps_have_mpi])
  AM_CONDITIONAL(HAVE_MPI, test "$ac_cv_alps_have_mpi" = yes)

  dnl check for LAPACK
  AC_MSG_CHECKING([for LAPACK library])
  AC_MSG_RESULT([$ac_cv_have_lapack])
  AM_CONDITIONAL(HAVE_LAPACK, test "$ac_cv_alps_have_lapack" = yes)
  ]
) 

dnl check for ALPS Library (light version)
AC_DEFUN([AC_ALPS_LIGHT],
  [
  AC_MSG_CHECKING([for ALPS Libraries])
  if test -f "$ac_cv_alps_includedir/alps/config.h"; then
    AC_MSG_RESULT([OK])
  else
    AC_MSG_RESULT([not found])
    AC_MSG_ERROR([ALPS Libraries not found.])
  fi
  AC_MSG_CHECKING([for ALPS version])
  AC_MSG_RESULT([$ac_cv_alps_version])
  CXXFLAGS="$ac_cv_alps_cxxflags"
  CPPFLAGS="$ac_cv_alps_cppflags"
  LDFLAGS="$ac_cv_alps_ldflags"
  LIBS=
  AC_SUBST([LDADD], "$ac_cv_alps_libs")
  if test "$ac_cv_alps_have_mpi" = yes; then
    AC_SUBST([LDFLAGS_MPI], "$ac_cv_alps_ldflags_mpi")
    AC_SUBST([LIBS_MPI], "$ac_cv_alps_libs_mpi")
  fi

  dnl check for MPI library
  AC_MSG_CHECKING([for MPI library])
  AC_MSG_RESULT([$ac_cv_alps_have_mpi])
  AM_CONDITIONAL(HAVE_MPI, test "$ac_cv_alps_have_mpi" = yes)

  dnl check for LAPACK
  AC_MSG_CHECKING([for LAPACK library])
  AC_MSG_RESULT([$ac_cv_have_lapack])
  AM_CONDITIONAL(HAVE_LAPACK, test "$ac_cv_alps_have_lapack" = yes)
  ]
) 
