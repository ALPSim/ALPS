AC_DEFUN([AC_BOOST],
  [
  AC_SUBST(BOOST_CPPFLAGS)
  AC_SUBST(BOOST_LDFLAGS)
  AC_SUBST(BOOST_LIBS)
  AC_SUBST(BOOST_INCDIR)
  AC_SUBST(BOOST_DIR)
  AC_SUBST(BOOST_SRCDIR)

  dnl check for the location of Boost

  AC_ARG_WITH(boost,
    AC_HELP_STRING([--with-boost=DIR],[Boost main tree]),
    [
    if test "x$withval" != "x"; then
      ac_cv_use_precompiled_boost=no
      boost_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      # Be sure to have absolute paths.
      case $boost_dir in
        [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
        *)  AC_MSG_ERROR([expected an absolute directory name for --with-boost : $boost_dir]);;
      esac
      if test -f "$boost_dir/boost/config.hpp"; then :; else
        AC_MSG_ERROR([Boost main tree not found])
      fi
      if test -f "$boost_dir/libs/program_options/src/cmdline.cpp"; then :; else
        AC_MSG_ERROR([Boost main tree not found])
      fi
      boost_incdir="$boost_dir"
      boost_srcdir="$boost_dir/libs"
    fi
    ]
  )
  AC_ARG_WITH(boost-incdir,
    AC_HELP_STRING([--with-boost-incdir=DIR],[Boost include directory (only for precompiled Boost library)]),
    [
    if test "x$withval" != "x"; then
      ac_cv_use_precompiled_boost=yes
      boost_incdir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      # Be sure to have absolute paths.
      case $boost_incdir in
        [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
        *)  AC_MSG_ERROR([expected an absolute directory name for --with-boost-incdir : $boost_incdir]);;
      esac
      if test -f "$boost_incdir/boost/config.hpp"; then :; else
        AC_MSG_ERROR([Boost include files not found])
      fi
    fi
    ]
  )
  AC_ARG_WITH(boost-libdir,
    AC_HELP_STRING([--with-boost-libdir=DIR],[Boost library directory (only for precompiled Boost library)]),
    [
    if test "x$withval" != "x"; then
      ac_cv_use_precompiled_boost=yes
      boost_libdir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      # Be sure to have absolute paths.
      case $boost_libdir in
        [[\\/$]]* | ?:[[\\/]]* | NONE | '' ) ;;
        *)  AC_MSG_ERROR([expected an absolute directory name for --with-boost-libdir : $boost_libdir]);;
      esac
    fi
    ]
  )
  AC_ARG_WITH(boost-libs,
    AC_HELP_STRING([--with-boost-libs=LIBS],[Boost libraries (only for precompiled Boost library)]),
    [
    if test "x$withval" != "x"; then
      ac_cv_use_precompiled_boost=yes
      boost_libs=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
    fi
    ]
  )
  AC_ARG_WITH(boost-toolset,
    AC_HELP_STRING([--with-boost-toolset=TOOLSET],[Boost toolset abbreviation, e.g. gcc-d for gcc with debugging, il-mt for intel-linux with multi-threading etc (default none) (only for precompiled Boost library)]),
    [
    if test "x$withval" != "x"; then
      ac_use_precompiled_boost=yes
      boost_toolset=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
    fi
    ]
  )
  AC_ARG_WITH(boost-all,
    AC_HELP_STRING([--with-boost-all],[compile all optional Boost libraries]), [
      if test "x$withval" = "xno"; then
        ac_cv_boost_signals=no
        ac_cv_boost_thread=no
      else
        ac_cv_boost_signals=yes
        ac_cv_boost_thread=yes
      fi
    ]
  )
  AC_ARG_WITH(boost-mpi,
    AC_HELP_STRING([--without-boost-mpi],[not compile Boost.MPI library]), [
      if test "x$withval" = "xno"; then
        ac_cv_boost_mpi=no
      else
        ac_cv_boost_mpi=yes
      fi
    ]
  )
  AC_ARG_WITH(boost-signals,
    AC_HELP_STRING([--with-boost-signals],[compile Boost.Signals library]), [
      if test "x$withval" = "xno"; then
        ac_cv_boost_signals=no
      else
        ac_cv_boost_signals=yes
      fi
    ]
  )
  AC_ARG_WITH(boost-thread,
    AC_HELP_STRING([--with-boost-thread],[compile Boost.Thread library]), [
      if test "x$withval" = "xno"; then
        ac_cv_boost_thread=no
      else
        ac_cv_boost_thread=yes
      fi
    ]
  )
  AC_ARG_WITH(boost-wchar,
    AC_HELP_STRING([--with-boost-wchar],[compile with wchar support in Boost.Regex and Boost.Serialization libraries]), [
      if test "x$withval" = "xno"; then
        ac_cv_boost_wchar=no
      else
        ac_cv_boost_wchar=yes
      fi
    ]
  )

  boost_dir_s=

  for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local/src; do
    for b in boost boost_1_41_0 boost_1_40_0 boost_1_39_0 boost_1_38_0 boost_1_37_0 boost_1_36_0 boost_1_35_0 boost_1_34_1 boost_1_34_0 boost_1_33_1 boost_1_33_0 boost_1_32_0; do
      if test -f "$d/$b/boost/config.hpp"; then
        if test -f "$d/$b/libs/filesystem/src/exception.cpp" || test -f "$d/$b/libs/filesystem/src/operations.cpp"; then
          boost_dir_s="$d/$b"
          boost_incdir_s="$boost_dir_s"
          boost_srcdir_s="$boost_dir_s/libs"
          break
        fi
      fi
    done
    if test -n "$boost_dir_s"; then
      break
    fi 
  done

  dnl check for Boost include files

  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  ac_save_CPPFLAGS=$CPPFLAGS
  ac_save_LDFLAGS=$LDFLAGS
  ac_save_LIBS=$LIBS

  AC_MSG_CHECKING([for Boost include files])
  found=no

  if test -n "$boost_incdir"; then
    ac_cv_boost_incdir="$boost_incdir"
    BOOST_CPPFLAGS="-I$boost_incdir"
    BOOST_INCDIR=$boost_incdir
    CPPFLAGS="$BOOST_CPPFLAGS $ac_save_CPPFLAGS"
  fi
  AC_TRY_COMPILE([#include <boost/config.hpp>],,found=yes)

  if test "$found" = no; then
    if test -n "$boost_incdir_s"; then
      ac_cv_boost_incdir="$boost_incdir_s"
      BOOST_CPPFLAGS="-I$boost_incdir_s"
      CPPFLAGS="$BOOST_CPPFLAGS $ac_save_CPPFLAGS"
      AC_TRY_COMPILE([#include <boost/config.hpp>],,found=yes)
      BOOST_INCDIR=$boost_incdir_s
    fi
  fi

  if test "$found" = yes; then
    if test -n "$BOOST_CPPFLAGS"; then
      AC_MSG_RESULT([$BOOST_CPPFLAGS])
    else
      AC_MSG_RESULT([yes])
    fi
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Boost include files not found])
  fi

dnl Boost version

  AC_MSG_CHECKING([for Boost version])
  AC_TRY_COMPILE([#include <boost/version.hpp>
#if BOOST_VERSION < 104200
#error
#endif],,ac_cv_boost_version=svn)
  if test -z "$ac_cv_boost_version"; then
    AC_TRY_COMPILE([#include <boost/version.hpp>
#if BOOST_VERSION < 104100
#error
#endif],,ac_cv_boost_version=1_41)
  fi
  if test -z "$ac_cv_boost_version"; then
    AC_TRY_COMPILE([#include <boost/version.hpp>
#if BOOST_VERSION < 104000
#error
#endif],,ac_cv_boost_version=1_40)
  fi
  case "x$ac_cv_boost_version" in
    xsvn )
      AC_MSG_RESULT([SVN])
      ;;
    x1_41 )
      AC_MSG_RESULT([1.41])
      ;;
    * )
      AC_MSG_RESULT([unknown])
      AC_MSG_ERROR([Boost library is too old])
      ;;
  esac
  AM_CONDITIONAL(BOOST_SVN, test "$ac_cv_boost_version" = svn)
  AM_CONDITIONAL(BOOST_1_41, test "$ac_cv_boost_version" = 1_41)

  dnl pre-compiled boost library

  if test "$ac_cv_use_precompiled_boost" = no; then :; else

    AC_MSG_CHECKING([for Boost toolset abbreviation])
    if test -n "$boost_toolset"; then
      AC_MSG_RESULT([$boost_toolset])
    else
      AC_MSG_RESULT([none])
    fi

    AC_MSG_CHECKING([for pre-compiled Boost library])
    found=no

    if test -n "$boost_libdir"; then
      BOOST_LDFLAGS="-L$boost_libdir"
      LDFLAGS="$BOOST_LDFLAGS $ac_save_LDFLAGS"
    fi
    if test -n "$boost_libs"; then
      BOOST_LIBS="$boost_libs"
      LIBS="$BOOST_LIBS $ac_save_LIBS"
    else
      if test -z "$boost_toolset"; then
        BOOST_LIBS="-lboost_date_time -lboost_filesystem -lboost_program_options -lboost_regex -lboost_serialization -lboost_system"
      else
        BOOST_LIBS="-lboost_date_time-$boost_toolset -lboost_filesystem-$boost_toolset -lboost_program_options-$boost_toolset -lboost_regex-$boost_toolset -lboost_serialization-$boost_toolset -lboost_system-$boost_toolset"
      fi
      LIBS="$BOOST_LIBS $ac_save_LIBS"
    fi

    AC_TRY_LINK([#include <boost/date_time/date_generators.hpp>],[const char *ptr = boost::date_time::nth_as_str(0);],found=yes)

    if test "$found" = yes; then
      found=no
      AC_TRY_LINK([#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>],[boost::filesystem::path p(boost::filesystem::initial_path());],
      found=yes)
    fi

    if test "$found" = yes; then
      found=no
      AC_TRY_LINK([#include <boost/regex.hpp>],[boost::regbase re;],
      found=yes)
    fi

    if test "$found" = yes; then
      AC_MSG_RESULT([$BOOST_LDFLAGS $BOOST_LIBS])
      ac_cv_use_precompiled_boost=yes
    else
      BOOST_LDFLAGS=
      BOOST_LIBS=
      AC_MSG_RESULT([no])
      if test "$ac_cv_use_precompiled_boost" = yes; then
        AC_MSG_ERROR([Boost library not found])
      fi
    fi
    
  fi

  dnl boost sources

  if test "$ac_cv_use_precompiled_boost" = yes; then :; else
    AC_MSG_CHECKING([for Boost source files])
    ac_cv_use_precompiled_boost=no
    if test -n "$boost_srcdir"; then
      AC_MSG_RESULT([$boost_srcdir])
      BOOST_DIR="$boost_dir"
      BOOST_SRCDIR="$boost_srcdir"
    else
      if test -n "$boost_srcdir_s"; then
        AC_MSG_RESULT([$boost_srcdir_s])
        BOOST_DIR="$boost_dir_s"
        BOOST_SRCDIR="$boost_srcdir_s"
      else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Boost source files not found])
      fi
    fi
  fi

  ac_cv_boost_cppflags="$BOOST_CPPFLAGS"
  ac_cv_boost_ldflags="$BOOST_LDFLAGS"
  ac_cv_boost_libs="$BOOST_LIBS"
  ac_cv_boost_srcdir="$BOOST_SRCDIR"

  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
  LIBS="$ac_save_LIBS"
  AC_LANG_RESTORE
  ]
)

AC_DEFUN([AC_BOOST_LIBS],
  [
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  ac_save_CPPFLAGS=$CPPFLAGS
  ac_save_LDFLAGS=$LDFLAGS
  ac_save_LIBS=$LIBS

  #
  # boost object libraries
  # 

  AM_CONDITIONAL(BUILD_BOOST, test "$ac_cv_use_precompiled_boost" != yes)

  if test "$ac_cv_use_precompiled_boost" = yes; then

    AC_MSG_CHECKING([for Boost.MPI library])
    ac_cv_boost_mpi=no
    if test -n "$boost_libs"; then
      BOOST_MPI_LIBS=
    else
      if test -z "$boost_toolset"; then
        BOOST_MPI_LIBS="-lboost_mpi"
      else
        BOOST_MPI_LIBS="-lboost_mpi-$boost_toolset"
      fi
    fi
    CPPFLAGS="$ac_save_CPPFLAGS $BOOST_CPPFLAGS $MPI_CPPFLAGS"
    LDFLAGS="$ac_save_LDFLAGS $BOOST_LDFLAGS $MPI_LDFLAGS"
    LIBS="$BOOST_MPI_LIBS $BOOST_LIBS $MPI_LIBS $ac_save_LIBS"
    AC_TRY_LINK([#include <boost/mpi/communicator.hpp>],[boost::mpi::communicator comm;],ac_cv_boost_mpi=yes)
    if test "$ac_cv_boost_mpi" = yes; then
      if test -n "$BOOST_MPI_LIBS"; then
        AC_MSG_RESULT([$BOOST_MPI_LIBS])
        BOOST_LIBS="$BOOST_MPI_LIBS $BOOST_LIBS"
      else
        AC_MSG_RESULT([yes])
      fi
    else
      AC_MSG_RESULT([no])
    fi

    AC_MSG_CHECKING([for Boost.Signals library])
    ac_cv_boost_signals=no
    if test -n "$boost_libs"; then
      BOOST_SIGNALS_LIBS=
    else
      if test -z "$boost_toolset"; then
        BOOST_SIGNALS_LIBS="-lboost_signals"
      else
        BOOST_SIGNALS_LIBS="-lboost_signals-$boost_toolset"
      fi
    fi
    LIBS="$BOOST_LIBS $BOOST_SIGNALS_LIBS $ac_save_LIBS"
    AC_TRY_LINK([#include <boost/signals/signal2.hpp>
#include <string>],[boost::signal2<void, int, int, boost::last_value<void>, std::string> sig;],ac_cv_boost_signals=yes)
    if test "$ac_cv_boost_signals" = yes; then
      if test -n "$BOOST_SIGNALS_LIBS"; then
        AC_MSG_RESULT([$BOOST_SIGNALS_LIBS])
        BOOST_LIBS="$BOOST_LIBS $BOOST_SIGNALS_LIBS"
      else
        AC_MSG_RESULT([yes])
      fi
    else
      AC_MSG_RESULT([no])
    fi

    AC_MSG_CHECKING([for Boost.Thread library])
    ac_cv_boost_thread=no
    if test -n "$boost_libs"; then
      BOOST_THREAD_LIBS=
    else
      if test -z "$boost_toolset"; then
        BOOST_THREAD_LIBS="-lboost_thread"
      else
        BOOST_THREAD_LIBS="-lboost_thread-$boost_toolset"
      fi
    fi
    LIBS="$BOOST_LIBS $BOOST_THREAD_LIBS $ac_save_LIBS"
    AC_TRY_LINK([#include <boost/thread/xtime.hpp>],[int x = boost::xtime_get(0, 0);],ac_cv_boost_thread=yes)
    if test "$ac_cv_boost_thread" = yes; then
      if test -n "$BOOST_THREAD_LIBS"; then
        AC_MSG_RESULT([$BOOST_THREAD_LIBS])
        BOOST_LIBS="$BOOST_LIBS $BOOST_THREAD_LIBS"
      else
        AC_MSG_RESULT([yes])
      fi
    else
      AC_MSG_RESULT([no])
    fi
    AM_CONDITIONAL(BUILD_BOOST_MPI, test "$ac_cv_boost_mpi" = yes)
    AM_CONDITIONAL(BUILD_BOOST_SIGNALS, false)
    AM_CONDITIONAL(BUILD_BOOST_THREAD, test "$ac_cv_boost_thread" = yes)
    AM_CONDITIONAL(BUILD_BOOST_WCHAR, test "$ac_cv_boost_wchar" = yes)
  else
    # check whether Boost.MPI is compiled
    AC_MSG_CHECKING([whether to build Boost.MPI])
    test -z "$ac_cv_boost_mpi" && ac_cv_boost_mpi=yes
    test "$ac_cv_have_mpi" = yes || ac_cv_boost_mpi=no
    if test "$ac_cv_boost_version" = svn || test "$ac_cv_boost_version" = 1_41 || test "$ac_cv_boost_version" = 1_40 || test "$ac_cv_boost_version" = 1_39 || test "$ac_cv_boost_version" = 1_38 || test "$ac_cv_boost_version" = 1_37 || test "$ac_cv_boost_version" = 1_36 || test "$ac_cv_boost_version" = 1_35; then
      :;
    else
      ac_cv_boost_mpi=no
    fi
    AC_MSG_RESULT([$ac_cv_boost_mpi])
    AM_CONDITIONAL(BUILD_BOOST_MPI, test "$ac_cv_boost_mpi" = yes)

    # check whether Boost.Signals is compiled
    AC_MSG_CHECKING([whether to build Boost.Signals])
    test -z "$ac_cv_boost_signals" && ac_cv_boost_signals=no
    if test "$ac_cv_boost_signals" = yes; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
    fi
    AM_CONDITIONAL(BUILD_BOOST_SIGNALS, test "$ac_cv_boost_signals" = yes)

    # check whether Boost.Thread is compiled
    AC_MSG_CHECKING([whether to build Boost.Thread])
    test -z "$ac_cv_boost_thread" && ac_cv_boost_thread=no
    if test "$ac_cv_have_pthread" = yes; then
      if test "$ac_cv_boost_thread" = yes; then
        AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
      fi
    else
      AC_MSG_RESULT([no])
    fi
    AM_CONDITIONAL(BUILD_BOOST_THREAD, test "$ac_cv_boost_thread" = yes)

    # check whether enable wchar support in Boost.Regex is compiled
    AC_MSG_CHECKING([whether to enable wchar support in Boost.Regex and Boost.Serialization])
    test -z "$ac_cv_boost_wchar" && ac_cv_boost_wchar=no
    AC_MSG_RESULT([$ac_cv_boost_wchar])
    AM_CONDITIONAL(BUILD_BOOST_WCHAR, test "$ac_cv_boost_wchar" = yes)
  fi

  AM_CONDITIONAL(HAVE_BOOST_MPI, test "$ac_cv_boost_mpi" = yes)
  test "$ac_cv_boost_mpi" = yes && AC_DEFINE(ALPS_HAVE_BOOST_MPI)
  AM_CONDITIONAL(HAVE_BOOST_SIGNALS, test "$ac_cv_boost_signals" = yes)
  test "$ac_cv_boost_signals" = yes && AC_DEFINE(ALPS_HAVE_BOOST_SIGNALS)
  AM_CONDITIONAL(HAVE_BOOST_THREAD, test "$ac_cv_boost_thread" = yes)
  test "$ac_cv_boost_thread" = yes && AC_DEFINE(ALPS_HAVE_BOOST_THREAD)
  AM_CONDITIONAL(HAVE_BOOST_WCHAR, test "$ac_cv_boost_wchar" = yes)
  test "$ac_cv_boost_wchar" = yes && AC_DEFINE(ALPS_HAVE_BOOST_WCHAR)

  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
  LIBS="$ac_save_LIBS"
  AC_LANG_RESTORE
  ]
)

AC_DEFUN([AC_BOOST_CONFIGURE],
  [
  # check whether Boost.Configure is used or not
  AC_ARG_ENABLE(boost-config,
    AC_HELP_STRING([--enable-boost-config],
      [setup Boost configuration header file]),
    [
    if test "x$enableval" != "xno"; then
      ac_cv_boost_config=yes
    fi
    ]
  )
  AC_SUBST(BOOST_USER_CONFIG_H)
  BOOST_USER_CONFIG_H=
  if test "$ac_cv_use_precompiled_boost" = yes; then
    ac_cv_boost_config=no
  else
    AC_MSG_CHECKING([whether to run Boost configure script])
    test -z "$ac_cv_boost_config" && ac_cv_boost_config=no
    if test "$ac_cv_boost_config" = yes; then
      ac_cv_boost_config="$includedir/$PACKAGE_NAME/boost-user.hpp"
      if test -f "$BOOST_SRCDIR/config/configure"; then :; else
        AC_MSG_RESULT
        AC_MSG_ERROR([Boost configure script was not found.])
      fi
    fi
    AC_MSG_RESULT([$ac_cv_boost_config])
  fi
  
  # configure boost if needed

  if test "$ac_cv_boost_config" != no; then
    AC_MSG_NOTICE([running Boost configuration script])
    boost_configure="sh $BOOST_SRCDIR/config/configure"
    mkdir -p src/boost
    command="(cd src/boost && export CXX=\"$CXX\" CXXFLAGS=\"$CXXFLAGS\" CPPFLAGS=\"$CPPFLAGS\" LDFLAGS=\"$LDFLAGS\" && $boost_configure --with-boost=$boost_incdir)"
    echo "$command"
    eval "$command"
    CPPFLAGS="$CPPFLAGS -DBOOST_USER_CONFIG=\<boost/user.hpp\>"
    BOOST_USER_CONFIG_H="boost/user.hpp"
  fi
  ]
)
