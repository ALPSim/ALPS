dnl configuration script for ALPS
dnl Process this file with autoconf to produce a configure script.

dnl disable loading config.site
define([AC_SITE_LOAD], )dnl
AC_PREREQ(2.59)

dnl
dnl Version Section ---------------------------------------------------------
dnl

dnl Version of ALPS Library

define([alps_version], [1.3])

dnl Latest publication year of ALPS

define([alps_year], [2008])

dnl Increment lib_revision if source code has changed at all.
dnl
dnl If the API has changed, increment lib_current and set lib_revision to 0.
dnl
dnl If the API changes compatibly (i.e. simply adding a new function
dnl without changing or removing earlier interfaces), then increment lib_age.
dnl 
dnl If the API changes incompatibly set lib_age back to 0
dnl
dnl cf. http://www.gnu.org/software/libtool/manual.html#SEC32

define([lib_current],  [1])
define([lib_revision], [5])
define([lib_age],      [1])

dnl
dnl Initialization ---------------------------------------------------------
dnl

AC_INIT([alps], alps_version)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT([$HOME/ALPS])
AM_MAINTAINER_MODE

dnl Full or Light version
if test -f "$srcdir/src/alps/lattice.h"; then
  ac_cv_alps_type="full"
  AC_MSG_NOTICE([starting configuration for ALPS Libraries version $PACKAGE_VERSION])
else
  ac_cv_alps_type="light"
  AC_MSG_NOTICE([starting configuration for ALPS Light Libraries version $PACKAGE_VERSION])
fi	

dnl version
ac_cv_alps_version=alps_version
AC_DEFINE_UNQUOTED([ALPS_VERSION], ["alps_version"])
if test "$ac_cv_alps_type" = full; then
  AC_DEFINE_UNQUOTED([ALPS_VERSION_STRING], ["ALPS Libraries version alps_version"])
else
  AC_DEFINE_UNQUOTED([ALPS_VERSION_STRING], ["ALPS Light Libraries version alps_version"])
fi
undefine([alps_version])

dnl latest publish year of ALPS
ac_cv_alps_year=alps_year
AC_DEFINE_UNQUOTED([ALPS_YEAR], ["alps_year"])
undefine([alps_year])

dnl shared library version
AC_SUBST(LIB_VERSION_INFO)
LIB_VERSION_INFO="-version-info lib_current:lib_revision:lib_age"
undefine([lib_current])
undefine([lib_revision])
undefine([lib_age])

dnl setting up cache file
cache_file='config/config.cache'
if test -d "config"; then :; else mkdir config; fi
if test -f $cache_file; then :; else touch $cache_file; fi

dnl
dnl Directories ------------------------------------------------------------
dnl

test "$prefix" = NONE && prefix=$ac_default_prefix
test "$exec_prefix" = NONE && exec_prefix=$prefix
test "$bindir" = '${exec_prefix}/bin' && bindir="$exec_prefix/bin"
test "$datadir" = '${prefix}/share' && datadir="$prefix/share"
test "$libdir" = '${exec_prefix}/lib' && libdir="$exec_prefix/lib"
test "$includedir" = '${prefix}/include' && includedir="$prefix/include"
AC_SUBST(docdir)
docdir="$prefix/doc"
AC_SUBST(SRC_DIR)
SRC_DIR="$srcdir"

ac_cv_alps_prefix="$prefix"
ac_cv_alps_exec_prefix="$exec_prefix"
ac_cv_alps_bindir="$bindir"
ac_cv_alps_datadir="$datadir"
ac_cv_alps_libdir="$libdir"
ac_cv_alps_xmldir="$libdir/xml"
ac_cv_alps_includedir="$includedir"
ac_cv_alps_docdir="$docdir"
ac_cv_alps_include_mk="$datadir/$PACKAGE_NAME/include.mk"
ac_cv_alps_srcdir="$srcdir"

AC_DEFINE_UNQUOTED([ALPS_PREFIX], ["$prefix"])
AC_DEFINE_UNQUOTED([ALPS_SRCDIR], ["$srcdir"])

dnl XSLT path
AC_DEFINE_UNQUOTED([ALPS_XML_DIR], ["$ac_cv_alps_xmldir"])
ALPS_XML_DIR=$ac_cv_alps_xmldir
AC_SUBST(ALPS_XML_DIR)

test -z "$ac_hostname" && ac_hostname="unknown"
AC_DEFINE_UNQUOTED([ALPS_CONFIG_HOST], ["$ac_hostname"])
ac_username=$USER
test -z "$ac_username" && ac_username="unknown"
AC_DEFINE_UNQUOTED([ALPS_CONFIG_USER], ["$ac_username"])

dnl
dnl Compilers --------------------------------------------------------
dnl

AC_CANONICAL_HOST

dnl C/C++ compiler
AC_COMPILER
ac_cv_alps_compiler="$ac_cv_compiler"
ac_cv_alps_compiler_optimization="$ac_cv_compiler_optimization"
ac_cv_alps_compiler_exceptions="$ac_cv_compiler_exceptions"
ac_cv_alps_compiler_warnings="$ac_cv_compiler_warnings"
ac_cv_alps_cc="$ac_cv_compiler_cc"
ac_cv_alps_cflags="$ac_cv_compiler_cflags"
ac_cv_alps_cxx="$ac_cv_compiler_cxx"
ac_cv_alps_cxxflags="$ac_cv_compiler_cxxflags"

dnl C/C++ preprocessor
AC_PROG_CPP
AC_PROG_CXXCPP
ac_cv_alps_cpp="$CPP"
ac_cv_alps_cxxcpp="$CXXCPP"

dnl disable F77
test -z "$F77" && F77=no

dnl
dnl System Headers ---------------------------------------------------------
dnl

AC_LANG_CPLUSPLUS

AC_CHECK_HEADER([inttypes.h])
test "$ac_cv_header_inttypes_h" = yes && AC_DEFINE(ALPS_HAVE_INTTYPES_H, [], [inttypes.h])

AC_CHECK_HEADER([pwd.h])
test "$ac_cv_header_pwd_h" = yes && AC_DEFINE(ALPS_HAVE_PWD_H, [], [pwd.h])

AC_CHECK_HEADER([stdarg.h])
test "$ac_cv_header_stdarg_h" = yes && AC_DEFINE(ALPS_HAVE_STDARG_H, [], [stdarg.h])

AC_CHECK_HEADER([stdint.h])
test "$ac_cv_header_stdint_h" = yes && AC_DEFINE(ALPS_HAVE_STDINT_H, [], [stdint.h])

AC_CHECK_HEADER([unistd.h])
test "$ac_cv_header_unistd_h" = yes && AC_DEFINE(ALPS_HAVE_UNISTD_H, [], [unistd.h])

AC_CHECK_HEADER([bind/bitypes.h])
test "$ac_cv_header_bind_bitypes_h" = yes && AC_DEFINE(ALPS_HAVE_BIND_BITYPES_H, [], [bind/bitypes.h])

AC_CHECK_HEADER([rpc/rpc.h])
test "$ac_cv_header_rpc_rpc_h" = yes && AC_DEFINE(ALPS_HAVE_RPC_XDR_H, [], [rpc/rpc.h])

AC_CHECK_HEADER([sys/int_types.h])
test "$ac_cv_header_sys_int_types_h" = yes && AC_DEFINE(ALPS_HAVE_SYS_INT_TYPES_H, [], [sys/int_types.h])

AC_CHECK_HEADER([sys/systeminfo.h])
test "$ac_cv_header_sys_systeminfo_h" = yes && AC_DEFINE(ALPS_HAVE_SYS_SYSTEMINFO_H, [], [sys/systeminfo.h])

AC_CHECK_HEADER([sys/time.h])
test "$ac_cv_header_sys_time_h" = yes && AC_DEFINE(ALPS_HAVE_SYS_TIME_H, [], [sys/time.h])

AC_CHECK_HEADER([sys/types.h])
test "$ac_cv_header_sys_types_h" = yes && AC_DEFINE(ALPS_HAVE_SYS_TYPES_H, [], [sys/types.h])

dnl
dnl Other tools --------------------------------------------------------
dnl

dnl make --disable-sahred the default for IBM compiler
if test x"$ac_cv_compiler" = xibm32 || test x"$ac_cv_compiler" = xibm64; then
  test -z "${enable_shared+set}" && enable_shared=no
fi

dnl set 64-bit flag for IBM compiler
dnl if test x"$ac_cv_compiler" = xibm64; then
dnl   test -z "$AR" && AR="ar -X64"
dnl   test -z "$NM" && NM="nm -X64 -B"
dnl fi

dnl libtool
AC_PROG_LIBTOOL
ac_cv_enable_shared="$enable_shared"
ac_cv_enable_static="$enable_static"

#
# checking boost directory
#

AC_SUBST(BOOST_CPPFLAGS)
AC_SUBST(BOOST_LDFLAGS)
AC_SUBST(BOOST_LIBS)

ac_cv_use_precompiled_boost=

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
    if test -f "$boost_dir/libs/filesystem/src/exception.cpp"; then :; else
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
AC_ARG_WITH(boost-signals,
  AC_HELP_STRING([--with-boost-signals],[compile Boost::signals library]), [
    if test "x$withval" = "xno"; then
      ac_cv_boost_signals=no
    else
      ac_cv_boost_signals=yes
    fi
  ]
)
AC_ARG_WITH(boost-thread,
  AC_HELP_STRING([--with-boost-thread],[compile boost::thread library]), [
    if test "x$withval" = "xno"; then
      ac_cv_boost_thread=no
    else
      ac_cv_boost_thread=yes
    fi
  ]
)
AC_ARG_WITH(boost-wchar,
  AC_HELP_STRING([--with-boost-wchar],[compile with wchar support in boost::regex and boost::serialization libraries]), [
    if test "x$withval" = "xno"; then
      ac_cv_boost_wchar=no
    else
      ac_cv_boost_wchar=yes
    fi
  ]
)

boost_dir_s=

for d in $HOME $HOME/src $prefix $prefix/src /usr/local /usr/local/src; do
  for b in boost boost_1_34_1 boost_1_34_0 boost_1_33_1 boost_1_33_0 boost_1_32_0; do
    if test -f "$d/$b/boost/config.hpp"; then
      if test -f "$d/$b/libs/filesystem/src/exception.cpp"; then
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

AC_LANG_SAVE
AC_LANG_CPLUSPLUS
ac_save_CPPFLAGS=$CPPFLAGS
ac_save_LDFLAGS=$LDFLAGS
ac_save_LIBS=$LIBS

AC_MSG_CHECKING([for Boost include files])
found=no

ac_cv_boost_incdir=

AC_SUBST(BOOST_INCDIR)

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
#if BOOST_VERSION < 103500
#error
#endif],,ac_cv_boost_version=cvs)
if test -z "$ac_cv_boost_version"; then
  AC_TRY_COMPILE([#include <boost/version.hpp>
#if BOOST_VERSION < 103400 
#error
#endif],,ac_cv_boost_version=1_34)
fi
if test -z "$ac_cv_boost_version"; then
  AC_TRY_COMPILE([#include <boost/version.hpp>
#if BOOST_VERSION < 103300
#error
#endif],,ac_cv_boost_version=1_33)
fi
if test -z "$ac_cv_boost_version"; then
  AC_TRY_COMPILE([#include <boost/version.hpp>
#if BOOST_VERSION < 103200
#error
#endif],,ac_cv_boost_version=1_32)
fi
case "x$ac_cv_boost_version" in
  xcvs )
    AC_MSG_RESULT([CVS])
    ;;
  x1_34 )
    AC_MSG_RESULT([1.34])
    ;;
  x1_33 )
    AC_MSG_RESULT([1.33])
    ;;
  x1_32 )
    AC_MSG_RESULT([1.32])
    ;;
  * )
    AC_MSG_RESULT([unknown])
    AC_MSG_ERROR([Boost library is too old])
    ;;
esac

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
      BOOST_LIBS="-lboost_date_time -lboost_filesystem -lboost_program_options -lboost_regex -lboost_serialization"
    else
      BOOST_LIBS="-lboost_date_time-$boost_toolset -lboost_filesystem-$boost_toolset -lboost_program_options-$boost_toolset -lboost_regex-$boost_toolset -lboost_serialization-$boost_toolset"
    fi
    LIBS="$BOOST_LIBS $ac_save_LIBS"
  fi

  AC_TRY_LINK([
    #include <boost/date_time/date_generators.hpp>
  ],[const char *ptr = boost::date_time::nth_as_str(0);],
  found=yes)

  if test "$found" = yes; then
    found=no
    AC_TRY_LINK([
      #include <boost/filesystem/operations.hpp>
      #include <boost/filesystem/path.hpp>
    ],[boost::filesystem::path p(boost::filesystem::initial_path());],
    found=yes)
  fi

  if test "$found" = yes; then
    found=no
    AC_TRY_LINK([
      #include <boost/regex.hpp>
    ],[boost::regbase re;],
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

  if test "$ac_cv_use_precompiled_boost" = yes; then
    #
    # TODO: check for Boost::program_options and Boost::serialization libraries
    #

    AC_MSG_CHECKING([for Boost::signals library])
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
    AC_TRY_LINK([
      #include <boost/signals/signal2.hpp>
      #include <string>
    ],[boost::signal2<void, int, int, boost::last_value<void>, std::string> sig;],
    ac_cv_boost_signals=yes)
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

    AC_MSG_CHECKING([for Boost::thread library])
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
    AC_TRY_LINK([#include <boost/thread/xtime.hpp>],
      [int x = boost::xtime_get(0, 0);],
      ac_cv_boost_thread=yes)
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
  fi
fi

dnl boost sources

AC_SUBST(BOOST_DIR)
AC_SUBST(BOOST_SRCDIR)

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

CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"
LIBS="$ac_save_LIBS"
AC_LANG_RESTORE

ac_cv_boost_cppflags="$BOOST_CPPFLAGS"
ac_cv_boost_ldflags="$BOOST_LDFLAGS"
ac_cv_boost_libs="$BOOST_LIBS"
ac_cv_boost_srcdir="$BOOST_SRCDIR"


# checking for other C++ headers

AC_CHECK_HEADER([cstdlib],,
  [
  if test -n "$ac_cv_boost_incdir"; then
    if test -f "$ac_cv_boost_incdir/boost/compatibility/cpp_c_headers/cstdlib"; then
      CPPFLAGS="-I$ac_cv_boost_incdir/boost/compatibility/cpp_c_headers $CPPFLAGS"
    fi
    AC_MSG_NOTICE([$ac_cv_boost_incdir/boost/compatibility/cpp_c_headers was added to include path])
  else
    AC_MSG_ERROR([cstdlib not fount])
  fi
  ]
)

AC_CXX_HAVE_VALARRAY
test "$ac_cv_cxx_have_valarray" = yes && AC_DEFINE(ALPS_HAVE_VALARRAY, [], [valarray])

#
# check for libm
#

AC_CHECK_LIB(m, sqrt)

#
# check for nanosleep
#

AC_CHECK_FUNC(nanosleep, ,[AC_CHECK_LIB(rt, nanosleep)])

#
# check for pthread
#

AC_PTHREAD
test "$ac_cv_have_pthread" = yes && AC_DEFINE(ALPS_HAVE_PTHREAD, [], [pthread])

#
# check for HDF5
#

AC_HDF5
test "$ac_cv_have_hdf5" = yes && AC_DEFINE(ALPS_HAVE_HDF5, [], [hdf5])

dnl
dnl check for LAPACK
dnl

AC_LAPACK

#
# check for SQLite
#

AC_SQLITE

#
# check for XML parser
#

AC_XMLPARSER
test "$ac_cv_xml_parser" = "expat" && AC_DEFINE(ALPS_HAVE_EXPAT_PARSER, [], [expat])
test "$ac_cv_xml_parser" = "xerces" && AC_DEFINE(ALPS_HAVE_XERCES_PARSER, [], [xerces])

#
# check for MPI
#

AC_SUBST(LIB_MPI)
ac_cv_alps_have_mpi=no
AC_MPI
if test "$ac_cv_have_mpi" = yes; then
  ac_cv_alps_have_mpi=yes
  AC_DEFINE(ALPS_HAVE_MPI, [], [mpi])
  LIB_MPI='$(LIB_MPI)'
fi

#
# boost object libraries
# 

# check whether Boost::signals is compiled
if test "$ac_cv_use_precompiled_boost" = yes; then :; else
  AC_MSG_CHECKING([whether to build boost::signals])
  test -z "$ac_cv_boost_signals" && ac_cv_boost_signals=no
  if test "$ac_cv_boost_signals" = yes; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
  fi
fi

# check whether boost::thread is compiled
if test "$ac_cv_use_precompiled_boost" = yes; then :; else
  AC_MSG_CHECKING([whether to build boost::thread])
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
fi

# check whether boost::thread is compiled
if test "$ac_cv_use_precompiled_boost" = yes; then :; else
  AC_MSG_CHECKING([whether to enable wchar support in boost::regex and boost::serialization])
  test -z "$ac_cv_boost_wchar" && ac_cv_boost_wchar=no
  AC_MSG_RESULT([$ac_cv_boost_wchar])
fi

# check whether Boost::configure is used or not
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
  
#
# configure boost if needed
#

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

dnl
dnl Doxygen
dnl

AC_PATH_PROG([ac_cv_prog_doxygen], [doxygen])
AC_SUBST(DOXYGEN)
test -n "$ac_cv_prog_doxygen" && DOXYGEN="$ac_cv_prog_doxygen"

dnl
dnl xsltproc
dnl

AC_XSLTPROC
AC_SUBST(XSLTPROC)
test -n "$ac_cv_have_xsltproc" && XSLTPROC="$ac_cv_prog_xsltproc"

dnl
dnl Java
dnl

AC_PATH_PROG([ac_cv_prog_java], [java])
AC_SUBST(JAVA)
test -n "$ac_cv_prog_java" && JAVA="$ac_cv_prog_java"

dnl
dnl for Generating Documentation (DocBook DTD, XSL, Apache FOP)
dnl

AC_ARG_ENABLE(boostbook,
  AC_HELP_STRING([--enable-boostbook],
    [setup boostbook documentation system]),
  [
  if test "x$enableval" != "xno"; then
    ac_cv_boostbook=yes
  fi
  ]
)
test -z "$ac_cv_boostbook" && ac_cv_boostbook=no
if test -z "$BOOST_SRCDIR" || test "$ac_cv_boost_version" = 1_32; then
  ac_cv_boostbook=no
fi

if test "$ac_cv_boostbook" = yes; then
  AC_DOCBOOK_DTD
  AC_SUBST(DOCBOOK_DTD_DIR)
  test "$ac_cv_have_docbook_dtd" = yes && DOCBOOK_DTD_DIR="$ac_cv_docbook_dtd_dir"

  AC_DOCBOOK_XSL
  AC_SUBST(DOCBOOK_XSL_DIR)
  test "$ac_cv_have_docbook_xsl" = yes && DOCBOOK_XSL_DIR="$ac_cv_docbook_xsl_dir"

  AC_FOP
  AC_SUBST(FOP_DIR)
  AC_SUBST(FOP_BIN)
  AC_SUBST(FOP_OPTS)
  AC_SUBST(JAVA_OPTS)
  if test "$ac_cv_have_fop" = yes; then
    FOP_DIR="$ac_cv_fop_dir"
    FOP_BIN="$ac_cv_fop_bin"
    case "x$build_os" in
      xdarwin*)
        FOP_OPTS="-Xmx128000000"
        ;;
      x*linux*)
        if test -f "/etc/debian_version"; then
          JAVA_OPTS="-Xmx128000000"
	fi
	;;
      *)
        ;;
    esac
  fi
fi

dnl check if all the required tools are present
if test "$ac_cv_boostbook" = yes && test -n "$ac_cv_prog_doxygen" && test -n "$ac_cv_prog_xsltproc" && test "$ac_cv_have_docbook_dtd" = yes && test "$ac_cv_have_docbook_xsl" = yes; then
  build_doc=yes
  if test "$ac_cv_have_fop" = yes; then
    build_pdf=yes
  else
    build_pdf=no
    AC_MSG_NOTICE([PDF documentation will not be built.])
  fi
else
  build_doc=no
  build_pdf=no
  AC_MSG_NOTICE([Documentation will not be built.])
fi

dnl
dnl # setup BASE_CPPFLAGS, BASE_LDFLAGS, BASE_LIBS
dnl

AC_SUBST(BASE_CPPFLAGS)
BASE_CPPFLAGS="$CPPFLAGS"
ac_cv_base_cppflags="$CPPFLAGS"
CPPFLAGS=
AC_SUBST(BASE_LDFLAGS)
BASE_LDFLAGS="$LDFLAGS"
ac_cv_base_ldflags="$LDFLAGS"
LDFLAGS=
AC_SUBST(BASE_LIBS)
BASE_LIBS="$LIBS"
ac_cv_base_libs="$LIBS"
LIBS=

ac_cv_alps_cppflags="$ac_cv_base_cppflags $ac_cv_hdf5_cppflags $ac_cv_lapack_cppflags $ac_cv_sqlite_cppflags $ac_cv_xml_cppflags -I$ac_cv_alps_includedir $ac_cv_boost_cppflags $ac_cv_mpi_cppflags"
ac_cv_alps_ldflags="$ac_cv_base_ldflags $ac_cv_hdf5_ldflags $ac_cv_lapack_ldflags $ac_cv_sqlite_ldflags $ac_cv_xml_ldflags -L$ac_cv_alps_libdir $ac_cv_boost_ldflags"
ac_cv_alps_libs="-lalps -lcomm-sgl $ac_cv_boost_libs $ac_cv_hdf5_libs $ac_cv_lapack_libs $ac_cv_sqlite_libs $ac_cv_xml_libs $ac_cv_base_libs"

ac_cv_alps_ldflags_mpi="$ac_cv_alps_ldflags $ac_cv_mpi_ldflags"
ac_cv_alps_libs_mpi="-lalps -lcomm-mpi $ac_cv_mpi_libs $ac_cv_boost_libs $ac_cv_hdf5_libs $ac_cv_lapack_libs $ac_cv_sqlite_libs $ac_cv_xml_libs $ac_cv_base_libs"

dnl
dnl Experimental Features
dnl

AC_ARG_ENABLE(libheap,
  AC_HELP_STRING([--enable-libheap], [building libheap wrapper @<:@default=yes@:>@]),
  [
  if test "x$enableval" = "xno"; then
    ac_cv_libheap=no
  fi
  ]
)
test -z "$ac_cv_libheap" && ac_cv_libheap=yes

dnl
dnl Absolute paths
dnl

AC_SUBST(abs_top_srcdir)
abs_top_srcdir=`(cd $top_srcdir && pwd)`

dnl
dnl Conditionals
dnl

AM_CONDITIONAL(ALPS_FULL, test "$ac_cv_alps_type" = "full")
AM_CONDITIONAL(ALPS_LIGHT, test "$ac_cv_alps_type" = "light")
AM_CONDITIONAL(BOOST_1_32, test "$ac_cv_boost_version" = 1_32)
AM_CONDITIONAL(BOOST_1_33, test "$ac_cv_boost_version" = 1_33)
AM_CONDITIONAL(BOOST_1_34, test "$ac_cv_boost_version" = 1_34)
AM_CONDITIONAL(BOOST_CVS, test "$ac_cv_boost_version" = cvs)
AM_CONDITIONAL(BUILD_BOOST, test "$ac_cv_use_precompiled_boost" != yes)
AM_CONDITIONAL(BUILD_BOOST_SIGNALS, test "$ac_cv_boost_signals" = yes)
AM_CONDITIONAL(BUILD_BOOST_THREAD, test "$ac_cv_boost_thread" = yes)
AM_CONDITIONAL(BUILD_BOOST_WCHAR, test "$ac_cv_boost_wchar" = yes)
AM_CONDITIONAL(BUILD_DOC, test "$build_doc" = yes)
AM_CONDITIONAL(BUILD_PDF, test "$build_pdf" = yes)
AM_CONDITIONAL(BUILD_LIBHEAP, test "$ac_cv_libheap" = yes)
AM_CONDITIONAL(BUILD_ARCHIVE, test "$ac_cv_have_sqlite" = yes)
AM_CONDITIONAL(HAVE_MPI, test "$ac_cv_have_mpi" = yes)

dnl
dnl Generate files
dnl

dnl alps/config.h
AC_CONFIG_HEADERS([src/alps/config.h])

dnl Makefiles
AC_CONFIG_FILES(Makefile config/config.site config/include.mk)
AC_CONFIG_FILES(src/Makefile)
AC_CONFIG_FILES(lib/Makefile lib/alps/Makefile lib/comm-sgl/Makefile lib/comm-mpi/Makefile lib/heap/Makefile)

AC_CONFIG_FILES(test/Makefile test/fixed_capacity/Makefile test/osiris/Makefile test/osiris/xdrdump2.input test/osiris/boostdump2.input test/parameter/Makefile test/parser/Makefile test/alea/Makefile)
AC_CONFIG_FILES(doc/Makefile)

if test "$ac_cv_alps_type" = "full"; then
  AC_CONFIG_FILES(test/lattice/Makefile test/lattice/parameters  test/lattice/parameters_complex)
  AC_CONFIG_FILES(test/model/Makefile)
  AC_CONFIG_FILES(example/Makefile example/model/Makefile example/scheduler/Makefile example/model/parm_symbolic example/model/parm_numeric)
  AC_CONFIG_FILES(lib/xml/Makefile lib/xml/changestylesheet.xsl lib/xml/lattices.xml lib/xml/models.xml lib/xml/plot2xsl.xsl)
  AC_CONFIG_FILES(tool/Makefile tool/convert2text tool/convert2html tool/extractgp tool/extracthtml tool/extracttext tool/extractxmgr tool/plot2gp tool/plot2html tool/plot2text tool/plot2xmgr tool/changestylesheet)
fi

if test "$ac_cv_use_precompiled_boost" != yes; then
  AC_DEFINE_UNQUOTED(PO_UTF8_CODECVT_FACET_CPP, "$BOOST_SRCDIR/program_options/src/utf8_codecvt_facet.cpp")
  AC_CONFIG_HEADERS(src/boost/po_utf8_codecvt_facet.cpp)
  AC_DEFINE_UNQUOTED(SR_UTF8_CODECVT_FACET_CPP, "$BOOST_SRCDIR/serialization/src/utf8_codecvt_facet.cpp")
  AC_CONFIG_HEADERS(src/boost/sr_utf8_codecvt_facet.cpp)
  if test "$ac_cv_boost_version" = 1_34 || test "$ac_cv_boost_version" = cvs; then
    AC_DEFINE_UNQUOTED(FS_UTF8_CODECVT_FACET_CPP, "$BOOST_SRCDIR/filesystem/src/utf8_codecvt_facet.cpp")
    AC_CONFIG_HEADERS(src/boost/fs_utf8_codecvt_facet.cpp)
  fi
fi

dnl doc generation
AC_CONFIG_FILES(doc/doxygen.conf)
if test "$build_doc" = yes; then
  AC_CONFIG_FILES(doc/autodoc)
  AC_CONFIG_FILES(doc/catalog.xml)
  AC_CONFIG_FILES(doc/bb2db)
  AC_CONFIG_FILES(doc/db2html)
  AC_CONFIG_FILES(doc/db2pdf)
  AC_CONFIG_FILES(doc/singledoc)
fi

AC_OUTPUT

dnl
dnl output summary
dnl

echo
echo "Summary of configuration:"
if test "$ac_cv_alps_type" = "full"; then
  echo "  ALPS version        $PACKAGE_VERSION"
else
  echo "  ALPS version        $PACKAGE_VERSION (Light)"
fi
echo "  prefix              $ac_cv_alps_prefix"
echo "  compiler type       $ac_cv_alps_compiler @<:@$ac_cv_alps_cxx@:>@"
echo "  optimization        $ac_cv_alps_compiler_optimization"
echo "  exceptions          $ac_cv_alps_compiler_exceptions"
echo "  warnings            $ac_cv_alps_compiler_warnings"
echo "  MPI support         $ac_cv_alps_have_mpi"
echo "  HDF5 support        $ac_cv_have_hdf5"
echo "  LAPACK support      $ac_cv_have_lapack"
echo "  SQLite support      $ac_cv_have_sqlite"
echo "  XML parser          $ac_cv_xml_parser"
echo "  Documentation       $build_doc"

if test "$ac_cv_alps_type" = "light"; then
  echo
  echo "*** This is the \"ALPS Light\" Libraries, public-domain part of the ALPS"
  echo "*** Libraries. If you need the full functionality of the ALPS Libraries,"
  echo "*** such as Lattice, Model, Scheduler, etc, please use the full version of"
  echo "*** ALPS Libraries, which is available from http://alps.comp-phys.org/."
fi
