dnl configuration script for ALPS
dnl Process this file with autoconf to produce a configure script.

dnl disable loading config.site
define([AC_SITE_LOAD], )dnl
AC_PREREQ(2.59)

dnl
dnl Version Section ---------------------------------------------------------
dnl

dnl Version of ALPS Library

define([alps_version], [1.3.3])

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

if test -n "$HOST"; then
  AC_DEFINE_UNQUOTED([ALPS_CONFIG_HOST], ["$HOST"])
else
  test -z "$ac_hostname" && ac_hostname="unknown"
  AC_DEFINE_UNQUOTED([ALPS_CONFIG_HOST], ["$ac_hostname"])
fi
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
# checking boost directory and include files
#

AC_BOOST

#
# checking for other C++ headers
#

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

# we could use this to check for the HDF5 C++ interface
#AC_HDF5_CPP
#test "$ac_cv_have_hdf5_cpp" = yes && AC_DEFINE(ALPS_HAVE_HDF5_CPP, [], [hdf5_cpp])

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

dnl
dnl boost object libraries
dnl

AC_BOOST_LIBS

dnl check whether Boost.Configure is used or not
AC_BOOST_CONFIGURE

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
ac_cv_alps_libs="-lparapack-sgl -lalps -lcomm-sgl $ac_cv_boost_libs $ac_cv_hdf5_libs $ac_cv_lapack_libs $ac_cv_sqlite_libs $ac_cv_xml_libs $ac_cv_base_libs"

ac_cv_alps_ldflags_mpi="$ac_cv_alps_ldflags $ac_cv_mpi_ldflags"
ac_cv_alps_libs_mpi="-lparapack-mpi -lalps -lcomm-mpi $ac_cv_mpi_libs $ac_cv_boost_libs $ac_cv_hdf5_libs $ac_cv_lapack_libs $ac_cv_sqlite_libs $ac_cv_xml_libs $ac_cv_base_libs"

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
  AC_CONFIG_FILES(test/parapack/Makefile)
  AC_CONFIG_FILES(example/Makefile example/model/Makefile example/scheduler/Makefile example/model/parm_symbolic example/model/parm_numeric example/parapack/Makefile example/parapack/single/Makefile example/parapack/multiple/Makefile example/parapack/loop/Makefile example/parapack/exchange/Makefile example/parapack/wanglandau/Makefile example/sampling/Makefile)
  AC_CONFIG_FILES(lib/xml/Makefile lib/xml/changestylesheet.xsl lib/xml/lattices.xml lib/xml/models.xml lib/xml/plot2xsl.xsl)
  AC_CONFIG_FILES(tool/Makefile tool/convert2text tool/convert2html tool/extractgp tool/extracthtml tool/extracttext tool/extractxmgr tool/plot2gp tool/plot2html tool/plot2text tool/plot2xmgr tool/changestylesheet)
fi

if test "$ac_cv_use_precompiled_boost" != yes; then
  AC_DEFINE_UNQUOTED(PO_UTF8_CODECVT_FACET_CPP, "$BOOST_SRCDIR/program_options/src/utf8_codecvt_facet.cpp")
  AC_CONFIG_HEADERS(src/boost/po_utf8_codecvt_facet.cpp)
  AC_DEFINE_UNQUOTED(SR_UTF8_CODECVT_FACET_CPP, "$BOOST_SRCDIR/serialization/src/utf8_codecvt_facet.cpp")
  AC_CONFIG_HEADERS(src/boost/sr_utf8_codecvt_facet.cpp)
  if test "$ac_cv_boost_version" = 1_34 || test "$ac_cv_boost_version" = 1_35 || test "$ac_cv_boost_version" = 1_36 || test "$ac_cv_boost_version" = 1_37 || test "$ac_cv_boost_version" = svn; then
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
echo "  Boost version       $ac_cv_boost_version"
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
