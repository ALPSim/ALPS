dnl Process this file with autoconf to produce a configure script.

dnl disable saving cache
define([AC_CACHE_SAVE], )dnl
AC_PREREQ(2.59)

dnl
dnl Version Section ---------------------------------------------------------
dnl

dnl name of package
define([package_name], [alps-applications])

dnl version of package
define([package_version], [1.4.0d])

dnl 
dnl Initialization ----------------------------------------------------------
dnl

AC_INIT(package_name, package_version)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT([$HOME/ALPS])
AM_MAINTAINER_MODE

AC_DEFINE_UNQUOTED([ALPS_APPLICATIONS_VERSION], ["package_version"],
                   [Version of ALPS Applications])

undefine([package_name])
undefine([package_version])

dnl
dnl check for ALPS
dnl

AC_ALPS

if test "$ac_cv_have_lapack" != yes; then
  AC_MSG_NOTICE([LAPACK library not found. diag/fulldiag will not be built.])
  AC_MSG_NOTICE([LAPACK library not found. diag/sparsediag will not be built.])
  AC_MSG_NOTICE([LAPACK library not found. dmrg/dmrg will not be built.])
  AC_MSG_NOTICE([LAPACK library not found. dmrg/interacting will not be built.])
  AC_MSG_NOTICE([LAPACK library not found. dmrg/particle_in_box will not be built.])
  AC_MSG_NOTICE([LAPACK library not found. dmft/DMFT will not be built.])
  AC_MSG_NOTICE([LAPACK library not found. mc/spins will not be built.])
  AC_MSG_NOTICE([LAPACK library not found. qmc/looper will not be built.])
fi

if test "$ac_cv_boost_mpi" != yes; then
  AC_MSG_NOTICE([Boost.MPI library not found. MPI version of qmc/looper will not be built.])
fi

dnl
dnl compiler
dnl

if test -z "${enable_shared+set}" && test -n "$ac_cv_enable_shared"; then
  enable_shared="$ac_cv_enable_shared"
fi
if test -z "${enable_static+set}" && test -n "$ac_cv_enable_static"; then
  enable_static="$ac_cv_enable_static"
fi

AC_PROG_CC
AC_PROG_CXX
F77=no
AC_PROG_LIBTOOL

dnl
dnl check for IETL library
dnl

if test "$ac_cv_have_lapack" = yes; then
  AC_IETL
  if test "$ac_cv_have_ietl" != yes; then
    AC_MSG_NOTICE([IETL library not found. sparsediag will not be built.])
    AC_MSG_NOTICE([IETL library not found. dmrg/interacting will not be built.])
  fi
fi
AM_CONDITIONAL(HAVE_IETL, test "$ac_cv_have_ietl" = yes)

dnl
dnl check for lp_solve
dnl

AC_LPSOLVE
if test "$ac_cv_have_lpsolve" != yes; then
  AC_MSG_NOTICE([lp_solve library not found. qmc/dirloop_sse will not be built.])
fi
AM_CONDITIONAL(HAVE_LPSOLVE, test "$ac_cv_have_lpsolve" = yes)

dnl
dnl compiler flags
dnl

CC="$ac_cv_alps_cc"
CFLAGS="$ac_cv_alps_cflags"
CPP="$ac_cv_alps_cpp"
CPPFLAGS="$ac_cv_alps_cppflags $IETL_CPPFLAGS $LPSOLVE_CPPFLAGS"
CXX="$ac_cv_alps_cxx"
CXXFLAGS="$ac_cv_alps_cxxflags"
CXXCPP="$ac_cv_alps_cxxcpp"
LDFLAGS="$ac_cv_alps_ldflags $LPSOLVE_LDFLAGS"
LIBS=

AC_SUBST([DEFAULT_INCLUDES], '-I. -I$(srcdir) -I$(top_builddir)')
AC_SUBST([LDADD], "$LPSOLVE_LIBS $ac_cv_alps_libs")

AC_SUBST([LDFLAGS_MPI], "$ac_cv_alps_ldflags_mpi $LPSOLVE_LDFLAGS")
AC_SUBST([LIBS_MPI], "$LPSOLVE_LIBS $ac_cv_alps_libs_mpi")

dnl
dnl Applications
dnl

dnl diag
exist_diag=no
exist_diag_fulldiag=no
build_diag_fulldiag=no
exist_diag_sparsediag=no
build_diag_sparsediag=no
if test -f "$srcdir/diag/Makefile.in"; then
  exist_diag=yes
  dnl diag/fulldiag
  if test -f "$srcdir/diag/fulldiag/Makefile.in"; then
    exist_diag_fulldiag=yes
    if test "$ac_cv_have_lapack" = yes; then
      build_diag_fulldiag=yes
    fi
  fi
  dnl diag/sparsediag
  if test -f "$srcdir/diag/sparsediag/Makefile.in"; then
    exist_diag_sparsediag=yes
    if test "$ac_cv_have_lapack" = yes && test "$ac_cv_have_ietl" = yes; then
      build_diag_sparsediag=yes
    fi
  fi
fi
AM_CONDITIONAL(EXIST_DIAG, test "$exist_diag" = yes)
AM_CONDITIONAL(EXIST_DIAG_FULLDIAG, test "$exist_diag_fulldiag" = yes)
AM_CONDITIONAL(BUILD_DIAG_FULLDIAG, test "$build_diag_fulldiag" = yes)
AM_CONDITIONAL(EXIST_DIAG_SPARSEDIAG, test "$exist_diag_sparsediag" = yes)
AM_CONDITIONAL(BUILD_DIAG_SPARSEDIAG, test "$build_diag_sparsediag" = yes)

dnl dmrg
exist_dmrg=no
exist_dmrg_dmrg=no
build_dmrg_dmrg=no
exist_dmrg_interacting=no
build_dmrg_interacting=no
exist_dmrg_particle_in_box=no
build_dmrg_particle_in_box=no
if test -f "$srcdir/dmrg/Makefile.in"; then
  exist_dmrg=yes
  dnl dmrg/interacting
  if test -f "$srcdir/dmrg/dmrg/Makefile.in"; then
    exist_dmrg_dmrg=yes
    if test "$ac_cv_have_lapack" = yes; then
      build_dmrg_dmrg=yes
    fi
  fi
  if test -f "$srcdir/dmrg/interacting/Makefile.in"; then
    exist_dmrg_interacting=yes
    if test "$ac_cv_have_lapack" = yes && test "$ac_cv_have_ietl" = yes; then
      build_dmrg_interacting=yes
    fi
  fi
  dnl dmrg/particle_in_box
  if test -f "$srcdir/dmrg/particle_in_box/Makefile.in"; then
    exist_dmrg_particle_in_box=yes
    if test "$ac_cv_have_lapack" = yes; then
      build_dmrg_particle_in_box=yes
    fi
  fi
fi
AM_CONDITIONAL(EXIST_DMRG, test "$exist_dmrg" = yes)
AM_CONDITIONAL(EXIST_DMRG_DMRG, test "$exist_dmrg_dmrg" = yes)
AM_CONDITIONAL(BUILD_DMRG_DMRG, test "$build_dmrg_dmrg" = yes)
AM_CONDITIONAL(EXIST_DMRG_INTERACTING, test "$exist_dmrg_interacting" = yes)
AM_CONDITIONAL(BUILD_DMRG_INTERACTING, test "$build_dmrg_interacting" = yes)
AM_CONDITIONAL(EXIST_DMRG_PARTICLE_IN_BOX, test "$exist_dmrg_particle_in_box" = yes)
AM_CONDITIONAL(BUILD_DMRG_PARTICLE_IN_BOX, test "$build_dmrg_particle_in_box" = yes)

dnl dmft
exist_dmft=no
exist_dmft_dmft=no
build_dmft_dmft=no
exist_dmft_solvers=no
build_dmft_solvers=no
exist_dmft_solvers_hybridization=no
build_dmft_solvers_hybridization=no

#the following line will switch on the DMFT code
#if test -f "$srcdir/dmft/Makefile.in"; then
if false; then
  exist_dmft=yes
  dnl dmft/DMFT
  if test -f "$srcdir/dmft/DMFT/Makefile.in"; then
    exist_dmft_dmft=yes
    if test "$ac_cv_have_lapack" = yes; then
      build_dmft_dmft=yes
    fi
  fi
  if test -f "$srcdir/dmft/DMFT-solvers/Makefile.in"; then
    exist_dmft_solvers=yes
    if test "$ac_cv_have_lapack" = yes; then
      build_dmft_solvers=yes
    fi
  fi
  dnl dmft/DMFT-solvers/Werner
  if test -f "$srcdir/dmft/DMFT-solvers/Werner/Makefile.in"; then
    exist_dmft_solvers_hybridization=yes
    if test "$ac_cv_have_lapack" = yes; then
      build_dmft_solvers_hybridization=yes
    fi
  fi
fi
AM_CONDITIONAL(EXIST_DMFT, test "$exist_dmft" = yes)
AM_CONDITIONAL(EXIST_DMFT_DMFT, test "$exist_dmft_dmft" = yes)
AM_CONDITIONAL(BUILD_DMFT_DMFT, test "$build_dmft_dmft" = yes)
AM_CONDITIONAL(EXIST_DMFT_SOLVERS, test "$exist_dmft_solver" = yes)
AM_CONDITIONAL(BUILD_DMFT_SOLVERS, test "$build_dmft_solver" = yes)
AM_CONDITIONAL(EXIST_DMFT_SOLVERS_HYBRIDIZATION, test "$exist_dmft_solvers_hybridization" = yes)
AM_CONDITIONAL(BUILD_DMFT_SOLVERS_HYBRIDIZATION, test "$build_dmft_solvers_hybridization" = yes)

dnl mc
exist_mc=no
exist_mc_spins=no
build_mc_spins=no
if test -f "$srcdir/mc/Makefile.in"; then
  exist_mc=yes
  dnl mc/spins
  if test -f "$srcdir/mc/spins/Makefile.in"; then
    exist_mc_spins=yes
    if test "$ac_cv_have_lapack" = yes; then
      build_mc_spins=yes
    fi
  fi
fi
AM_CONDITIONAL(EXIST_MC, test "$exist_mc" = yes)
AM_CONDITIONAL(EXIST_MC_SPINS, test "$exist_mc_spins" = yes)
AM_CONDITIONAL(BUILD_MC_SPINS, test "$build_mc_spins" = yes)

dnl qmc
exist_qmc=no
exist_qmc_checksign=no
build_qmc_checksign=no
exist_qmc_looper=no
build_qmc_looper=no
exist_qmc_sse=no
exist_qmc_sse2=no
build_qmc_sse=no
exist_qmc_sse_loop=no
build_qmc_sse_loop=no
exist_qmc_worms=no
build_qmc_worms=no
exist_qmc_qwl=no
build_qmc_qwl=no
if test -f "$srcdir/qmc/Makefile.in"; then
  exist_qmc=yes
  dnl qmc/checksign
  if test -f "$srcdir/qmc/checksign/Makefile.in"; then
    exist_qmc_checksign=yes
    build_qmc_checksign=yes
  fi
  dnl qmc/looper
  if test -f "$srcdir/qmc/looper/Makefile.in"; then
    exist_qmc_looper=yes
    if test "$ac_cv_have_lapack" = yes; then
      build_qmc_looper=yes
    fi
  fi
  dnl qmc/sse
  if test -f "$srcdir/qmc/sse/Makefile.in"; then
    exist_qmc_sse=yes
    if test "$ac_cv_have_lpsolve" = yes; then
      build_qmc_sse=yes
    fi
  fi
#  dnl qmc/sse2
#  if test -f "$srcdir/qmc/sse2/Makefile.in"; then
#    exist_qmc_sse2=yes
#    if test "$ac_cv_have_lpsolve" = yes; then
#      build_qmc_sse2=yes
#    fi
#  fi
  dnl qmc/sse-loop
  if test -f "$srcdir/qmc/sse-loop/Makefile.in"; then
    exist_qmc_sse_loop=yes
    build_qmc_sse_loop=yes
  fi
  dnl qmc/worms
  if test -f "$srcdir/qmc/worms/Makefile.in"; then
    exist_qmc_worms=yes
    build_qmc_worms=yes
  fi
  dnl qmc/qwl
  if test -f "$srcdir/qmc/qwl/Makefile.in"; then
    exist_qmc_qwl=yes
    build_qmc_qwl=yes
  fi
fi
AM_CONDITIONAL(EXIST_QMC, test "$exist_qmc" = yes)
AM_CONDITIONAL(EXIST_QMC_CHECKSIGN, test "$exist_qmc_checksign" = yes)
AM_CONDITIONAL(BUILD_QMC_CHECKSIGN, test "$build_qmc_checksign" = yes)
AM_CONDITIONAL(EXIST_QMC_LOOPER, test "$exist_qmc_looper" = yes)
AM_CONDITIONAL(BUILD_QMC_LOOPER, test "$build_qmc_looper" = yes)
AM_CONDITIONAL(EXIST_QMC_SSE, test "$exist_qmc_sse" = yes)
AM_CONDITIONAL(BUILD_QMC_SSE, test "$build_qmc_sse" = yes)
AM_CONDITIONAL(EXIST_QMC_SSE2, test "$exist_qmc_sse2" = yes)
AM_CONDITIONAL(BUILD_QMC_SSE2, test "$build_qmc_sse2" = yes)
AM_CONDITIONAL(EXIST_QMC_SSE_LOOP, test "$exist_qmc_sse_loop" = yes)
AM_CONDITIONAL(BUILD_QMC_SSE_LOOP, test "$build_qmc_sse_loop" = yes)
AM_CONDITIONAL(EXIST_QMC_WORMS, test "$exist_qmc_worms" = yes)
AM_CONDITIONAL(BUILD_QMC_WORMS, test "$build_qmc_worms" = yes)
AM_CONDITIONAL(EXIST_QMC_QWL, test "$exist_qmc_qwl" = yes)
AM_CONDITIONAL(BUILD_QMC_QWL, test "$build_qmc_qwl" = yes)

dnl
dnl Tests
dnl

dnl test
if test -f "$srcdir/test/Makefile.in"; then
  exist_test=yes
  if test -f "$srcdir/test/fulldiag/Makefile.in"; then
    exist_test_fulldiag=yes
  fi
  if test -f "$srcdir/test/sparsediag/Makefile.in"; then
    exist_test_sparsediag=yes
  fi
  if test -f "$srcdir/test/checksign/Makefile.in"; then
    exist_test_checksign=yes
  fi
  if test -f "$srcdir/test/looper/Makefile.in"; then
    exist_test_looper=yes
  fi
  if test -f "$srcdir/test/sse/Makefile.in"; then
    exist_test_sse=yes
  fi
  if test -f "$srcdir/test/sse-loop/Makefile.in"; then
    exist_test_sse_loop=yes
  fi
  if test -f "$srcdir/test/worms/Makefile.in"; then
    exist_test_worms=yes
  fi
  if test -f "$srcdir/test/qwl/Makefile.in"; then
    exist_test_qwl=yes
  fi
  if test -f "$srcdir/test/mc/Makefile.in"; then
    exist_test_mc=yes
    if test -f "$srcdir/test/mc/spins/Makefile.in"; then
      exist_test_mc_spins=yes
    fi
  fi
fi
AM_CONDITIONAL(EXIST_TEST, test "$exist_test" = yes)
AM_CONDITIONAL(EXIST_TEST_FULLDIAG, test "$exist_test_fulldiag" = yes)
AM_CONDITIONAL(EXIST_TEST_SPARSEDIAG, test "$exist_test_sparsediag" = yes)
AM_CONDITIONAL(EXIST_TEST_CHECKSIGN, test "$exist_test_checksign" = yes)
AM_CONDITIONAL(EXIST_TEST_LOOPER, test "$exist_test_looper" = yes)
AM_CONDITIONAL(EXIST_TEST_SSE, test "$exist_test_sse" = yes)
AM_CONDITIONAL(EXIST_TEST_SSE_LOOP, test "$exist_test_sse_loop" = yes)
AM_CONDITIONAL(EXIST_TEST_WORMS, test "$exist_test_worms" = yes)
AM_CONDITIONAL(EXIST_TEST_QWL, test "$exist_test_qwl" = yes)
AM_CONDITIONAL(EXIST_TEST_MC, test "$exist_test_mc" = yes)
AM_CONDITIONAL(EXIST_TEST_MC_SPINS, test "$exist_test_mc_spins" = yes)

dnl
dnl generate files
dnl

AC_CONFIG_HEADER(config/config.h)
AC_CONFIG_FILES(Makefile)
test "$exist_diag" = yes && AC_CONFIG_FILES(diag/Makefile)
test "$exist_diag_fulldiag" = yes && AC_CONFIG_FILES(diag/fulldiag/Makefile)
test "$exist_diag_sparsediag" = yes && AC_CONFIG_FILES(diag/sparsediag/Makefile)
test "$exist_dmrg" = yes && AC_CONFIG_FILES(dmrg/Makefile)
test "$exist_dmrg_dmrg" = yes && AC_CONFIG_FILES(dmrg/dmrg/Makefile)
test "$exist_dmrg_interacting" = yes && AC_CONFIG_FILES(dmrg/interacting/Makefile)
test "$exist_dmrg_particle_in_box" = yes && AC_CONFIG_FILES(dmrg/particle_in_box/Makefile)
test "$exist_dmft" = yes && AC_CONFIG_FILES(dmft/Makefile)
test "$exist_dmft_dmft" = yes && AC_CONFIG_FILES(dmft/DMFT/Makefile)
test "$exist_dmft_solver" = yes && AC_CONFIG_FILES(dmft/DMFT-Solvers/Makefile)
test "$exist_dmft_solver_hybridization" = yes && AC_CONFIG_FILES(dmft/DMFT-Solvers/Werner/Makefile)
test "$exist_mc" = yes && AC_CONFIG_FILES(mc/Makefile)
test "$exist_mc_spins" = yes && AC_CONFIG_FILES(mc/spins/Makefile)
test "$exist_qmc" = yes && AC_CONFIG_FILES(qmc/Makefile)
test "$exist_qmc_checksign" = yes && AC_CONFIG_FILES(qmc/checksign/Makefile)
test "$exist_qmc_looper" = yes && AC_CONFIG_FILES(qmc/looper/Makefile)
test "$exist_qmc_sse" = yes && AC_CONFIG_FILES(qmc/sse/Makefile)
#test "$exist_qmc_sse2" = yes && AC_CONFIG_FILES(qmc/sse2/Makefile)
test "$exist_qmc_sse_loop" = yes && AC_CONFIG_FILES(qmc/sse-loop/Makefile)
test "$exist_qmc_worms" = yes && AC_CONFIG_FILES(qmc/worms/Makefile)
test "$exist_qmc_qwl" = yes && AC_CONFIG_FILES(qmc/qwl/Makefile)
test "$exist_test" = yes && AC_CONFIG_FILES(test/Makefile test/parm test/parm_fermion test/parm_fermion2 test/parm_boson test/energy.xsl)
test "$exist_test_fulldiag" = yes && AC_CONFIG_FILES(test/fulldiag/Makefile)
test "$exist_test_sparsediag" = yes && AC_CONFIG_FILES(test/sparsediag/Makefile test/sparsediag/groundstate.xsl)
test "$exist_test_checksign" = yes && AC_CONFIG_FILES(test/checksign/Makefile)
test "$exist_test_looper" = yes && AC_CONFIG_FILES(test/looper/Makefile)
test "$exist_test_sse" = yes && AC_CONFIG_FILES(test/sse/Makefile)
test "$exist_test_sse_loop" = yes && AC_CONFIG_FILES(test/sse-loop/Makefile)
test "$exist_test_worms" = yes && AC_CONFIG_FILES(test/worms/Makefile)
test "$exist_test_qwl" = yes && AC_CONFIG_FILES(test/qwl/Makefile test/qwl/parm)
if test "$exist_test_mc" = yes; then
  AC_CONFIG_FILES([test/mc/Makefile])
  if test "$exist_test_mc_spins" = yes; then
    AC_CONFIG_FILES([test/mc/spins/Makefile])
    AC_CONFIG_FILES([test/mc/spins/check.sh])
    AC_CONFIG_FILES([test/mc/spins/ising_local.parm])
    AC_CONFIG_FILES([test/mc/spins/ising_cluster.parm])
    AC_CONFIG_FILES([test/mc/spins/heisenberg_chain.parm])
    AC_CONFIG_FILES([test/mc/spins/heisenberg_ladder.parm])
  fi
fi
AC_OUTPUT

dnl
dnl output summary
dnl

echo
echo "Summary of configuration:"
if test "x$exist_diag" = xyes; then
  echo "  Diagonalization"
  if test "x$exist_diag_fulldiag" = xyes; then
    echo "    Full                 $build_diag_fulldiag"
  fi
  if test "x$exist_diag_sparsediag" = xyes; then
    echo "    Sparse               $build_diag_sparsediag"
  fi
fi
if test "x$exist_dmrg" = xyes; then
  echo "  DMRG"
  if test "x$exist_dmrg_dmrg" = xyes; then
    echo "    DMRG                 $build_dmrg_dmrg"
  fi
  if test "x$exist_dmrg_interacting" = xyes; then
    echo "    Interacting          $build_dmrg_interacting"
  fi
  if test "x$exist_dmrg_particle_in_box" = xyes; then
    echo "    Particle in box      $build_dmrg_particle_in_box"
  fi
fi
if test "x$exist_dmft" = xyes; then
  echo "  DMFT"
  if test "x$exist_dmft_dmft" = xyes; then
    echo "    DMFT                 $build_dmft_dmft"
  fi
  if test "x$exist_dmft_solvers" = xyes; then
    echo "    Solvers              $build_dmft_solvers"
  fi
  if test "x$exist_dmft_solvers" = xyes; then
    echo "     Hirsch Fye          $build_dmft_solvers"
  fi
  if test "x$exist_dmft_solvers_hybridization" = xyes; then
    echo "     Hybridization       $build_dmft_solvers_hybridization"
  fi
fi
if test "x$exist_mc" = xyes; then
  echo "  Classical Monte Carlo"
  if test "x$exist_mc_spins" = xyes; then
    echo "    Spin                 $build_mc_spins"
  fi
fi
if test "x$exist_qmc" = xyes; then
  echo "  Quantum Monte Carlo"
  if test "x$exist_qmc_checksign" = xyes; then
    echo "    Checksign            $build_qmc_checksign"
  fi
  if test "x$exist_qmc_looper" = xyes; then
    echo "    Looper               $build_qmc_looper"
  fi
  if test "x$exist_qmc_sse" = xyes; then
    echo "    SSE-dirloop          $build_qmc_sse"
  fi
#  if test "x$exist_qmc_sse2" = xyes; then
#    echo "    SSE-dirloop-2        $build_qmc_sse2"
#  fi
  if test "x$exist_qmc_sse_loop" = xyes; then
    echo "    SSE-loop             $build_qmc_sse_loop"
  fi
  if test "x$exist_qmc_worms" = xyes; then
    echo "    Worms                $build_qmc_worms"
  fi
  if test "x$exist_qmc_qwl" = xyes; then
    echo "    Quantum Wang-Landau  $build_qmc_qwl"
  fi
fi
