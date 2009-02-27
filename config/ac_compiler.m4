AC_DEFUN([AC_COMPILER],
  [
  AC_SUBST(COMPILER)
  AC_SUBST(CPPFLAGS)
  AC_SUBST(CFLAGS)
  AC_SUBST(CXXFLAGS)

  AC_MSG_CHECKING([for optimization])
  AC_ARG_ENABLE(optimization,
    AC_HELP_STRING([--enable-optimization],
      [enable optimization @<:@default=yes@:>@]),
    [
    if test "x$enableval" = "xno"; then
      ac_cv_compiler_optimization=no
    fi
    ]
  )
  test -z "$ac_cv_compiler_optimization" && ac_cv_compiler_optimization=yes
  AC_MSG_RESULT($ac_cv_compiler_optimization)

  AC_MSG_CHECKING([for exception handling])
  AC_ARG_ENABLE(exceptions,
    AC_HELP_STRING([--enable-exceptions],
      [enable exception handling @<:@default=yes@:>@]),
    [
    if test "x$enableval" = "xno"; then
      ac_cv_compiler_exceptions=no
    fi
    ]
  )
  test -z "$ac_cv_compiler_exceptions" && ac_cv_compiler_exceptions=yes
  AC_MSG_RESULT($ac_cv_compiler_exceptions)

  AC_MSG_CHECKING([for warning messages])
  AC_ARG_ENABLE(warnings,
    AC_HELP_STRING([--enable-warnings],
      [enable warning messages]),
    [
    if test "x$enableval" = "xyes"; then
      ac_cv_compiler_warnings=yes
    fi
    ]
  )
  test -z "$ac_cv_compiler_warnings" && ac_cv_compiler_warnings=no
  AC_MSG_RESULT($ac_cv_compiler_warnings)

  AC_ARG_WITH(compiler,
    AC_HELP_STRING([--with-compiler=MODE],
      [set compiler mode (MODE = gnu, kai, intel, como, hp32, hp64, dec, sgi32, sgi64, cray, cray-gcc, ibm32, ibm64, macos, macos-gcc-3, macos-gcc-3.3, macos-gcc-4.0, cygwin, pgi, fcc, generic)]),
    [
    case "x$withval" in
      xgnu* | xGNU* | xgcc* | xg++* )
        COMPILER="gnu"
        ;;
      xcygwin* )
        COMPILER="cygwin"
        ;;
      xkai* | xKAI*)
        COMPILER="kai"
        ;;
      xintel* | xicc* )
        COMPILER="intel"
        ;;
      xcomo*)
        COMPILER="como"
        ;;
      xhp32* | xaCC32* | xhp | xaCC)
        COMPILER="hp32"
        ;;
      xhp64* | xaCC64*)
        COMPILER="hp64"
        ;;
      xdec* | xcxx*)
        COMPILER="dec"
        ;;
      xsgi32* | xsgi | xirix32* | xirix*)
        COMPILER="sgi32"
        ;;
      xsgi64* | xirix64*)
        COMPILER="sgi64"
        ;;
      xcray-gcc* | xCRAY-gcc*)
        COMPILER="cray-gcc"
        ;;
      xcray* | xCRAY*)
        COMPILER="cray"
        ;;
      xibm32* | xIBM32*)
        COMPILER="ibm32"
        ;;
      xibm64* | xIBM64* | xibm* | xIBM* | xvacpp* | xxlC*)
        COMPILER="ibm64"
        ;;
      xmacos-gcc-3 | xmacos-gnu-3 | xmacos-g++-3)
        COMPILER="macos-gcc-3"
        ;;
      xmacos-gcc-3.3 | xmacos-gnu-3.3 | xmacos-g++-3.3)
        COMPILER="macos-gcc-3.3"
        ;;
      xmacos-gcc-4.0 | xmacos-gnu-4.0 | xmacos-g++-4.0)
        COMPILER="macos-gcc-4.0"
        ;;
      xmacos* | xmac* | xosx*)
        COMPILER="macos"
        ;;
      xpgi)
      	COMPILER="pgi"
	;;
      xfcc)
      	COMPILER="fcc"
	;;
      xgeneric*)
        COMPILER="generic"
        ;;
      *)
        AC_MSG_ERROR([unknown mode $withval])
        ;;
    esac
    ]
  )

  AC_MSG_CHECKING([for compiler mode])
  test -z "$COMPILER" && COMPILER="generic"
  AC_MSG_RESULT([$COMPILER])

  case "$COMPILER" in
    gnu)
      try_CC="gcc"
      try_CXX="g++"
      ;;
    cygwin)
      try_CC="gcc"
      try_CXX="g++"
      ;;
    kai)
      try_CC="cc"
      try_CXX="KCC"
      ;;
    intel*)
      try_CC="icc"
      try_CXX="icpc"
      ;;
    como)
      try_CC="como"
      try_CXX="como"
      ;;
    hp*)
      try_CC="cc"
      try_CXX="aCC"
      ;;
    dec)
      try_CC="cc"
      try_CXX="cxx"
      ;;
    sgi*)
      try_CC="cc"
      try_CXX="CC"
      ;;
    cray)
      try_CC="cc"
      try_CXX="CC"
      ;;
    cray-gcc)
      try_CC="cc"
      try_CXX="CC"
      ;;
    ibm*)
      try_CC="xlc_r"
      try_CXX="xlC_r"
      ;;
    macos)
      try_CC="gcc"
      try_CXX="g++"
      ;;
    macos-gcc-3)
      try_CC="gcc3"
      try_CXX="g++3"
      ;;
    macos-gcc-3.3)
      try_CC="gcc-3.3"
      try_CXX="g++-3.3"
      ;;
    macos-gcc-4.0)
      try_CC="gcc-4.0"
      try_CXX="g++-4.0"
      ;;
    pgi)
      try_CC="pgcc"
      try_CXX="pgCC"
      ;;
    fcc)
      try_CC="fcc"
      try_CXX="FCC"
      ;;
    generic)
      # nothing to do
      ;;
    *)
      AC_MSG_ERROR([unknown mode $COMPILER])
      ;;
  esac

  if test "$COMPILER" != "generic"; then
    if test -z "$CC"; then
      CC=$try_CC
    fi
    if test -z "$CXX"; then
      CXX=$try_CXX
    fi
  fi

  save_CFLAGS="$CFLAGS"
  save_CPPFLAGS="$CPPFLAGS"
  save_CXXFLAGS="$CXXFLAGS"

  AC_PROG_CC
  AC_PROG_CXX

  if test "$COMPILER" = generic && test "$ac_cv_cxx_compiler_gnu" = yes; then
    COMPILER=gnu
    AC_MSG_NOTICE([compiler mode is reset to $COMPILER])
  fi

  test -z "$COMPILER" && AC_MSG_ERROR([compiler mode is not set])

  # default options
  try_CFLAGS=
  try_CFLAGS_OPT="-O3"
  try_CFLAGS_DEBUG="-g -O0"
  try_CFLAGS_WARN=
  try_CFLAGS_NOWARN=
  try_CXXFLAGS=
  try_CXXFLAGS_OPT="-O3 -DBOOST_DISABLE_ASSERTS"
  try_CXXFLAGS_DEBUG="-g -O0"
  try_CXXFLAGS_WARN=
  try_CXXFLAGS_NOWARN=
  try_CXXFLAGS_EH=
  try_CXXFLAGS_NOEH=
  case "$COMPILER" in
    gnu)
      try_CFLAGS="-pthread"
      try_CFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare -Wno-deprecated"
      try_CFLAGS_NOWARN=""
      try_CXXFLAGS="-pthread -ftemplate-depth-150"
      try_CXXFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare -Wno-deprecated"
      try_CXXFLAGS_NOWARN=""
      try_CXXFLAGS_EH="-fexceptions"
      try_CXXFLAGS_NOEH="-fno-exceptions"
      ;;
    cray-gcc)
      try_CFLAGS=""
      try_CFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare"
      try_CFLAGS_NOWARN="-w"
      try_CXXFLAGS="-fabi-version=0 -DMPICH_SKIP_MPICXX -ftemplate-depth-150"
      try_CXXFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare"
      try_CXXFLAGS_NOWARN="-w"
      try_CXXFLAGS_EH="-fexceptions"
      try_CXXFLAGS_NOEH="-fno-exceptions"
      ;;
    cygwin)
      try_CFLAGS="-mthreads"
      try_CFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare"
      try_CFLAGS_NOWARN="-w"
      try_CXXFLAGS="-mthreads -ftemplate-depth-150 -DBOOST_POSIX -DBOOST_POSIX_API -DBOOST_POSIX_PATH"
      # -DBOOST_POSIX for boost 1.33
      # -DBOOST_POSIX_API -DBOOST_POSIX_PATH for boost 1.34 or later
      try_CXXFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare"
      try_CXXFLAGS_NOWARN="-w"
      try_CXXFLAGS_EH="-fexceptions"
      try_CXXFLAGS_NOEH="-fno-exceptions"
      ;;
    kai)
      try_CXXFLAGS="--restrict --one_instantiation_per_object --thread_safe -DBOOST_REGEX_NO_EXTERNAL_TEMPLATES"
      try_CXXFLAGS_OPT="+K3 -O3"
      try_CXXFLAGS_DEBUG="-g +K0 -O0"
      try_CXXFLAGS_EH="--exceptions"
      try_CXXFLAGS_NOEH="--no_exceptions"
      ;;
    intel)
      try_CFLAGS_WARN=
      try_CFLAGS_NOWARN="-w"
      try_CXXFLAGS_OPT="-O2 -DBOOST_DISABLE_ASSERTS"
      try_CXXFLAGS="-D_REENTRANT -restrict"
      try_CXXFLAGS_WARN=
      try_CXXFLAGS_NOWARN="-w"
      ;;
    como)
      ;;
    hp32)
      try_CFLAGS="-Aa +DA1.1 -D_HPUX_SOURCE"
      try_CXXFLAGS="-Aa +DA1.1 -D_HPUX_SOURCE"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="+noeh"
      ;;
    hp64)
      try_CFLAGS="-Aa +DA2.0W -D_HPUX_SOURCE"
      try_CXXFLAGS="-Aa +DA2.0W -D_HPUX_SOURCE"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="+noeh"
      ;;
    dec)
      try_CFLAGS="-pthread"
      try_CFLAGS_OPT="-O4 "
      try_CFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS="-std ansi -model ansi -noimplicit_include -nousing_std -tweak -pthread -D__USE_STD_IOSTREAM"
      try_CXXFLAGS_OPT="-O4 -DBOOST_DISABLE_ASSERTS"
      try_CXXFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="-noexceptions"
      ;;
    sgi32)
      try_CFLAGS="-n32 -diag_error 1035"
      try_CFLAGS_OPT="-Ofast -INLINE"
      try_CFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS="-n32 -LANG:std -diag_error 1035"
      try_CXXFLAGS_OPT="-Ofast -INLINE -DBOOST_DISABLE_ASSERTS"
      try_CXXFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS_EH="-LANG:exceptions=ON"
      try_CXXFLAGS_NOEH="-LANG:exceptions=OFF"
      CPPFLAGS="$CPPFLAGS -DBOOST_UBLAS_REVERSE_ITERATOR_OVERLOADS"
      if test -n "$TOOLROOT"; then
        CPPFLAGS="$CPPFLAGS -I$TOOLROOT/usr/include/CC"
      fi
      ;;
    sgi64)
      try_CFLAGS="-64 -diag_error 1035"
      try_CFLAGS_OPT="-Ofast -INLINE"
      try_CFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS="-64 -LANG:std -diag_error 1035"
      try_CXXFLAGS_OPT="-Ofast -INLINE -DBOOST_DISABLE_ASSERTS"
      try_CXXFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS_EH="-LANG:exceptions=ON"
      try_CXXFLAGS_NOEH="-LANG:exceptions=OFF"
      CPPFLAGS="$CPPFLAGS -DBOOST_UBLAS_REVERSE_ITERATOR_OVERLOADS"
      if test -n "$TOOLROOT"; then
        CPPFLAGS="$CPPFLAGS -I$TOOLROOT/usr/include/CC"
      fi
      ;;
    cray)
      try_CFLAGS="-h conform"
      try_CFLAGS_OPT="-O2"
      try_CFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS="-h one_instantiation_per_object -h new_for_init -h nodep_name -h parse_templates"
      try_CXXFLAGS_OPT="-O2 -DBOOST_DISABLE_ASSERTS"
      try_CXXFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS_EH="-h exceptions"
      try_CXXFLAGS_NOEH="-h noexceptions"
      ;;
    ibm32)
      try_CFLAGS="-q32"
      try_CFLAGS_NOWARN="-w"
      try_CFLAGS_OPT="-O2"
      try_CFLAGS_DEBUG="-g"
      try_CXXFLAGS="-q32"
      try_CXXFLAGS_NOWARN="-w"
      try_CXXFLAGS_OPT="-O2 -qrtti -DBOOST_DISABLE_ASSERTS"
      try_CXXFLAGS_DEBUG="-g -qrtti"
      ;;
    ibm64)
      try_CFLAGS="-q64"
      try_CFLAGS_NOWARN="-w"
      try_CFLAGS_OPT="-O2"
      try_CFLAGS_DEBUG="-g"
      try_CXXFLAGS="-q64"
      try_CXXFLAGS_NOWARN="-w"
      try_CXXFLAGS_OPT="-O2 -qrtti -DBOOST_DISABLE_ASSERTS"
      try_CXXFLAGS_DEBUG="-g -qrtti"
      ;;
    macos-gcc-3)
      try_CFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare -Wno-long-double"
      try_CFLAGS_NOWARN="-w"
      try_CXXFLAGS="-DUSE_DATE_TIME_PRE_1_33_FACET_IO -DBOOST_DATE_TIME_NO_LOCALE -ftemplate-depth-150"
      try_CXXFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare -Wno-long-double"
      try_CXXFLAGS_NOWARN="-w"
      try_CXXFLAGS_EH="-fexceptions"
      try_CXXFLAGS_NOEH="-fno-exceptions"
      ;;
    macos-gcc-3.3)
      try_CFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare"
      try_CFLAGS_NOWARN="-w"
      try_CXXFLAGS="-ftemplate-depth-150"
      try_CXXFLAGS_OPT="-fabi-version=0 -O3"
      try_CXXFLAGS_DEBUG="-fabi-version=0 -g -O0"
      try_CXXFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare"
      try_CXXFLAGS_NOWARN="-w"
      try_CXXFLAGS_EH="-fexceptions"
      try_CXXFLAGS_NOEH="-fno-exceptions"
      ;;
    macos*)
      try_CFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare -Wno-deprecated"
      try_CFLAGS_NOWARN="-w"
      try_CXXFLAGS="-ftemplate-depth-150"
      try_CXXFLAGS_WARN="-W -Wall -Wno-comment -Wno-sign-compare -Wno-deprecated"
      try_CXXFLAGS_NOWARN="-w"
      try_CXXFLAGS_EH="-fexceptions"
      try_CXXFLAGS_NOEH="-fno-exceptions"
      ;;
    pgi*)
      ;;
    fcc*)
      try_CXXFLAGS_OPT="-Kfast -DBOOST_DISABLE_ASSERTS"
      try_CXXFLAGS_DEBUG="-O0"
      try_CXXFLAGS_NOWARN="-w"
      ;;
    generic)
      try_CFLAGS_OPT="$CFLAGS"
      if test $ac_cv_prog_cc_g = yes; then
        try_CFLAGS_DEBUG="-g"
      else
        try_CFLAGS_DEBUG="$CFLAGS"
      fi
      try_CXXFLAGS_OPT="$CXXFLAGS"
      if test $ac_cv_prog_cxx_g = yes; then
        try_CXXFLAGS_DEBUG="-g"
      else
        try_CXXFLAGS_DEBUG="$CXXFLAGS"
      fi
      ;;
    *)
      AC_MSG_ERROR([unknown mode $COMPILER])
      ;;
  esac

  if test "$ac_cv_compiler_optimization" = yes; then
    CPPFLAGS="$CPPFLAGS -DNDEBUG"
  fi
  if test "$ac_cv_compiler_exceptions" = no; then
    CPPFLAGS="$CPPFLAGS -DBOOST_NO_EXCEPTIONS"
  fi
  CPPFLAGS=`echo $CPPFLAGS | sed 's/^ *//' | sed 's/ *$//' | sed 's/  */ /'`

  if test -n "$save_CFLAGS"; then
    CFLAGS="$save_CFLAGS"
  else
    if test "$ac_cv_compiler_optimization" = yes; then
      CFLAGS="$try_CFLAGS_OPT"
    else
      CFLAGS="$try_CFLAGS_DEBUG"
    fi
    if test "$ac_cv_compiler_warnings" = yes; then
      CFLAGS="$CFLAGS $try_CFLAGS_WARN"
    else
      CFLAGS="$CFLAGS $try_CFLAGS_NOWARN"
    fi
  fi
  CFLAGS=`echo $CFLAGS | sed 's/^ *//' | sed 's/ *$//' | sed 's/  */ /'`
  AC_LANG_SAVE
  AC_LANG_C
  AC_MSG_CHECKING([whether $CC accepts $CFLAGS])
  AC_TRY_COMPILE([],[],
    AC_MSG_RESULT(yes),
    [
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([compiler flags check failed.  Please set CFLAGS explicitly.])
    ]
  )
  AC_LANG_RESTORE

  if test -n "$save_CXXFLAGS"; then
    CXXFLAGS="$save_CXXFLAGS"
  else
    CXXFLAGS="$try_CXXFLAGS"
    if test "$ac_cv_compiler_optimization" = yes; then
      CXXFLAGS="$CXXFLAGS $try_CXXFLAGS_OPT"
    else
      CXXFLAGS="$CXXFLAGS $try_CXXFLAGS_DEBUG"
    fi
    if test "$ac_cv_compiler_exceptions" = yes; then
      CXXFLAGS="$CXXFLAGS $try_CXXFLAGS_EH"
    else
      CXXFLAGS="$CXXFLAGS $try_CXXFLAGS_NOEH"
    fi
    if test "$ac_cv_compiler_warnings" = yes; then
      CXXFLAGS="$CXXFLAGS $try_CXXFLAGS_WARN"
    else
      CXXFLAGS="$CXXFLAGS $try_CXXFLAGS_NOWARN"
    fi
  fi
  CXXFLAGS=`echo $CXXFLAGS | sed 's/^ *//' | sed 's/ *$//' | sed 's/  */ /'`
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_MSG_CHECKING([whether $CXX accepts $CXXFLAGS])
  AC_TRY_COMPILE([],[],
    AC_MSG_RESULT(yes),
    [
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([compiler flags check failed.  Please set CXXFLAGS explicitly.])
    ]
  )
  AC_LANG_RESTORE
  ]

  ac_cv_compiler="$COMPILER"
  ac_cv_compiler_cc="$CC"
  ac_cv_compiler_cflags="$CFLAGS"
  ac_cv_compiler_cxx="$CXX"
  ac_cv_compiler_cxxflags="$CXXFLAGS"

  ac_cv_prog_ac_ct_CC="$CC"
  ac_cv_prog_ac_ct_CXX="$CXX"
)
