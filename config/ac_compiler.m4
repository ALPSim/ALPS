AC_DEFUN([AC_COMPILER],
  [
  AC_SUBST(COMPILER)
  AC_SUBST(CPPFLAGS)
  AC_SUBST(CFLAGS)
  AC_SUBST(CXXFLAGS)

  AC_MSG_CHECKING([for exception handling])
  AC_ARG_ENABLE(exceptions,
    AC_HELP_STRING([--disable-exceptions],
      [disable exception handling]),
    [
    if test "x$enableval" = "xno"; then
      ac_cv_compiler_exceptions=no
    fi
    ]
  )
  test -z "$ac_cv_compiler_exceptions" && ac_cv_compiler_exceptions=yes
  AC_MSG_RESULT($ac_cv_compiler_exceptions)

  AC_MSG_CHECKING([for optimization])
  AC_ARG_ENABLE(optimization,
    AC_HELP_STRING([--disable-optimization],
      [disable optimization flags]),
    [
    if test "x$enableval" = "xno"; then
      ac_cv_compiler_optimization=no
    fi
    ]
  )
  test -z "$ac_cv_compiler_optimization" && ac_cv_compiler_optimization=yes
  AC_MSG_RESULT($ac_cv_compiler_optimization)

  AC_ARG_WITH(compiler,
    AC_HELP_STRING([--with-compiler=MODE],
      [set compiler mode (MODE = gnu, kai, intel, intel64, como, hp, dec, sgi32, sgi64, cray, macos, generic)]),
    [
    case "x$withval" in
      xgnu* | xGNU* | xgcc* | xg++* )
        COMPILER="gnu"
        ;;
      xkai* | xKAI*)
        COMPILER="kai"
        ;;
      xintel32* | xintel | xicc32* | xicc)
        COMPILER="intel"
        ;;
      xintel64* | xicc64*)
        COMPILER="intel64"
        ;;
      xcomo*)
        COMPILER="como"
        ;;
      xhp32* | xaCC32*)
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
      xcray* | xCRAY*)
        COMPILER="cray"
        ;;
      xmacos* | xmac* | xosx*)
        COMPILER="macos"
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
    kai)
      try_CC="KCC"
      try_CXX="KCC"
      ;;
    intel)
      try_CC="icc"
      try_CXX="icc"
      ;;
    intel64)
      try_CC="ecc"
      try_CXX="ecc"
      ;;
    como)
      try_CC="como"
      try_CXX="como"
      ;;
    hp32)
      try_CC="cc"
      try_CXX="aCC"
      ;;
    hp64)
      try_CC="cc"
      try_CXX="aCC"
      ;;
    dec)
      try_CC="cc"
      try_CXX="cxx"
      ;;
    sgi32)
      try_CC="cc"
      try_CXX="CC"
      ;;
    sgi64)
      try_CC="cc"
      try_CXX="CC"
      ;;
    cray)
      try_CC="cc"
      try_CXX="CC"
      ;;
    macos)
      try_CC="cc"
      try_CXX="c++"
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

  case "$COMPILER" in
    gnu)
      try_CFLAGS_OPT="-pthread -w -O3"
      try_CFLAGS_DEBUG="-pthread -W -Wall -Wno-sign-compare -g -O0"
      try_CXXFLAGS_OPT="-pthread -w -ftemplate-depth-45 -O3"
      try_CXXFLAGS_DEBUG="-pthread -W -Wall -Wno-sign-compare -ftemplate-depth-45 -g -O0"
      try_CXXFLAGS_EH="-fexceptions"
      try_CXXFLAGS_NOEH="-fno-exceptions"
      ;;
    kai)
      try_CFLAGS_OPT="--c --thread_safe -w +K3 -O3"
      try_CFLAGS_DEBUG="--c --thread_safe -g +K0 -O0"
      try_CXXFLAGS_OPT="--restrict --one_instantiation_per_object --thread_safe -w +K3 -O3 -DBOOST_REGEX_NO_EXTERNAL_TEMPLATES"
      try_CXXFLAGS_DEBUG="--restrict --one_instantiation_per_object --thread_safe --display_error_number -g +K0 -O0 -DBOOST_REGEX_NO_EXTERNAL_TEMPLATES"
      try_CXXFLAGS_EH="--exceptions"
      try_CXXFLAGS_NOEH="--no_exceptions"
      ;;
    intel)
      try_CFLAGS_OPT="-w -O3"
      try_CFLAGS_DEBUG=" -O0"
      try_CXXFLAGS_OPT="-D_REENTRANT -restrict -w -O3"
      try_CXXFLAGS_DEBUG="-D_REENTRANT -restrict -g -O0"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="-DBOOST_NO_EXCEPTIONS"
      ;;
    intel64)
      try_CFLAGS_OPT="-w -O3"
      try_CFLAGS_DEBUG=" -O0"
      try_CXXFLAGS_OPT="-D_REENTRANT -restrict -w -O3"
      try_CXXFLAGS_DEBUG="-D_REENTRANT -restrict -g -O0"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="-DBOOST_NO_EXCEPTIONS"
      ;;
    como)
      try_CFLAGS_OPT="-O3"
      try_CFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS_OPT="-O3"
      try_CXXFLAGS_DEBUG="-g -O0"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH=
      ;;
    hp32)
      try_CFLAGS_OPT="-Aa +DA1.1 -O3 -D_HPUX_SOURCE"
      try_CFLAGS_DEBUG="-Aa +DA1.1 -g -O0 -D_HPUX_SOURCE"
      try_CXXFLAGS_OPT="-Aa +DA1.1 -O3 -D_HPUX_SOURCE"
      try_CXXFLAGS_DEBUG="-Aa +DA1.1 -g -O0 -D_HPUX_SOURCE"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="+noeh"
      ;;
    hp64)
      try_CFLAGS_OPT="-Aa +DA2.0W -O3 -D_HPUX_SOURCE"
      try_CFLAGS_DEBUG="-Aa +DA2.0W -g -O0 -D_HPUX_SOURCE"
      try_CXXFLAGS_OPT="-Aa +DA2.0W -O3 -D_HPUX_SOURCE"
      try_CXXFLAGS_DEBUG="-Aa +DA2.0W -g -O0 -D_HPUX_SOURCE"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="+noeh"
      ;;
    dec)
      try_CFLAGS_OPT="-pthread -O4"
      try_CFLAGS_DEBUG="-pthread -g -O0"
      try_CXXFLAGS_OPT="-std ansi -pthread -D__USE_STD_IOSTREAM -O4"
      try_CXXFLAGS_DEBUG="-std ansi -pthread -D__USE_STD_IOSTREAM -g -O0"
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="-noexceptions"
      ;;
    sgi32)
      # try_CFLAGS_OPT="-n32 -D_POSIX_C_SOURCE=199506L -D_BSD_COMPAT -Ofast -INLINE"
      # try_CFLAGS_DEBUG="-n32 -D_POSIX_C_SOURCE=199506L -D_BSD_COMPAT -g -O0"
      # try_CXXFLAGS_OPT="-n32 -LANG:std -D_POSIX_C_SOURCE=199506L -D_BSD_COMPAT -w -Ofast -INLINE"
      # try_CXXFLAGS_DEBUG="-n32 -LANG:std -D_POSIX_C_SOURCE=199506L -D_BSD_COMPAT -g -O0"
      try_CFLAGS_OPT="-n32 -Ofast -INLINE"
      try_CFLAGS_DEBUG="-n32 -g -O0"
      try_CXXFLAGS_OPT="-n32 -LANG:std -w -Ofast -INLINE"
      try_CXXFLAGS_DEBUG="-n32 -LANG:std -g -O0"
      try_CXXFLAGS_EH="-LANG:exceptions=ON"
      try_CXXFLAGS_NOEH="-LANG:exceptions=OFF"
      ;;
    sgi64)
      # try_CFLAGS_OPT="-64 -D_POSIX_C_SOURCE=199506L -D_BSD_COMPAT -Ofast -INLINE"
      # try_CFLAGS_DEBUG="-64 -D_POSIX_C_SOURCE=199506L -D_BSD_COMPAT -g -O0"
      # try_CXXFLAGS_OPT="-64 -LANG:std -D_POSIX_C_SOURCE=199506L -D_BSD_COMPAT -w -Ofast -INLINE"
      # try_CXXFLAGS_DEBUG="-64 -LANG:std -D_POSIX_C_SOURCE=199506L -D_BSD_COMPAT -g -O0"
      try_CFLAGS_OPT="-64 -Ofast -INLINE"
      try_CFLAGS_DEBUG="-64 -g -O0"
      try_CXXFLAGS_OPT="-64 -LANG:std -w -Ofast -INLINE"
      try_CXXFLAGS_DEBUG="-64 -LANG:std -g -O0"
      try_CXXFLAGS_EH="-LANG:exceptions=ON"
      try_CXXFLAGS_NOEH="-LANG:exceptions=OFF"
      ;;
    cray)
      try_CFLAGS_OPT="-O2 -h conform"
      try_CFLAGS_DEBUG="-g -O0 -h conform"
      try_CXXFLAGS_OPT="-O2 -h one_instantiation_per_object -h new_for_init -h nodep_name -h parse_templates"
      try_CXXFLAGS_DEBUG="-g -O0 -h one_instantiation_per_object -h new_for_init -h nodep_name -h parse_templates"
      try_CXXFLAGS_EH="-h exceptions"
      try_CXXFLAGS_NOEH="-h noexceptions"
      ;;
    macos)
      #try_CFLAGS_OPT="-w -pthread -O3"
      #try_CFLAGS_DEBUG="-pthread -g -O0"
      #try_CXXFLAGS_OPT="-w -pthread -ftemplate-depth-45 -O3"
      #try_CXXFLAGS_DEBUG="-pthread -ftemplate-depth-45 -g -O0"
      try_CFLAGS_OPT="-w -O3 -Wno-long-double"
      try_CFLAGS_DEBUG="-g -O0 -Wno-long-double"
      try_CXXFLAGS_OPT="-w -ftemplate-depth-45 -O3 -Wno-long-double"
      try_CXXFLAGS_DEBUG="-ftemplate-depth-45 -g -O0 -Wno-long-double"
      try_CXXFLAGS_EH="-fexceptions"
      try_CXXFLAGS_NOEH="-fno-exceptions"
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
      try_CXXFLAGS_EH=
      try_CXXFLAGS_NOEH="-DBOOST_NO_EXCEPTIONS"
      ;;
    *)
      AC_MSG_ERROR([unknown mode $COMPILER])
      ;;
  esac

  if test -z "$CPPFLAGS"; then
    if test "$ac_cv_compiler_optimization" = yes; then
      CPPFLAGS="-DNDEBUG"
    fi
  fi

  if test -n "$save_CFLAGS"; then
    CFLAGS="$save_CFLAGS"
  else
    if test "$ac_cv_compiler_optimization" = yes; then
      CFLAGS="$try_CFLAGS_OPT"
    else
      CFLAGS="$try_CFLAGS_DEBUG"
    fi
  fi

  if test "$CFLAGS" != "$save_CFLAGS"; then
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
  fi

  if test -n "$save_CXXFLAGS"; then
    CXXFLAGS="$save_CXXFLAGS"
  else
    if test "$ac_cv_compiler_optimization" = yes; then
      if test "$ac_cv_compiler_exceptions" = yes; then
        if test -n "$try_CXXFLAGS_OPT"; then
          if test -n "$try_CXXFLAGS_EH"; then
            CXXFLAGS="$try_CXXFLAGS_OPT $try_CXXFLAGS_EH"
          else
            CXXFLAGS="$try_CXXFLAGS_OPT"
          fi
        else
          if test -n "$try_CXXFLAGS_EH"; then
            CXXFLAGS="$try_CXXFLAGS_EH"
          else
            CXXFLAGS=
          fi
        fi
      else
        if test -n "$try_CXXFLAGS_OPT"; then
          if test -n "$try_CXXFLAGS_NOEH"; then
            CXXFLAGS="$try_CXXFLAGS_OPT $try_CXXFLAGS_NOEH"
          else
            CXXFLAGS="$try_CXXFLAGS_OPT"
          fi
        else
          if test -n "$try_CXXFLAGS_NOEH"; then
            CXXFLAGS="$try_CXXFLAGS_NOEH"
          else
            CXXFLAGS=
          fi
        fi
      fi
    else
      if test "$ac_cv_compiler_exceptions" = yes; then
        if test -n "$try_CXXFLAGS_DEBUG"; then
          if test -n "$try_CXXFLAGS_EH"; then
            CXXFLAGS="$try_CXXFLAGS_DEBUG $try_CXXFLAGS_EH"
          else
            CXXFLAGS="$try_CXXFLAGS_DEBUG"
          fi
        else
          if test -n "$try_CXXFLAGS_EH"; then
            CXXFLAGS="$try_CXXFLAGS_EH"
          else
            CXXFLAGS=
          fi
        fi
      else
        if test -n "$try_CXXFLAGS_DEBUG"; then
          if test -n "$try_CXXFLAGS_NOEH"; then
            CXXFLAGS="$try_CXXFLAGS_DEBUG $try_CXXFLAGS_NOEH"
          else
            CXXFLAGS="$try_CXXFLAGS_DEBUG"
          fi
        else
          if test -n "$try_CXXFLAGS_NOEH"; then
            CXXFLAGS="$try_CXXFLAGS_NOEH"
          else
            CXXFLAGS=
          fi
        fi
      fi
    fi
  fi

  if test "$CXXFLAGS" != "$save_CXXFLAGS"; then
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    AC_MSG_CHECKING([whether $CXX accepts $CXXFLAGS])
    AC_TRY_COMPILE([],[],
      AC_MSG_RESULT(yes),
      [
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([compiler flags check failed.  Please set CXXFLAGS and explicitly.])
      ]
    )
    AC_LANG_RESTORE
    ]
  fi

  ac_cv_compiler="$COMPILER"
  ac_cv_compiler_cc="$CC"
  ac_cv_compiler_cflags="$CFLAGS"
  ac_cv_compiler_cxx="$CXX"
  ac_cv_compiler_cxxflags="$CXXFLAGS"
)
