AC_DEFUN([AC_LAPACK],
  [
  AC_REQUIRE([AC_CANONICAL_HOST])

  AC_SUBST(LAPACK_CPPFLAGS)
  AC_SUBST(LAPACK_LDFLAGS)
  AC_SUBST(LAPACK_LIBS)
    
  AC_ARG_WITH(atlas,
    AC_HELP_STRING([--with-atlas=LIBS],[use ATLAS Library]),
    [
    if test "x$withval" = xno; then
      atlas=no
    else
      atlas=yes
      if test "x$withval" != xyes; then
        atlas_flags=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )
  AC_ARG_WITH(atlas-dir,
    AC_HELP_STRING([--with-atlas-dir=DIR],[ATLAS lib directory]),
    [atlas_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`]
  )
  
  AC_ARG_WITH(blas,
    AC_HELP_STRING([--with-blas=LIBS],[use BLAS Library]),
    [
    if test "x$withval" = xno; then
      blas=no
    else
      blas=yes
      if test "x$withval" != xyes; then
        blas_flags=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )
  AC_ARG_WITH(blas-dir,
    AC_HELP_STRING([--with-blas-dir=DIR],[BLAS lib directory]),
    [blas_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`]
  )

  AC_ARG_WITH(dxml,
    AC_HELP_STRING([--with-dxml], [use Digital Extended Math Library]),
    [
    if test "x$withval" = xno; then
      dxml=no
    else
      dxml=yes
    fi
    ]
  )

  AC_ARG_WITH(lapack,
    AC_HELP_STRING([--with-lapack=LIBS],[use LAPACK library]),
    [
    if test "x$withval" = xno; then
      lapack=no
    else
      lapack=yes
      if test "x$withval" != xyes; then
        lapack_flags=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`
      fi
    fi
    ]
  )
  AC_ARG_WITH(lapack-dir,
    AC_HELP_STRING([--with-lapack-dir=DIR],[LAPACK lib directory]),
    [lapack_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`]
  )

  AC_ARG_WITH(mkl,
    AC_HELP_STRING([--with-mkl=PROC],
    [use Intel Math Kernel Library (PROC=p3, p4, itp)]),
    [
    if test "x$withval" = xno; then
      mkl=no
    else
      case "x$withval" in
        xp3)
          mkl=yes
          mkl_proc="p3"
          ;;
        xp4)
          mkl=yes
          mkl_proc="p4"
          ;;
        xitp)
          mkl=yes
          mkl_proc="itp"
          ;;
        *)
          mkl=yes
          mkl_proc=
      esac
    fi
    ]
  )
  AC_ARG_WITH(mkl-dir,
    AC_HELP_STRING([--with-mkl-dir=DIR],[MKL lib directory]),
    [mkl_dir=`echo "$withval" | sed 's,//*,/,g' | sed 's,/$,,'`]
  )

  AC_ARG_WITH(scsl,
    AC_HELP_STRING([--with-scsl], [use SGI SCSL Scientific Library]),
    [
    if test "x$withval" = xno; then
      scsl=no
    else
      scsl=yes
    fi
    ]
  )

  AC_ARG_WITH(essl,
    AC_HELP_STRING([--with-essl], [use ESSL library on IBM AIX]),
    [
    if test "x$withval" = xno; then
      essl=no
    else
      essl=yes
    fi
    ]
  )

  AC_ARG_WITH(veclib,
    AC_HELP_STRING([--with-veclib], [use vecLib framework on Mac OS X]),
    [
    if test "x$withval" = xno; then
      veclib=no
    else
      veclib=yes
    fi
    ]
  )

  if test -n "$atlas_dir"; then
    if test "x$atlas" = xno; then
      AC_MSG_ERROR([inconsistent options for ATLAS library.])
    fi
    atlas=yes
  fi

  if test -n "$blas_dir"; then
    if test "x$blas" = xno; then
      AC_MSG_ERROR([inconsistent options for BLAS library.])
    fi
    blas=yes
  fi

  if test -n "$lapack_dir"; then
    if test "x$lapack" = xno; then
      AC_MSG_ERROR([inconsistent options for LAPACK library.])
    fi
    lapack=yes
  fi

  if test -n "$mkl_dir"; then
    if test "x$mkl" = xno; then
      AC_MSG_ERROR([inconsistent options for MKL library.])
    fi
    mkl=yes
  fi

  if test "x$atlas" = xyes; then
    if test "x$blas" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$dxml" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$mkl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$scsl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$essl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$veclib" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
  fi
  if test "x$blas" = xyes; then
    if test "x$dxml" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$mkl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$scsl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$veclib" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
  fi
  if test "x$dxml" = xyes; then
    if test "x$mkl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$scsl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$veclib" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$essl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
  fi
  if test "x$mkl" = xyes; then
    if test "x$scsl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$veclib" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$essl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
  fi
  if test "x$scsl" = xyes; then
    if test "x$veclib" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
    if test "x$essl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
  fi
  if test "x$veclib" = xyes; then
    if test "x$essl" = xyes; then
      AC_MSG_ERROR([more than one BLAS-like libraries specified.])
    fi
  fi
    
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_LIBS="$LIBS"

  found_blas=no
  found_lapack=no

  # for Math Kernel Library
  if test "$found_blas" = no; then
    if test "x$mkl" != xno; then

      AC_MSG_NOTICE([checking for Math Kernel Library])
      AC_MSG_CHECKING([MKL directory])
      if test -z "$mkl_dir"; then
        for d in mkl721 mkl72 mkl70 mkl61 mkl; do
          if test -d "/opt/intel/$d"; then
            mkl_dir="/opt/intel/$d"
            AC_MSG_RESULT([$mkl_dir])
            break;
          fi
        done
        if test -z "$mkl_dir"; then
          AC_MSG_RESULT([not found])
        fi
      else
        if test -d "$mkl_dir"; then
          AC_MSG_RESULT([$mkl_dir])
        else
          AC_MSG_RESULT([$not found])
          AC_MSG_ERROR([MKL directory not found.])
        fi
      fi

      if test -z "$mkl_proc"; then
        case "x$host_cpu" in
        xi686)
          mkl_proc="p4"
          ;;
        xx86_64)
          mkl_proc="x86_64"
          ;;
        xia64)
          mkl_proc="itp"
          ;;
        *)
          ;;
        esac
      fi

      AC_MSG_CHECKING([for processor type])
      case "x$mkl_proc" in
      xp3)
        mkl_procname="Pentium 3"
        mkl_bits="32"
        mkl_lib="mkl_p3"
        AC_MSG_RESULT([$mkl_procname])
        ;;
      xp4)
        mkl_procname="Pentium 4"
        mkl_bits="32"
        mkl_lib="mkl_p4"
        AC_MSG_RESULT([$mkl_procname])
        ;;
      xx86_64)
        mkl_procname="Pentium 4 with EM64T"
        mkl_bits="em64t"
        mkl_lib="mkl_em64t"
        AC_MSG_RESULT([$mkl_procname])
        ;;
      xitp)
        mkl_procname="Itanium"
        mkl_bits="64"
        mkl_lib="mkl_itp"
        AC_MSG_RESULT([$mkl_procname])
        ;;
      *)
        mkl_proc=
        AC_MSG_RESULT([unknown processor type])
      esac

      if test -n "$mkl_proc"; then
        if test "$mkl_bits" = "32"; then
          if test -f "$mkl_dir/lib/$mkl_bits/libmkl_ia32.a"; then
            # version 6.1 & 7
            mkl_lib="mkl_ia32"
          fi
        fi
        if test "$mkl_bits" = "64"; then
          if test -f "$mkl_dir/lib/$mkl_bits/libmkl_ipf.a"; then
            # version 6.1 & 7
            mkl_lib="mkl_ipf"
          fi
        fi

        if test -d "$mkl_dir/lib/$mkl_bits"; then
          mkl_ldflags="-L$mkl_dir/lib/$mkl_bits"
        else
          mkl_ldflags=
        fi
        LDFLAGS="$mkl_ldflags $ac_save_LDFLAGS"

        mkl_libs=
        found=no

        if test "$mkl_proc" != "itp"; then
          AC_CHECK_FUNC(d_abs,,
            [AC_CHECK_LIB(g2c, d_abs, 
              [mkl_libs="-lg2c"; LIBS="$mkl_libs $ac_save_LIBS"],
              [AC_CHECK_LIB(F90, d_abs, 
                [mkl_libs="-lF90 $mkl_libs"; LIBS="$mkl_libs $ac_save_LIBS"])
              ])
            ]
          )
          AC_CHECK_LIB(guide, omp_in_parallel_, 
            [mkl_libs="-lguide $mkl_libs"
             LIBS="$mkl_libs $ac_save_LIBS"
             AC_CHECK_LIB($mkl_lib, dgemm_, 
              [mkl_libs="-l$mkl_lib $mkl_libs"
               LIBS="$mkl_libs $ac_save_LIBS"
               AC_CHECK_LIB(mkl_lapack, dsyev_,
                [mkl_libs="-lmkl_lapack $mkl_libs"; found=yes])
              ])
            ]
          )
        else
          AC_CHECK_FUNC(d_abs,,
            [AC_CHECK_LIB(g2c, d_abs, 
              [mkl_libs="-lg2c $mkl_libs"; LIBS="$mkl_libs $ac_save_LIBS"],
              [AC_CHECK_LIB(F90, d_abs, 
                [mkl_libs="-lF90 $mkl_libs"; LIBS="$mkl_libs $ac_save_LIBS"])
              ])
            ]
          )
          AC_CHECK_FUNC(f_concat,,
            [AC_CHECK_LIB(IEPCF90, f_concat, 
              [mkl_libs="-lIEPCF90 $mkl_libs"; LIBS="$mkl_libs $ac_save_LIBS"])
            ]
          )
          AC_CHECK_LIB(guide, omp_in_parallel_, 
            [mkl_libs="-lguide $mkl_libs"
             LIBS="$mkl_libs $ac_save_LIBS"
             AC_CHECK_LIB($mkl_lib, dgemm_, 
              [mkl_libs="-l$mkl_lib $mkl_libs"
               LIBS="$mkl_libs $ac_save_LIBS"
               AC_CHECK_LIB(mkl_lapack, dsyev_,
                [mkl_libs="-lmkl_lapack $mkl_libs"; found=yes])
              ])
            ]
          )
        fi
        if test "$found" = yes; then
          LAPACK_LDFLAGS="$mkl_ldflags"
          LAPACK_LIBS="$mkl_libs"
          found_blas=yes
          found_lapack=yes
        fi
      fi
    fi
    if test "x$mkl" = xyes; then
      if test "$found_blas" = no; then
        AC_MSG_ERROR([Math Kernel Library not found.])
      fi
    fi
  fi

  # for Digital Extended Math Library
  if test "$found_blas" = no; then
    if test "x$dxml" != xno; then
      AC_MSG_NOTICE([checking for Digital Extended Math Library])
      LDFLAGS="$ac_save_LDFLAGS"
      LIBS="$ac_save_LIBS"
      AC_CHECK_LIB(dxml, dgemm_, 
        [AC_CHECK_LIB(dxml, dsyev_,
          [LAPACK_LDFLAGS=; LAPACK_LIBS="-ldxml";
           found_blas=yes; found_lapack=yes])
        ]
      )
    fi
    if test "x$dxml" = xyes; then
      if test "$found_blas" = no; then
        AC_MSG_ERROR([Digital Extended Math Library not found.])
      fi
    fi
  fi

  # for SGI SCSL
  if test "$found_blas" = no; then
    if test "x$scsl" != xno; then
      AC_MSG_NOTICE([checking for SGI SCSL Library])
      LDFLAGS="$ac_save_LDFLAGS"
      LIBS="$ac_save_LIBS"
      AC_CHECK_LIB(scs, dgemm_, 
        [AC_CHECK_LIB(scs, dsyev_,
          [LAPACK_LDFLAGS=; LAPACK_LIBS="-lscs";
           found_blas=yes; found_lapack=yes])
        ]
      )
    fi
    if test "x$scsl" = xyes; then
      if test "$found_blas" = no; then
        AC_MSG_ERROR([SCSL Library not found.])
      fi
    fi
  fi

  # for vecLib on Mac OS X
  if test "$found_blas" = no; then
    if test "x$veclib" != xno; then
      AC_MSG_NOTICE([checking for vecLib framework on Mac OS X])
      veclib_ldflags="-faltivec -framework vecLib"
      LDFLAGS="$veclib_ldflags $ac_save_LDFLAGS"
      LIBS="$ac_save_LIBS"
      AC_MSG_CHECKING([for dgemm_ in $veclib_ldflags])
      AC_TRY_LINK([extern "C" char dgemm_();],[dgemm_();],
        [AC_MSG_RESULT(yes)
         AC_MSG_CHECKING([for dsyev_ in $veclib_ldflags])
         AC_TRY_LINK([extern "C" char dsyev_();],[dsyev_();],
                     [AC_MSG_RESULT(yes)
                      LAPACK_LDFLAGS="$veclib_ldflags"; LAPACK_LIBS=
                      found_blas=yes; found_lapack=yes],
                     [AC_MSG_RESULT(no)])],
        [AC_MSG_RESULT(no)]
      )
    fi
    if test "x$veclib" = xyes; then
      if test "$found_blas" = no; then
        AC_MSG_ERROR([Mac OS X vecLib framework not found.])
      fi
    fi
  fi

  # for ESSL on IBM AIX
  if test "$found_blas" = no; then
    if test "x$essl" != xno; then
      AC_MSG_NOTICE([checking for ESSL on IBM AIX])
      LDFLAGS="$ac_save_LDFLAGS"
      LIBS="$ac_save_LIBS"
      AC_CHECK_LIB(essl, dgemm_, 
          [LAPACK_LDFLAGS=; LAPACK_LIBS="-lessl";
           found_blas=yes;])
    fi
    if test "x$essl" = xyes; then
      if test "$found_blas" = no; then
        AC_MSG_ERROR([ESSL not found.])
      fi
    fi
  fi

  # for ATLAS
  if test "$found_blas" = no; then
    if test "x$atlas" != xno; then
      AC_MSG_NOTICE([checking for ATLAS Library])
      found=no
      if test -n "$atlas_dir"; then
        AC_MSG_CHECKING([for ATLAS library directory])
        if test -d "$atlas_dir"; then 
          AC_MSG_RESULT([$atlas_dir])
          atlas_ldflags="-L$atlas_dir"
        else
          AC_MSG_RESULT([no])
          AC_MSG_ERROR([$atlas_dir not found.])
        fi
      fi
      if test -n "$atlas_flags"; then
        atlas_libs="$atlas_flags"
        LDFLAGS="$atlas_ldflags $ac_save_LDFLAGS"
        LIBS="$atlas_libs $ac_save_LIBS"
        AC_MSG_CHECKING([for dgemm_ in $atlas_ldflags $atlas_libs])
        AC_TRY_LINK([extern "C" char dgemm_();],[dgemm_();],
                    [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
      fi
      if test "$found" = no; then
        for atlas_libs in '-latlas' '-lf77blas -latlas' '-lf77blas -latlas -lg2c' '-lf77blas -latlas -lI77 -lF77'; do
          LDFLAGS="$atlas_ldflags $ac_save_LDFLAGS"
          LIBS="$atlas_libs $ac_save_LIBS"
          AC_MSG_CHECKING([for dgemm_ in $atlas_ldflags $atlas_libs])
          AC_TRY_LINK([extern "C" char dgemm_();],[dgemm_();],
                      [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
          if test "$found" = no; then
            if test -z "$atlas_dir"; then
              for d in $HOME $HOME/src $prefix /usr/local /usr/local/src; do
                if test -d "$d/lib"; then
                  atlas_ldflags="-L$d/lib"
                  LDFLAGS="$atlas_ldflags $ac_save_LDFLAGS"
                  AC_MSG_CHECKING([for dgemm_ in $atlas_ldflags $atlas_libs])
                  AC_TRY_LINK([extern "C" char dgemm_();],[dgemm_();],
                              [AC_MSG_RESULT(yes); found=yes],
                              [AC_MSG_RESULT(no)])
                  if test "$found" = yes; then
                    break
                  fi
                  atlas_ldflags=
                fi
              done
            fi
          fi
          if test "$found" = yes; then
            break
          fi
        done
      fi
    
      if test "$found" = yes; then
        found_blas=yes
        LAPACK_LDFLAGS="$atlas_ldflags"
        LAPACK_LIBS="$atlas_libs"
      fi
    fi
  fi 

  # for BLAS
  if test "$found_blas" = no; then
    if test "x$blas" != xno ; then
      AC_MSG_NOTICE([checking for BLAS Library])
      found=no
      if test -n "$blas_dir"; then
        AC_MSG_CHECKING([for BLAS library directory])
        if test -d "$blas_dir"; then 
          AC_MSG_RESULT([$blas_dir])
          blas_ldflags="-L$blas_dir"
        else
          AC_MSG_RESULT([no])
          AC_MSG_ERROR([$blas_dir not found.])
        fi
      fi
      if test -n "$blas_flags"; then
        blas_libs="$blas_flags"
        LDFLAGS="$blas_ldflags $ac_save_LDFLAGS"
        LIBS="$blas_libs $ac_save_LIBS"
        AC_MSG_CHECKING([for dgemm_ in $blas_ldflags $blas_libs])
        AC_TRY_LINK([extern "C" char dgemm_();],[dgemm_();],
                    [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
      fi
      if test "$found" = no ;  then
        for blas_libs in '-lsci' '-lblas' '-lblas -lf2c' '-lblas -lg2c' '-lblas -lpgftnrtl -lpgc'; do
          LDFLAGS="$blas_ldflags $ac_save_LDFLAGS"
          LIBS="$blas_libs $ac_save_LIBS"
          AC_MSG_CHECKING([for dgemm_ in $blas_ldflags $blas_libs])
          AC_TRY_LINK([extern "C" char dgemm_();],[dgemm_();],
                      [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
          if test "$found" = no; then
            if test -z "$blas_dir"; then
              for d in $HOME $HOME/src $prefix /usr/local /usr/local/src; do
                if test -d "$d/lib"; then
                  blas_ldflags="-L$d/lib"
                  LDFLAGS="$blas_ldflags $ac_save_LDFLAGS"
                  AC_MSG_CHECKING([for dgemm_ in $blas_ldflags $blas_libs])
                  AC_TRY_LINK([extern "C" char dgemm_();],[dgemm_();],
                              [AC_MSG_RESULT(yes); found=yes],
                              [AC_MSG_RESULT(no)])
                  if test "$found" = yes; then
                    break
                  fi
                  blas_ldflags=
                fi
              done
            fi
          fi
          if test "$found" = yes; then
            break
          fi
        done
      fi
    
      if test "$found" = yes; then
        found_blas=yes
        LAPACK_LDFLAGS="$blas_ldflags"
        LAPACK_LIBS="$blas_libs"
      fi
    fi
  fi

  # for LAPACK
  if test "$found_blas" = yes; then
    if test "$found_lapack" = no; then
      if test "$lapack" != no ; then
        AC_MSG_NOTICE([checking for LAPACK Library])
        found=no
        if test -n "$lapack_dir"; then
          AC_MSG_CHECKING([for LAPACK library directory])
          if test -d "$lapack_dir"; then 
            AC_MSG_RESULT([$lapack_dir])
            lapack_ldflags="-L$lapack_dir"
          else
            AC_MSG_RESULT([no])
            AC_MSG_ERROR([$lapack_dir not found.])
          fi
        fi
        if test -n "$lapack_flags"; then
          lapack_libs="$lapack_flags"
          LDFLAGS="$lapack_ldflags $LAPACK_LDFLAGS $ac_save_LDFLAGS"
          LIBS="$lapack_libs $LAPACK_LIBS $ac_save_LIBS"
          if test "$ac_cv_compiler" = ibm32; then
	    AC_MSG_CHECKING([for dsyev in $lapack_ldflags $lapack_libs])
	    AC_TRY_LINK([extern "C" char dsyev();],[dsyev();],
                        [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
	  else
	    if test "$ac_cv_compiler" = ibm64; then
	      AC_MSG_CHECKING([for dsyev in $lapack_ldflags $lapack_libs])
	      AC_TRY_LINK([extern "C" char dsyev();],[dsyev();],
                          [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
	    else
              AC_MSG_CHECKING([for dsyev_ in $lapack_ldflags $lapack_libs])
              AC_TRY_LINK([extern "C" char dsyev_();],[dsyev_();],
                          [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
	    fi
	  fi
        fi
        if test "$found" = no; then
          for lapack_libs in '-llapack' '-llapack -lf2c' '-llapack -lg2c'; do
            LDFLAGS="$lapack_ldflags $LAPACK_LDFLAGS $ac_save_LDFLAGS"
            LIBS="$lapack_libs $LAPACK_LIBS $ac_save_LIBS"
            if test "$ac_cv_compiler" = ibm32; then
              AC_MSG_CHECKING([for dsyev in $lapack_ldflags $lapack_libs])
              AC_TRY_LINK([extern "C" char dsyev();],[dsyev();],
                          [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
	    else 
	      if test "$ac_cv_compiler" = ibm64; then
                AC_MSG_CHECKING([for dsyev in $lapack_ldflags $lapack_libs])
                AC_TRY_LINK([extern "C" char dsyev();],[dsyev();],
                            [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
	      else
                AC_MSG_CHECKING([for dsyev_ in $lapack_ldflags $lapack_libs])
                AC_TRY_LINK([extern "C" char dsyev_();],[dsyev_();],
                            [AC_MSG_RESULT(yes); found=yes],[AC_MSG_RESULT(no)])
	      fi
	    fi
            if test "$found" = no; then
              if test -z "$lapack_dir"; then
                for d in $HOME $HOME/src $prefix /usr/local /usr/local/src; do
                  if test -d "$d/lib"; then
                    lapack_ldflags="-L$d/lib"
                    LDFLAGS="$lapack_ldflags $LAPACK_LDFLAGS $ac_save_LDFLAGS"
                    if test "$ac_cv_compiler" = ibm32; then
                      AC_MSG_CHECKING([for dsyev in $lapack_ldflags $lapack_libs])
                      AC_TRY_LINK([extern "C" char dsyev();],[dsyev();],
                                  [AC_MSG_RESULT(yes); found=yes],
                                  [AC_MSG_RESULT(no)])
		    else
		      if test "$ac_cv_compiler" = ibm64; then
                        AC_MSG_CHECKING([for dsyev in $lapack_ldflags $lapack_libs])
                        AC_TRY_LINK([extern "C" char dsyev();],[dsyev();],
                                    [AC_MSG_RESULT(yes); found=yes],
                                    [AC_MSG_RESULT(no)])
		      else
                        AC_MSG_CHECKING([for dsyev_ in $lapack_ldflags $lapack_libs])
                        AC_TRY_LINK([extern "C" char dsyev_();],[dsyev_();],
                                    [AC_MSG_RESULT(yes); found=yes],
                                    [AC_MSG_RESULT(no)])
		      fi
		    fi
                    if test "$found" = yes; then
                      break
                    fi
                    lapack_ldflags=
                  fi
                done
              fi
            fi
            if test "$found" = yes; then
              break
            fi
          done
        fi
      
        if test "$found" = yes; then
          found_lapack=yes
          LAPACK_LDFLAGS="$lapack_ldflags $LAPACK_LDFLAGS"
          LAPACK_LIBS="$lapack_libs $LAPACK_LIBS"
        fi
      fi
    fi
  fi 

  LDFLAGS=$ac_save_LDFLAGS
  LIBS=$ac_save_LIBS
  AC_LANG_RESTORE

  if test "$found_blas" = yes; then
    ac_cv_have_blas=yes
    AC_MSG_NOTICE([enabling BLAS support])
    AC_DEFINE(HAVE_BLAS, [], [Define if you have a BLAS library.])
  else
    ac_cv_have_blas=no
    AC_MSG_NOTICE([disabling BLAS support])
  fi

  if test "$found_lapack" = yes; then
    ac_cv_have_lapack=yes
    ac_cv_lapack_cppflags="$LAPACK_CPPFLAGS"
    ac_cv_lapack_ldflags="$LAPACK_LDFLAGS"
    ac_cv_lapack_libs="$LAPACK_LIBS"
    AC_MSG_NOTICE([enabling LAPACK support])
    AC_DEFINE(HAVE_LAPACK, [], [Define if you have a LAPACK library.])
  else
    ac_cv_have_lapack=no
    AC_MSG_NOTICE([disabling LAPACK support])
  fi
  ]
)
