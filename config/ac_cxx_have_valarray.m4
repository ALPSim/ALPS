dnl @synopsis AC_CXX_HAVE_VALARRAY
dnl
dnl If the compiler has valarray<T>, define HAVE_VALARRAY.
dnl
dnl @author Synge Todo
dnl
AC_DEFUN([AC_CXX_HAVE_VALARRAY],
[AC_CACHE_CHECK(whether the compiler has valarray<T>,
ac_cv_cxx_have_valarray,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <valarray>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[valarray<float> x(100); return 0;],
 ac_cv_cxx_have_valarray=yes, ac_cv_cxx_have_valarray=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_valarray" = yes; then
  AC_DEFINE(HAVE_VALARRAY,,[define if the compiler has valarray<T>])
fi
])
