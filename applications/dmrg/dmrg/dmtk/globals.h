#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#ifndef MIN_VECTOR_SIZE
#define MIN_VECTOR_SIZE 1000
#endif
#ifndef MIN_MATRIX_SIZE
#define MIN_MATRIX_SIZE 100
#endif

#include "vector.h"
#include "matrix.h"

namespace dmtk
{

template<class T>
class DMTKglobals
{
  public:
    Matrix<T> m1;
    Matrix<T> m2;
    Matrix<T> m3;
    Matrix<T> m4;
    Vector<T> v1;
    Vector<T> v2;
    Vector<T> v3;
    Vector<T> v4;

    DMTKglobals() { m1 = m2 = m3 = m4 = Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE); v1 = v2 = v3 = v4 = Vector<T>(MIN_VECTOR_SIZE); }
};

static DMTKglobals<double> globals_double;
static DMTKglobals<complex<double> > globals_complex;

DMTKglobals<double> &
get_globals(const double&)
{
  return globals_double;
}

DMTKglobals<complex<double> >&
get_globals(const complex<double>&)
{
  return globals_complex;
}

} //namespace dmtk

#endif // __GLOBALS_H__
