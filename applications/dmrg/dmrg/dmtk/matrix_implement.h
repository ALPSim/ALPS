#ifndef __DMTK_MATRIX_IMPLEMENT_H__
#define __DMTK_MATRIX_IMPLEMENT_H__

/////////////////////////////////////////////////////////////////

template<class T>
slice_iter<T>
Matrix<T>::operator()(Range col_range, size_t row)
{
  std::slice s(col_range.start()*num_rows + row, 
          col_range.size(),
          col_range.stride()*num_rows);

  slice_iter<T> iter(this, s);
  return iter; 
}

template<class T>
cslice_iter<T>
Matrix<T>::operator()(Range col_range, size_t row) const
{
  std::slice s(col_range.start()*num_rows + row, 
          col_range.size(),
          col_range.stride()*num_rows);

  cslice_iter<T> iter(this, s);
  return iter; 
}

template<class T>
slice_iter<T>
Matrix<T>::operator()(size_t col, Range row_range)
{
  std::slice s(col*num_rows+row_range.start(), 
          row_range.size(),
          row_range.stride());

  slice_iter<T> iter(this, s);
  return iter; 
}

template<class T>
cslice_iter<T>
Matrix<T>::operator()(size_t col, Range row_range) const
{
  std::slice s(col*num_rows+row_range.start(), 
          row_range.size(),
          row_range.stride());

  cslice_iter<T> iter(this, s);
  return iter; 
}

template<class T>
gslice_iter<T>
Matrix<T>::operator()(Range col_range, Range row_range)
{
  std::slice s_col(col_range.start()*num_rows, 
           col_range.size(),
           col_range.stride()*num_rows);
  std::slice s_row(row_range.start(), 
           row_range.size(),
           row_range.stride());
  gslice_iter<T> iter(this, s_col, s_row);
  return iter; 
}

template<class T>
cgslice_iter<T>
Matrix<T>::operator()(Range col_range, Range row_range) const
{
  std::slice s_col(col_range.start()*num_rows, 
           col_range.size(),
           col_range.stride()*num_rows);
  std::slice s_row(row_range.start(), 
           row_range.size(),
           row_range.stride());
  cgslice_iter<T> iter(this, s_col, s_row);
  return iter; 
}

template<class T>
slice_iter<T>
Matrix<T>::diagonal()
{
  std::slice s(0, std::min(num_rows,num_cols), num_rows+1);

  slice_iter<T> iter(this, s);
  return iter; 
}

template<class T>
cslice_iter<T>
Matrix<T>::diagonal() const
{
  std::slice s(0, std::min(num_rows,num_cols), num_rows+1);

  cslice_iter<T> iter(this, s);
  return iter; 
}

template<class T>
slice_iter<T>
Matrix<T>::as_vector()
{
  std::slice s(0, size1()*size2(), 1);

  slice_iter<T> iter(this, s);
  return iter; 
}

template<class T>
cslice_iter<T>
Matrix<T>::as_vector() const
{
  std::slice s(0, size1()*size2(), 1);

  cslice_iter<T> iter(this, s);
  return iter; 
}

template<class T>
slice_iter<T>
inline diagonal(Matrix<T> &m)
{
  return m.diagonal();
}

template<class T>
inline T trace(Matrix<T> &m) 
{ 
  return m.trace();
}

/////////////////////////////////////////////////////////////////


#define MATRIX_MATRIX_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, ConstRef<T,Matrix<T> >, ConstRef<T,Matrix<T> >, ap<T, ConstRef<T,Matrix<T> >, ConstRef<T,Matrix<T> > > > > \
op(Matrix<T>& a, Matrix<T>& b) \
{ \
  typedef ConstRef<T,Matrix<T> > IT; \
  typedef IterBinOp<T, IT, IT, ap<T, IT, IT> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b.ref())); \
}

MATRIX_MATRIX_OP(operator+, DMApAdd)
MATRIX_MATRIX_OP(operator-, DMApSubs)
MATRIX_MATRIX_OP(operator*, DMApMul)
MATRIX_MATRIX_OP(operator/, DMApDiv)
#undef MATRIX_MATRIX_OP

/////////////////////////////////////////////////////////////////

#define EXPR_MATRIX_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, IterExpr<T, Expr>, ConstRef<T,Matrix<T> >, ap<T, IterExpr<T, Expr>, ConstRef<T,Matrix<T> > > > > \
op(IterExpr<T, Expr> a, Matrix<T>& b) \
{ \
  typedef ConstRef<T,Matrix<T> > IT; \
  typedef IterBinOp<T, IterExpr<T, Expr>, IT, ap<T, IterExpr<T, Expr>, IT > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b.ref())); \
}

EXPR_MATRIX_OP(operator+, DMApAdd)
EXPR_MATRIX_OP(operator-, DMApSubs)
EXPR_MATRIX_OP(operator*, DMApMul)
EXPR_MATRIX_OP(operator/, DMApDiv)
#undef EXPR_MATRIX_OP

/////////////////////////////////////////////////////////////////

#define MATRIX_EXPR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, ConstRef<T,Matrix<T> >, IterExpr<T, Expr> , ap<T, ConstRef<T,Matrix<T> >, IterExpr<T, Expr> > > > \
op(Matrix<T>& a, IterExpr<T, Expr> b) \
{ \
  typedef ConstRef<T,Matrix<T> > IT; \
  typedef IterBinOp<T, IT, IterExpr<T, Expr>, ap<T, IT, IterExpr<T, Expr> > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b)); \
}

MATRIX_EXPR_OP(operator+, DMApAdd)
MATRIX_EXPR_OP(operator-, DMApSubs)
MATRIX_EXPR_OP(operator*, DMApMul)
MATRIX_EXPR_OP(operator/, DMApDiv)
#undef MATRIX_EXPR_OP

/////////////////////////////////////////////////////////////////

#define MATRIX_SCALAR_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, ConstRef<T,Matrix<T> >, T, ap<T, ConstRef<T,Matrix<T> >, T> > > \
op(Matrix<T>& a, const T& s) \
{ \
  typedef ConstRef<T,Matrix<T> > IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), s)); \
}

MATRIX_SCALAR_OP(operator+, DMApAdd)
MATRIX_SCALAR_OP(operator-, DMApSubs)
MATRIX_SCALAR_OP(operator*, DMApMul)
MATRIX_SCALAR_OP(operator/, DMApDiv)
#undef MATRIX_SCALAR_OP

/////////////////////////////////////////////////////////////////

#define SCALAR_MATRIX_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, ConstRef<T,Matrix<T> >, T, ap<T, ConstRef<T,Matrix<T> >, T> > > \
op(const T& s, Matrix<T>& a) \
{ \
  typedef ConstRef<T,Matrix<T> > IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), s)); \
}

SCALAR_MATRIX_OP(operator+, DMApAdd)
SCALAR_MATRIX_OP(operator-, DMApSubs)
SCALAR_MATRIX_OP(operator*, DMApMul)
SCALAR_MATRIX_OP(operator/, DMApDiv)
#undef SCALAR_MATRIX_OP

/////////////////////////////////////////////////////////////////

template<class T>
class MatrixIdentity 
{
  public:
    static inline T apply(size_t i, size_t j) 
    { return (i == j) ? T(1) : T(0); }

    static inline size_t size1() { return 0; }
    static inline size_t size2() { return 0; }
};

template<class T>
inline IterExpr<T, IdentityOp<T, MatrixIdentity<T> > > 
I()
{
  return IterExpr<T, IdentityOp<T, MatrixIdentity<T> > >(IdentityOp<T, MatrixIdentity<T> >());
}

/////////////////////////////////////////////////////////////////

template<class T>
class MatrixVectorProduct
{
  public:
    static inline T apply(ConstRef<T,Matrix<T> > m, ConstRef<T,Vector<T> > v, size_t n)
    {
      const Matrix<T> &mm = m.ref();
      const Vector<T> &vv = v.ref();
      cslice_iter<T> row = mm.row(n);
     
      T r = T(0); 
      for(size_t i = 0; i < row.size(); i++)
        r += row[i]*vv[i];

      return r; 
    }
    
    static inline size_t size1(ConstRef<T,Matrix<T> > m, ConstRef<T,Vector<T> > v) { return 1; }
    static inline size_t size2(ConstRef<T,Matrix<T> > m, ConstRef<T,Vector<T> > v) { return m.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, ConstRef<T,Matrix<T> >, ConstRef<T, Vector<T> >, MatrixVectorProduct<T> > > 
product(const Matrix<T>& m, const Vector<T>& v) 
{ 
  typedef IterBinOp<T, ConstRef<T,Matrix<T> >, ConstRef<T,Vector<T> >, MatrixVectorProduct<T> > BinExpr; 

  return IterExpr<T, BinExpr > (BinExpr(m.ref(), v.ref())); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class MatrixSliceProduct
{
  public:
    static inline T apply(ConstRef<T,Matrix<T> > m, cslice_iter<T> v, size_t n)
    {
      const Matrix<T> &mm = m.ref();
      cslice_iter<T> row = mm.row(n);

      T r = T(0); 
      for(size_t i = 0; i < row.size(); i++)
        r += row[i]*v[i];

      return r;
    }
    static inline size_t size1(ConstRef<T,Matrix<T> > m, cslice_iter<T> v) { return 1; }
    static inline size_t size2(ConstRef<T,Matrix<T> > m, cslice_iter<T> v) { return m.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, ConstRef<T,Matrix<T> >, cslice_iter<T>, MatrixSliceProduct<T> > > 
product(const Matrix<T>& m, cslice_iter<T> v) 
{ 
  typedef IterBinOp<T, ConstRef<T,Matrix<T> >, cslice_iter<T>, MatrixSliceProduct<T> > BinExpr; 

  return IterExpr<T, BinExpr > (BinExpr(m.ref(), v)); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class MatrixMatrixProduct
{
  public:
    static inline T apply(ConstRef<T,Matrix<T> > m1, ConstRef<T,Matrix<T> > m2, size_t i, size_t j)
    {
      const Matrix<T> &mm1 = m1.ref();
      const Matrix<T> &mm2 = m2.ref();
      cslice_iter<T> row = mm1.row(j);
      cslice_iter<T> column = mm2.column(i);
//      return product(row,column); 

      T r = T(0);
      for(size_t i = 0; i < std::min(row.size(), column.size()); i++)
        r += (*row++)*(*column++);

      return r;
    }
    static inline size_t size1(ConstRef<T,Matrix<T> > m1, ConstRef<T,Matrix<T> > m2) { return m2.size1(); }
    static inline size_t size2(ConstRef<T,Matrix<T> > m1, ConstRef<T,Matrix<T> > m2) { return m1.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, ConstRef<T,Matrix<T> >, ConstRef<T,Matrix<T> >, MatrixMatrixProduct<T> > > 
product(const Matrix<T>& m1, const Matrix<T>& m2) 
{ 
  typedef ConstRef<T,Matrix<T> > IT; 
  typedef IterBinOp<T, IT, IT, MatrixMatrixProduct<T> > BinExpr; 
  return IterExpr<T, BinExpr > (BinExpr(m1.ref(), m2.ref())); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class MatrixCGSliceProduct
{
  public:
    static inline T apply(ConstRef<T,Matrix<T> >_m1, cgslice_iter<T> m2, size_t i, size_t j)
    {
      const Matrix<T>&m1 = _m1.ref();
      cslice_iter<T> row = m1.row(j);
      cslice_iter<T> column = m2.column(i);
//      return product(row,column); 

      T r = T(0);
      for(size_t i = 0; i < std::min(row.size(), column.size()); i++)
        r += (*row++)*(*column++);

      return r;
    }
    static inline size_t size1(ConstRef<T,Matrix<T> >m1, cgslice_iter<T> m2) { return m2.size1(); }
    static inline size_t size2(ConstRef<T,Matrix<T> >m1, cgslice_iter<T> m2) { return m1.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, ConstRef<T,Matrix<T> >, cgslice_iter<T>, MatrixCGSliceProduct<T> > >
product(const Matrix<T> &m1, cgslice_iter<T> m2) 
{ 
  typedef ConstRef<T, Matrix<T> > IT1; 
  typedef cgslice_iter<T> IT2; 
  typedef IterBinOp<T, IT1, IT2, MatrixCGSliceProduct<T> > BinExpr; 
  return IterExpr<T, BinExpr > (BinExpr(m1.ref(),m2)); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class CGSliceMatrixProduct
{
  public:
    static inline T apply(cgslice_iter<T> m1, ConstRef<T,Matrix<T> >_m2, size_t i, size_t j)
    {
      const Matrix<T>&m2 = _m2.ref();
      cslice_iter<T> row = m1.row(j);
      cslice_iter<T> column = m2.column(i);
//      return product(row,column); 

      T r = T(0);
      for(size_t i = 0; i < std::min(row.size(), column.size()); i++)
        r += (*row++)*(*column++);

      return r;
    }
    static inline size_t size1(cgslice_iter<T> m1, ConstRef<T,Matrix<T> >m2) { return m2.size1(); }
    static inline size_t size2(cgslice_iter<T> m1, ConstRef<T,Matrix<T> >m2) { return m1.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, cgslice_iter<T>, ConstRef<T,Matrix<T> >, CGSliceMatrixProduct<T> > >
product(cgslice_iter<T> m1, const Matrix<T> &m2) 
{ 
  typedef cgslice_iter<T> IT1; 
  typedef ConstRef<T, Matrix<T> > IT2; 
  typedef IterBinOp<T, IT1, IT2, CGSliceMatrixProduct<T> > BinExpr; 
  return IterExpr<T, BinExpr > (BinExpr(m1, m2.ref())); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class CGSliceCGSliceProduct
{
  public:
    static inline T apply(cgslice_iter<T> m1, cgslice_iter<T> m2, size_t i, size_t j)
    {
      cslice_iter<T> row = m1.row(j);
      cslice_iter<T> column = m2.column(i);
//      return product(row,column); 

      T r = T(0);
      for(size_t i = 0; i < std::min(row.size(), column.size()); i++)
        r += (*row++)*(*column++);

      return r;
    }
    static inline size_t size1(cgslice_iter<T> m1, cgslice_iter<T> m2) { return m2.size1(); }
    static inline size_t size2(cgslice_iter<T> m1, cgslice_iter<T> m2) { return m1.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, cgslice_iter<T>, cgslice_iter<T>, CGSliceCGSliceProduct<T> > > 
product(cgslice_iter<T> m1, cgslice_iter<T> m2) 
{ 
  typedef cgslice_iter<T> IT; 
  typedef IterBinOp<T, IT, IT, CGSliceCGSliceProduct<T> > BinExpr; 
  return IterExpr<T, BinExpr > (BinExpr(m1,m2)); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class GSliceGSliceProduct
{
  public:
    static inline T apply(gslice_iter<T> m1, gslice_iter<T> m2, size_t i, size_t j)
    {
      slice_iter<T> row = m1.row(j);
      slice_iter<T> column = m2.column(i);

      T r = T(0);
      for(size_t i = 0; i < std::min(row.size(), column.size()); i++)
        r += (*row++)*(*column++);

      return r;
    }
    static inline size_t size1(gslice_iter<T> m1, gslice_iter<T> m2) { return m2.size1(); }
    static inline size_t size2(gslice_iter<T> m1, gslice_iter<T> m2) { return m1.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, gslice_iter<T>, gslice_iter<T>, GSliceGSliceProduct<T> > > 
product(gslice_iter<T> m1, gslice_iter<T> m2) 
{ 
  typedef gslice_iter<T> IT; 
  typedef IterBinOp<T, IT, IT, GSliceGSliceProduct<T> > BinExpr; 
  return IterExpr<T, BinExpr > (BinExpr(m1,m2)); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class MatrixMatrixTensor
{
  public:
    static inline T apply(ConstRef<T,Matrix<T> > m1, ConstRef<T,Matrix<T> > m2, size_t i, size_t j)
    {
//      const Matrix<T> &mm1 = m1.ref();
      const Matrix<T> &mm2 = m2.ref();
      size_t col1 = i / mm2.cols();
      size_t col2 = i - col1 * mm2.cols();
//      size_t col2 = i % mm2.cols();
      size_t row1 = j / mm2.rows();
      size_t row2 = j - row1 * mm2.rows();
//      size_t row2 = j % mm2.rows();
// cout << i << " " << j << " " << col1 << " " << row1 << " " <<col2 << " " <<row2<< endl;
// cout << m2(col2,row2) << endl;
// cout << m1(col1,row1) << endl;
      return m1(col1,row1)*m2(col2,row2); 
    }
    static inline size_t size1(ConstRef<T,Matrix<T> > m1, ConstRef<T,Matrix<T> > m2) { return m2.size1()*m1.size1(); }
    static inline size_t size2(ConstRef<T,Matrix<T> > m1, ConstRef<T,Matrix<T> > m2) { return m2.size2()*m1.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, ConstRef<T,Matrix<T> >, ConstRef<T,Matrix<T> >, MatrixMatrixTensor<T> > > 
tensor(const Matrix<T>& a, const Matrix<T>& b) 
{ 
  typedef ConstRef<T,Matrix<T> > IT; 
  typedef IterBinOp<T, IT, IT, MatrixMatrixTensor<T> > BinExpr; 
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b.ref())); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class VectorVectorTensor
{
  public:
    static inline T apply(ConstRef<T,Vector<T> > v1, ConstRef<T,Vector<T> > v2, size_t i, size_t j)
    {
      return v1[j]*v2[i]; 
    }
    static inline size_t size1(ConstRef<T,Vector<T> > v1, ConstRef<T,Vector<T> > v2) { return v2.size2(); }
    static inline size_t size2(ConstRef<T,Vector<T> > v1, ConstRef<T,Vector<T> > v2) { return v1.size2(); }
};

template<class T> 
IterExpr<T, IterBinOp<T, typename Vector<T>::const_iterator, typename Vector<T>::const_iterator, VectorVectorTensor<T> > > 
tensor(const Vector<T>& a, const Vector<T>& b) 
{ 
  typedef ConstRef<T,Vector<T> > IT; 
  typedef IterBinOp<T, IT, IT, VectorVectorTensor<T> > BinExpr; 
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b.ref())); 
}

/////////////////////////////////////////////////////////////////

template<class T>
class MatrixTranspose
{
  public:
    static inline T apply(ConstRef<T,Matrix<T> > m, size_t i, size_t j)
    {
      return m(j, i); 
    }
    static inline size_t size1(ConstRef<T,Matrix<T> > m) { return m.size2(); }
    static inline size_t size2(ConstRef<T,Matrix<T> > m) { return m.size1(); }
};

template<class T>
class MatrixCTranspose
{
  public:
    static inline T apply(ConstRef<T,Matrix<T> > m, size_t i, size_t j)
    {
      return std::conj(m(j, i)); 
    }
    static inline size_t size1(ConstRef<T,Matrix<T> > m) { return m.size2(); }
    static inline size_t size2(ConstRef<T,Matrix<T> > m) { return m.size1(); }
};

template<class T>
class MatrixConj
{
  public:
    static inline T apply(ConstRef<T,Matrix<T> > m, size_t i, size_t j)
    {
      return std::conj(m(i, j)); 
    }
    static inline size_t size1(ConstRef<T,Matrix<T> > m) { return m.size1(); }
    static inline size_t size2(ConstRef<T,Matrix<T> > m) { return m.size2(); }
};

#define MATRIX_METHOD(op, ap) \
template<class T> \
IterExpr<T, IterOp<T, ConstRef<T,Matrix<T> >, ap<T> > > \
op(const Matrix<T> &a) \
{ \
  typedef ConstRef<T,Matrix<T> > IT; \
  typedef IterOp<T, IT, ap<T> > Expr; \
  return IterExpr<T, Expr > (Expr(a.ref())); \
}

MATRIX_METHOD(ctranspose, MatrixCTranspose);
MATRIX_METHOD(transpose, MatrixTranspose);
MATRIX_METHOD(conj, MatrixConj);
#undef MATRIX_METHOD

/////////////////////////////////////////////////////////////////

template<class T>
inline Matrix<T>
inverse(Matrix<T>& a)
{
  if(a.num_rows != a.num_cols)
    cerr << "** Warning: Matrix dimensions are different\n";

  Matrix<T> aux(a);
  int l, j, i, m = std::min(a.rows(), a.num());
  T r;

  for(i = 0; i < m; i++){
    r = (T)1.0f / aux[i][i];
    aux[i][i] = (T)0.0f;
    for(j = 0; j < m; j++){
      aux[j][i] *= r;
      for(l = 0; l < m; l++) aux[j][l] = aux[j][l] - aux[i][l]*aux[j][i];
    }
    for(j = 0; j < m; j++) aux[i][j] *= (-r);
    aux[i][i] = r;
  }

  return aux;
}

/////////////////////////////////////////////////////////////////

#define MATRIX_GSLICE_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, ConstRef<T, Matrix<T> >, gslice_iter<T>, ap<T, ConstRef<T, Matrix<T> >, gslice_iter<T> > > > \
op(Matrix<T>& a, gslice_iter<T>& b) \
{ \
  typedef ConstRef<T, Matrix<T> > IT1; \
  typedef gslice_iter<T> IT2; \
  typedef IterBinOp<T, IT1, IT2, ap<T, IT1, IT2> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b.begin())); \
}

MATRIX_GSLICE_OP(operator+, DMApAdd)
MATRIX_GSLICE_OP(operator-, DMApSubs)
MATRIX_GSLICE_OP(operator*, DMApMul)
MATRIX_GSLICE_OP(operator/, DMApDiv)
#undef MATRIX_GSLICE_OP

/////////////////////////////////////////////////////////////////

#define GSLICE_MATRIX_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, ConstRef<T, Matrix<T> >, gslice_iter<T>, ap<T, ConstRef<T, Matrix<T> >, gslice_iter<T> > > > \
op(gslice_iter<T>& a, Matrix<T>& b) \
{ \
  typedef ConstRef<T, Matrix<T> > IT1; \
  typedef gslice_iter<T> IT2; \
  typedef IterBinOp<T, IT1, IT2, ap<T, IT1, IT2> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(b, a.begin())); \
}

GSLICE_MATRIX_OP(operator+, DMApAdd)
GSLICE_MATRIX_OP(operator-, DMApSubs)
GSLICE_MATRIX_OP(operator*, DMApMul)
GSLICE_MATRIX_OP(operator/, DMApDiv)
#undef GSLICE_MATRIX_OP

/////////////////////////////////////////////////////////////////

#endif /* __DMTK_MATRIX_IMPLEMENT_H__ */
