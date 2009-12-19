#ifndef __DMTK_VECTOR_IMPLEMENT_H__
#define __DMTK_VECTOR_IMPLEMENT_H__


template<class T>
class VectorConj
{
  public:
    static inline T apply(ConstRef<T,Vector<T> > v, size_t i)
    {
      return std::conj(v(i));
    }
    static inline size_t size1(ConstRef<T, Vector<T> > v) { return v.size1(); }
    static inline size_t size2(ConstRef<T, Vector<T> > v) { return v.size2(); }
};

#define VECTOR_METHOD(op, ap) \
template<class T> \
IterExpr<T, IterOp<T, ConstRef<T,Vector<T> >, ap<T> > > \
op(const Vector<T> &a) \
{ \
  typedef ConstRef<T,Vector<T> > IT; \
  typedef IterOp<T, IT, ap<T> > Expr; \
  return IterExpr<T, Expr > (Expr(ConstRef<T,Vector<T> >(a))); \
}

VECTOR_METHOD(conj, VectorConj);
#undef VECTOR_METHOD

/////////////////////////////////////////////////////////////////

#define VECTOR_VECTOR_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, ConstRef<T,Vector<T> >, ConstRef<T,Vector<T> >, ap<T, ConstRef<T,Vector<T> >, ConstRef<T,Vector<T> > > > >\
op(const Vector<T>& a, const Vector<T>& b) \
{ \
  typedef ConstRef<T,Vector<T> > IT; \
  typedef IterBinOp<T, IT, IT, ap<T, IT, IT> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b.ref())); \
}

VECTOR_VECTOR_OP(operator+, DMApAdd)
VECTOR_VECTOR_OP(operator-, DMApSubs)
VECTOR_VECTOR_OP(operator*, DMApMul)
VECTOR_VECTOR_OP(operator/, DMApDiv)
#undef VECTOR_VECTOR_OP

/////////////////////////////////////////////////////////////////

#define EXPR_VECTOR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, IterExpr<T, Expr>, ConstRef<T,Vector<T> >, ap<T, IterExpr<T, Expr>, ConstRef<T,Vector<T> > > > > \
op(IterExpr<T, Expr> a, Vector<T>& b) \
{ \
  typedef ConstRef<T,Vector<T> > IT; \
  typedef IterBinOp<T, IterExpr<T, Expr>, IT, ap<T, IterExpr<T, Expr>, IT> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b.ref())); \
}

EXPR_VECTOR_OP(operator+, DMApAdd)
EXPR_VECTOR_OP(operator-, DMApSubs)
EXPR_VECTOR_OP(operator*, DMApMul)
EXPR_VECTOR_OP(operator/, DMApDiv)
#undef EXPR_VECTOR_OP

/////////////////////////////////////////////////////////////////

#define VECTOR_EXPR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, ConstRef<T,Vector<T> >, IterExpr<T, Expr> , ap<T, ConstRef<T,Vector<T> >, IterExpr<T, Expr> > > > \
op(Vector<T>& a, IterExpr<T, Expr> b) \
{ \
  typedef ConstRef<T,Vector<T> > IT; \
  typedef IterBinOp<T, IT, IterExpr<T, Expr>, ap<T, IT, IterExpr<T, Expr> > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), b)); \
}

VECTOR_EXPR_OP(operator+, DMApAdd)
VECTOR_EXPR_OP(operator-, DMApSubs)
VECTOR_EXPR_OP(operator*, DMApMul)
VECTOR_EXPR_OP(operator/, DMApDiv)
#undef VECTOR_EXPR_OP

/////////////////////////////////////////////////////////////////

#define VECTOR_SCALAR_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, ConstRef<T,Vector<T> >, T, ap<T, ConstRef<T,Vector<T> >, T> > > \
op(Vector<T>& a, const T& s) \
{ \
  typedef ConstRef<T,Vector<T> > IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), s)); \
}

VECTOR_SCALAR_OP(operator+, DMApAdd)
VECTOR_SCALAR_OP(operator-, DMApSubs)
VECTOR_SCALAR_OP(operator*, DMApMul)
VECTOR_SCALAR_OP(operator/, DMApDiv)
#undef VECTOR_SCALAR_OP

/////////////////////////////////////////////////////////////////

#define SCALAR_VECTOR_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, ConstRef<T,Vector<T> >, T, ap<T, ConstRef<T,Vector<T> >, T> > > \
op(const T& s, Vector<T>& a) \
{ \
  typedef ConstRef<T,Vector<T> > IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.ref(), s)); \
}

SCALAR_VECTOR_OP(operator+, DMApAdd)
SCALAR_VECTOR_OP(operator-, DMApSubs)
SCALAR_VECTOR_OP(operator*, DMApMul)
SCALAR_VECTOR_OP(operator/, DMApDiv)
#undef SCALAR_VECTOR_OP

/////////////////////////////////////////////////////////////////

#define SLICE_VECTOR_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, cslice_iter<T>, typename Vector<T>::const_iterator, ap<T, cslice_iter<T>, typename Vector<T>::const_iterator > > > \
op(cslice_iter<T>& a, typename Vector<T>::const_iterator& b) \
{ \
  typedef cslice_iter<T> IT1; \
  typedef typename Vector<T>::const_iterator IT2; \
  typedef IterBinOp<T, IT1, IT2, ap<T, IT1, IT2> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), b)); \
}

SLICE_VECTOR_OP(operator+, DMApAdd)
SLICE_VECTOR_OP(operator-, DMApSubs)
SLICE_VECTOR_OP(operator*, DMApMul)
SLICE_VECTOR_OP(operator/, DMApDiv)
#undef SLICE_VECTOR_OP

/////////////////////////////////////////////////////////////////

#define VECTOR_SLICE_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, typename Vector<T>::const_iterator, cslice_iter<T>, ap<T, cslice_iter<T>, typename Vector<T>::const_iterator > > > \
op(typename Vector<T>::const_iterator& a, cslice_iter<T>& b) \
{ \
  typedef typename Vector<T>::const_iterator IT1; \
  typedef cslice_iter<T> IT2; \
  typedef IterBinOp<T, IT1, IT2, ap<T, IT1, IT2> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b.begin())); \
}

VECTOR_SLICE_OP(operator+, DMApAdd)
VECTOR_SLICE_OP(operator-, DMApSubs)
VECTOR_SLICE_OP(operator*, DMApMul)
VECTOR_SLICE_OP(operator/, DMApDiv)
#undef VECTOR_SLICE_OP

/////////////////////////////////////////////////////////////////

template<class T>
inline T
product(cslice_iter<T> a, typename Vector<T>::const_iterator b)
{
  T r = T(0);

  a = a.begin();
  typename Vector<T>::const_iterator _b = b;

  for(size_t i = 0; i < a.size(); i++)
    r += std::conj(*a++) * (*_b++);

  return r;
}

template<class T>
inline T
product(typename Vector<T>::const_iterator a, cslice_iter<T> b)
{
  T r = T(0);

  b = b.begin();
  typename Vector<T>::const_iterator _a = a;

  for(size_t i = 0; i < a.size(); i++)
    r += std::conj(*_a++) * (*b++);

  return r;
}

/////////////////////////////////////////////////////////////////

template<class T>
class GSliceVectorProduct
{
  public:
    static inline T apply(cgslice_iter<T> m, typename Vector<T>::const_iterator v, size_t n)
    {
      return product(m[n],v);
    }
    static inline size_t size1(cgslice_iter<T> m, typename Vector<T>::const_iterator v) { return 1; }
    static inline size_t size2(cgslice_iter<T> m, typename Vector<T>::const_iterator v) { return m.size2(); }
};

template<class T>
IterExpr<T, IterBinOp<T, cgslice_iter<T>, typename Vector<T>::const_iterator, GSliceVectorProduct<T> > >
product(cslice_iter<T> m, const Vector<T>& v)
{
  typedef IterBinOp<T, cgslice_iter<T>, typename Vector<T>::const_iterator, GSliceVectorProduct<T> > BinExpr;

  return IterExpr<T, BinExpr > (BinExpr(m, v.begin()));
}

/////////////////////////////////////////////////////////////////

template<class T>
gslice_iter<T>
Vector<T>::as_matrix(size_t ncols, size_t nrows)
{
  std::slice s_col(0, ncols, nrows);
  std::slice s_row(0, nrows, 1);
                                                                                
  gslice_iter<T> iter(this, s_col, s_row);
  return iter;
}
                                                                                
template<class T>
cgslice_iter<T>
Vector<T>::as_matrix(size_t ncols, size_t nrows) const
{
  std::slice s_col(0, ncols, nrows);
  std::slice s_row(0, nrows, 1);
                                                                                
  cgslice_iter<T> iter(this, s_col, s_row);
  return iter;
}

/////////////////////////////////////////////////////////////////

#endif /* __DMTK_VECTOR_IMPLEMENT_H__ */
