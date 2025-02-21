/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef __DMTK_META_IMPLEMENT_H__
#define __DMTK_META_IMPLEMENT_H__

// Iterator-Iterator binary operations
// This class contains the references to the objects that
// participate in the operation.

namespace dmtk
{

template<class T, class Iter1, class Iter2, class Op>
class IterBinOp
{
  private:

    Iter1 _i1; 
    Iter2 _i2;

  public:

    IterBinOp(const Iter1 i1, const Iter2 i2): _i1(i1), _i2(i2){}

    T operator[](size_t n) const { return Op::apply(_i1, _i2, n); } 

    T operator()(size_t n) const { return Op::apply(_i1, _i2, n); } 

    T operator()(size_t i, size_t j) const { return Op::apply(_i1, _i2, i, j); }

    T operator*() { return Op::apply(*_i1, *_i2); } 

    IterBinOp<T,Iter1,Iter2,Op>& operator++()
    {
      ++_i1; ++_i2;
      return *this;
    }
    size_t size1() const { return Op::size1(_i1,_i2); }
    size_t size2() const { return Op::size2(_i1,_i2); }
};

// Scalar-Iterator operations
// This class contains the references to the objects that
// participate in the operation.

template<class T, class Iter, class Scalar, class Op>
class IterScalarOp
{
  private:

    Iter _i; 
    Scalar _s;

  public:

    IterScalarOp(Iter i, const Scalar s): _i(i), _s(s){}

    T operator[](size_t n) const { return Op::apply(_i[n], _s); }

    T operator()(size_t i, size_t j) const { return Op::apply(_i(i,j), _s); }

    T operator*() const { return Op::apply(*_i, _s); }

    IterScalarOp<T,Iter,Scalar,Op>& operator++()
    {
      ++_i; 
      return *this;
    }
    size_t size1() const { return _i.size1(); }
    size_t size2() const { return _i.size2(); }
};

// Unary operations on one single Iterator
// This class contains the references to the Iterator subexpression

template<class T, class Iter, class Op>
class IterOp
{
  private:

    Iter _i; 

  public:

    IterOp(Iter i): _i(i){}

    T operator[](size_t n) const { return Op::apply(_i, n); }

    T operator()(size_t n) const { return Op::apply(_i, n); }

    T operator()(size_t i, size_t j) const { return Op::apply(_i, i, j); }

    T operator*() const { return Op::apply(*_i); }

    IterOp<T,Iter,Op>& operator++()
    {
      ++_i; 
      return *this;
    }
    size_t size1() const { return Op::size1(_i); }
    size_t size2() const { return Op::size2(_i); }
};

// Virtual operator 
// It returns the value of Op

template<class T, class Op>
class IdentityOp
{
  public:
    T operator[](size_t n) const { return Op::apply(n); }

    T operator()(size_t i, size_t j) const { return Op::apply(i, j); } 

    size_t size1() const { return Op::size1(); }
    size_t size2() const { return Op::size2(); }
};

// Container for a general subexpression

template<class T, class Child>
class ConstRef
{
  private:
    const Child& _i;

  public:
    ConstRef(const Child &a) : _i(a){}

    T operator[](size_t n) const { return _i.operator[](n); } 
    T operator()(size_t i, size_t j) const { return _i.operator()(i, j); }   

    const Child& ref() const { return _i; }

    size_t size1() const { return _i.size1(); }
    size_t size2() const { return _i.size2(); }
};

// Expression-wrapper class, used as container interface for an object.
// It can be IterOp, IdentityOp, IterBinOp, or IterScalarOp.

template<class T, class Expr>
class IterExpr
{
  private:
    Expr iter;

  public:
    IterExpr(const Expr &a) : iter(a){}

    T operator*() const { return *iter; }  
    T operator[](size_t n) const { return iter[n]; } 
    T operator()(size_t n) const { return iter(n); }  
    T operator()(size_t i, size_t j) const { return iter(i, j); }  

    IterExpr<T,Expr>& operator++() { ++iter; return *this; } 
    size_t size1() const { return iter.size1(); }
    size_t size2() const { return iter.size2(); }
};

template<class T>
class DMApAdd0
{
  public:
    static inline T apply(T a, T b) { return a + b; }
};
 
template<class T>
class DMApSubs0
{
  public:
    static inline T apply(T a, T b) { return a - b; }
};
 
template<class T, class A, class B>
class DMApAdd
{
  public:
    static inline T apply(T a, T b) { return a + b; }
    static inline T apply(A a, B b, size_t n) { return a[n] + b[n]; }
    static inline T apply(A a, B b, size_t i, size_t j) { return a(i, j) + b(i, j); }
    static inline size_t size1(T a, T b) { return 1; }
    static inline size_t size2(T a, T b) { return 1; }
    static inline size_t size1(A a, B b) { return std::min(a.size1(),b.size1()); }
    static inline size_t size2(A a, B b) { return std::min(a.size2(),b.size2()); }
};

template<class T, class A, class B>
class DMApSubs
{
  public:
    static inline T apply(T a, T b) { return a - b; }
    static inline T apply(A a, B b, size_t n) { return a[n] - b[n]; }
    static inline T apply(A a, B b, size_t i, size_t j) { return a(i, j) - b(i, j); }
    static inline size_t size1(T a, T b) { return 1; }
    static inline size_t size2(T a, T b) { return 1; }
    static inline size_t size1(A a, B b) { return std::min(a.size1(),b.size1()); }
    static inline size_t size2(A a, B b) { return std::min(a.size2(),b.size2()); }
};

template<class T, class A, class B>
class DMApMul
{
  public:
    static inline T apply(T a, T b) { return a * b; }
    static inline T apply(A a, B b, size_t n) { return a[n] * b[n]; }
    static inline T apply(A a, B b, size_t i, size_t j) { return a(i, j) * b(i, j); }
    static inline size_t size1(T a, T b) { return 1; }
    static inline size_t size2(T a, T b) { return 1; }
    static inline size_t size1(A a, B b) { return std::min(a.size1(),b.size1()); }
    static inline size_t size2(A a, B b) { return std::min(a.size2(),b.size2()); }
};

template<class T, class A, class B>
class DMApDiv
{
  public:
    static inline T apply(T a, T b) { return a / b; }
    static inline T apply(A a, B b, size_t n) { return a[n] / b[n]; }
    static inline T apply(A a, B b, size_t i, size_t j) { return a(i, j) / b(i, j); }
    static inline size_t size1(T a, T b) { return 1; }
    static inline size_t size2(T a, T b) { return 1; }
    static inline size_t size1(A a, B b) { return std::min(a.size1(),b.size1()); }
    static inline size_t size2(A a, B b) { return std::min(a.size2(),b.size2()); }
};

template<class T, class S>
class DMApIdentity
{
  public:
    static inline T apply(S m, size_t i) { return m[i]; }
    static inline T apply(S m, size_t i, size_t j) { return m(i, j); }
    static inline size_t size1(T a) { return 1; }
    static inline size_t size2(T a) { return 1; }
    static inline size_t size1(S a) { return a.size1(); }
    static inline size_t size2(S a) { return a.size2(); }
};

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

#define SCALAR_EXPR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterScalarOp<T, IterExpr<T, Expr>, T, ap<T, IterExpr<T, Expr>, T> > > \
op(const T& s, IterExpr<T, Expr> a) \
{ \
  typedef IterScalarOp<T, IterExpr<T, Expr>, T, ap<T, IterExpr<T, Expr>, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, s)); \
}

SCALAR_EXPR_OP(operator+, DMApAdd)
SCALAR_EXPR_OP(operator-, DMApSubs)
SCALAR_EXPR_OP(operator*, DMApMul)
SCALAR_EXPR_OP(operator/, DMApDiv)
#undef SCALAR_EXPR_OP

/////////////////////////////////////////////////////////////////

#define EXPR_SCALAR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterScalarOp<T, IterExpr<T, Expr>, T, ap<T, IterExpr<T, Expr>, T> > > \
op(IterExpr<T, Expr> a, const T&s) \
{ \
  typedef IterScalarOp<T, IterExpr<T, Expr>, T, ap<T, IterExpr<T, Expr>, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, s)); \
}

EXPR_SCALAR_OP(operator+, DMApAdd)
EXPR_SCALAR_OP(operator-, DMApSubs)
EXPR_SCALAR_OP(operator*, DMApMul)
EXPR_SCALAR_OP(operator/, DMApDiv)
#undef EXPR_SCALAR_OP

/////////////////////////////////////////////////////////////////

#define EXPR_EXPR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, IterExpr<T, Expr>, IterExpr<T, Expr>, ap<T, IterExpr<T, Expr>, IterExpr<T, Expr> > > > \
op(IterExpr<T, Expr> a, IterExpr<T, Expr> b) \
{ \
  typedef IterBinOp<T, IterExpr<T, Expr>, IterExpr<T, Expr>, ap<T, IterExpr<T, Expr>, IterExpr<T, Expr> > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b)); \
}

EXPR_EXPR_OP(operator+, DMApAdd)
EXPR_EXPR_OP(operator-, DMApSubs)
EXPR_EXPR_OP(operator*, DMApMul)
EXPR_EXPR_OP(operator/, DMApDiv)
#undef EXPR_EXPR_OP

#define EXPR_EXPR_OP(op, ap) \
template<class T, class Expr1, class Expr2> \
IterExpr<T, IterBinOp<T, IterExpr<T, Expr1>, IterExpr<T, Expr2>, ap<T, IterExpr<T, Expr1>, IterExpr<T, Expr2> > > > \
op(IterExpr<T, Expr1> a, IterExpr<T, Expr2> b) \
{ \
  typedef IterBinOp<T, IterExpr<T, Expr1>, IterExpr<T, Expr2>, ap<T, IterExpr<T, Expr1>, IterExpr<T, Expr2> > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b)); \
}

EXPR_EXPR_OP(operator+, DMApAdd)
EXPR_EXPR_OP(operator-, DMApSubs)
EXPR_EXPR_OP(operator*, DMApMul)
EXPR_EXPR_OP(operator/, DMApDiv)
#undef EXPR_EXPR_OP

/////////////////////////////////////////////////////////////////

} //namespace dmtk

#endif /* __DMTK_META_IMPLEMENT_H__ */
