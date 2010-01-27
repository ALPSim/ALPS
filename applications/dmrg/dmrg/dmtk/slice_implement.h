/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef __DMTK_SLICE_IMPLEMENT_H__
#define __DMTK_SLICE_IMPLEMENT_H__

/////////////////////////////////////////////////////////////////

#define SLICE_SLICE_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, cslice_iter<T>, cslice_iter<T>, ap<T, cslice_iter<T>, cslice_iter<T> > > > \
op(cslice_iter<T> a, cslice_iter<T> b) \
{ \
  typedef cslice_iter<T> IT; \
  typedef IterBinOp<T, IT, IT, ap<T, IT, IT> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), b.begin())); \
}

SLICE_SLICE_OP(operator+, DMApAdd)
SLICE_SLICE_OP(operator-, DMApSubs)
SLICE_SLICE_OP(operator*, DMApMul)
SLICE_SLICE_OP(operator/, DMApDiv)
#undef SLICE_SLICE_OP

/////////////////////////////////////////////////////////////////

#define EXPR_SLICE_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, IterExpr<T, Expr>, cslice_iter<T>, ap<T, IterExpr<T, Expr>, cslice_iter<T> > > > \
op(IterExpr<T, Expr> a, cslice_iter<T> b) \
{ \
  typedef cslice_iter<T> IT; \
  typedef IterBinOp<T, IterExpr<T, Expr>, IT, ap<T, IterExpr<T, Expr>, IT > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b.begin())); \
}

EXPR_SLICE_OP(operator+, DMApAdd)
EXPR_SLICE_OP(operator-, DMApSubs)
EXPR_SLICE_OP(operator*, DMApMul)
EXPR_SLICE_OP(operator/, DMApDiv)
#undef EXPR_SLICE_OP

/////////////////////////////////////////////////////////////////

#define SLICE_EXPR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, cslice_iter<T>, IterExpr<T, Expr> , ap<T, cslice_iter<T>, IterExpr<T, Expr> > > > \
op(cslice_iter<T> a, IterExpr<T, Expr> b) \
{ \
  typedef cslice_iter<T> IT; \
  typedef IterBinOp<T, IT, IterExpr<T, Expr>, ap<T, IT, IterExpr<T, Expr> > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), b)); \
}

SLICE_EXPR_OP(operator+, DMApAdd)
SLICE_EXPR_OP(operator-, DMApSubs)
SLICE_EXPR_OP(operator*, DMApMul)
SLICE_EXPR_OP(operator/, DMApDiv)
#undef SLICE_EXPR_OP

/////////////////////////////////////////////////////////////////

#define SLICE_SCALAR_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, cslice_iter<T>, T, ap<T, cslice_iter<T>, T> > > \
op(cslice_iter<T> a, const T& s) \
{ \
  typedef cslice_iter<T> IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), s)); \
}

SLICE_SCALAR_OP(operator+, DMApAdd)
SLICE_SCALAR_OP(operator-, DMApSubs)
SLICE_SCALAR_OP(operator*, DMApMul)
SLICE_SCALAR_OP(operator/, DMApDiv)
#undef SLICE_SCALAR_OP

/////////////////////////////////////////////////////////////////

#define SCALAR_SLICE_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, cslice_iter<T>, T, ap<T, cslice_iter<T>, T> > > \
op(const T& s, cslice_iter<T> a) \
{ \
  typedef cslice_iter<T> IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), s)); \
}

SCALAR_SLICE_OP(operator+, DMApAdd)
SCALAR_SLICE_OP(operator-, DMApSubs)
SCALAR_SLICE_OP(operator*, DMApMul)
SCALAR_SLICE_OP(operator/, DMApDiv)
#undef SCALAR_SLICE_OP

/////////////////////////////////////////////////////////////////

#define SLICE_SCALAR_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, slice_iter<T>, T, ap<T, slice_iter<T>, T> > > \
op(slice_iter<T> a, const T& s) \
{ \
  typedef slice_iter<T> IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), s)); \
}

SLICE_SCALAR_OP(operator+, DMApAdd)
SLICE_SCALAR_OP(operator-, DMApSubs)
SLICE_SCALAR_OP(operator*, DMApMul)
SLICE_SCALAR_OP(operator/, DMApDiv)
#undef SLICE_SCALAR_OP

/////////////////////////////////////////////////////////////////

#define SCALAR_SLICE_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, slice_iter<T>, T, ap<T, slice_iter<T>, T> > > \
op(const T& s, slice_iter<T> a) \
{ \
  typedef slice_iter<T> IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), s)); \
}

SCALAR_SLICE_OP(operator+, DMApAdd)
SCALAR_SLICE_OP(operator-, DMApSubs)
SCALAR_SLICE_OP(operator*, DMApMul)
SCALAR_SLICE_OP(operator/, DMApDiv)
#undef SCALAR_SLICE_OP

/////////////////////////////////////////////////////////////////

#endif /* __DMTK_SLICE_IMPLEMENT_H__ */
