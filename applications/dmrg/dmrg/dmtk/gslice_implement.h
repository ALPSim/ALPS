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

#ifndef __DMTK_GSLICE_IMPLEMENT_H__
#define __DMTK_GSLICE_IMPLEMENT_H__

/////////////////////////////////////////////////////////////////

#define GSLICE_GSLICE_OP(op, ap) \
template<class T> \
IterExpr<T, IterBinOp<T, cgslice_iter<T>, cgslice_iter<T>, ap<T, cgslice_iter<T>, cgslice_iter<T> > > > \
op(cgslice_iter<T>& a, cgslice_iter<T>& b) \
{ \
  typedef cgslice_iter<T> IT; \
  typedef IterBinOp<T, IT, IT, ap<T, IT, IT> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), b.begin())); \
}

GSLICE_GSLICE_OP(operator+, DMApAdd)
GSLICE_GSLICE_OP(operator-, DMApSubs)
GSLICE_GSLICE_OP(operator*, DMApMul)
GSLICE_GSLICE_OP(operator/, DMApDiv)
#undef GSLICE_GSLICE_OP

/////////////////////////////////////////////////////////////////

#define EXPR_GSLICE_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, IterExpr<T, Expr>, cgslice_iter<T>, ap<T, IterExpr<T, Expr>, cgslice_iter<T> > > > \
op(IterExpr<T, Expr> a, cgslice_iter<T>& b) \
{ \
  typedef cgslice_iter<T> IT; \
  typedef IterBinOp<T, IterExpr<T, Expr>, IT, ap<T, IterExpr<T, Expr>, IT > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a, b.begin())); \
}

EXPR_GSLICE_OP(operator+, DMApAdd)
EXPR_GSLICE_OP(operator-, DMApSubs)
EXPR_GSLICE_OP(operator*, DMApMul)
EXPR_GSLICE_OP(operator/, DMApDiv)
#undef EXPR_GSLICE_OP

/////////////////////////////////////////////////////////////////

#define GSLICE_EXPR_OP(op, ap) \
template<class T, class Expr> \
IterExpr<T, IterBinOp<T, cgslice_iter<T>, IterExpr<T, Expr> , ap<T, cgslice_iter<T>, IterExpr<T, Expr> > > > \
op(cgslice_iter<T>& a, IterExpr<T, Expr> b) \
{ \
  typedef cgslice_iter<T> IT; \
  typedef IterBinOp<T, IT, IterExpr<T, Expr>, ap<T, IT, IterExpr<T, Expr> > > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), b)); \
}

GSLICE_EXPR_OP(operator+, DMApAdd)
GSLICE_EXPR_OP(operator-, DMApSubs)
GSLICE_EXPR_OP(operator*, DMApMul)
GSLICE_EXPR_OP(operator/, DMApDiv)
#undef GSLICE_EXPR_OP

/////////////////////////////////////////////////////////////////

#define GSLICE_SCALAR_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, cgslice_iter<T>, T, ap<T, cgslice_iter<T>, T> > > \
op(cgslice_iter<T>& a, const T& s) \
{ \
  typedef cgslice_iter<T> IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), s)); \
}

GSLICE_SCALAR_OP(operator+, DMApAdd)
GSLICE_SCALAR_OP(operator-, DMApSubs)
GSLICE_SCALAR_OP(operator*, DMApMul)
GSLICE_SCALAR_OP(operator/, DMApDiv)
#undef GSLICE_SCALAR_OP

/////////////////////////////////////////////////////////////////

#define SCALAR_GSLICE_OP(op, ap) \
template<class T> \
IterExpr<T, IterScalarOp<T, cgslice_iter<T>, T, ap<T, cgslice_iter<T>, T> > > \
op(const T& s, cgslice_iter<T>& a) \
{ \
  typedef cgslice_iter<T> IT; \
  typedef IterScalarOp<T, IT, T, ap<T, IT, T> > BinExpr; \
  return IterExpr<T, BinExpr > (BinExpr(a.begin(), s)); \
}

SCALAR_GSLICE_OP(operator+, DMApAdd)
SCALAR_GSLICE_OP(operator-, DMApSubs)
SCALAR_GSLICE_OP(operator*, DMApMul)
SCALAR_GSLICE_OP(operator/, DMApDiv)
#undef SCALAR_GSLICE_OP

/////////////////////////////////////////////////////////////////

#endif /* __DMTK_GSLICE_IMPLEMENT_H__ */
