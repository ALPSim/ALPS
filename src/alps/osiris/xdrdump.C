/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/osiris/xdrdump.h>
#include <alps/osiris/std/string.h>

#include <boost/throw_exception.hpp>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace alps {

static_assert(sizeof(long long) == 8, "XDR requires 64-bit long long");

namespace detail {

bool xdr_bool(XDR *xdrs, bool *bp)
{
  if (xdrs->x_op == XDR_ENCODE) {
    bool_t b = *bp;
    return ::xdr_bool(xdrs, &b);
  } else if (xdrs->x_op == XDR_DECODE) {
    bool_t b;
    bool retval = ::xdr_bool(xdrs, &b);
    *bp=b;
    return retval;
  } else if (xdrs->x_op == XDR_FREE) {
    return true;
  }
  return false;
}

static bool xdr_s_char(XDR *xdrs, signed char *scp) 
{
  if (xdrs->x_op == XDR_ENCODE) {
    char c = *scp;
    return xdr_char(xdrs, &c);
  } else if (xdrs->x_op == XDR_DECODE) {
    char c;
    bool retval = ::xdr_char(xdrs, &c);
    *scp = c;
    return retval;
  } else if (xdrs->x_op == XDR_FREE) {
    return true;
  }
  return false;
}

bool xdr_u_hyper(XDR *xdrs, unsigned long long *value)
{
  uint32_t high = 0, low = 0;
  if (xdrs->x_op == XDR_ENCODE) {
    high = static_cast<uint32_t>(*value >> 32);
    low = static_cast<uint32_t>(*value);
  }
  if (!::xdr_uint32_t(xdrs, &high) || !::xdr_uint32_t(xdrs, &low)) return false;
  if (xdrs->x_op == XDR_DECODE)
    *value = (static_cast<unsigned long long>(high) << 32) | low;
  return true;
}

bool xdr_hyper(XDR *xdrs, long long *value)
{
  unsigned long long bits = 0;
  if (xdrs->x_op == XDR_ENCODE) std::memcpy(&bits, value, sizeof(bits));
  if (!alps::detail::xdr_u_hyper(xdrs, &bits)) return false;
  if (xdrs->x_op == XDR_DECODE) std::memcpy(value, &bits, sizeof(bits));
  return true;
}

bool xdr_long_8(XDR *xdrs, long *lp)
{
  long long t;
  if (xdrs->x_op == XDR_ENCODE) {
    t = (long long)(*lp);
    return alps::detail::xdr_hyper(xdrs, &t);
  } else if (xdrs->x_op == XDR_DECODE) {
    if (!alps::detail::xdr_hyper(xdrs, &t)) return false;
    *lp = (long)t;
    return true;
  } else if (xdrs->x_op == XDR_FREE) {
    return true;
  }
  return false;
}

bool xdr_u_long_8(XDR *xdrs, unsigned long *lp)
{
  unsigned long long t;
  if (xdrs->x_op == XDR_ENCODE) {
    t = (unsigned long long)(*lp);
    return alps::detail::xdr_u_hyper(xdrs, &t);
  } else if (xdrs->x_op == XDR_DECODE) {
    if (!alps::detail::xdr_u_hyper(xdrs, &t)) return false;
    *lp = (unsigned long)t;
    return true;
  } else if (xdrs->x_op == XDR_FREE) {
    return true;
  }
  return false;
}

bool xdr_long_double(XDR *xdrs, long double *ldp) 
{
  if (xdrs->x_op == XDR_ENCODE) {
    double high = *ldp;
    double low  = (*ldp-high);
    return xdr_double(xdrs, &high) && xdr_double(xdrs, &low);
  } else if (xdrs->x_op == XDR_DECODE) {
    double high = 0.;
    double low  = 0.;
    bool retval = xdr_double(xdrs, &high) && xdr_double(xdrs, &low); 
    *ldp = low + high;
    return retval;
  } else if (xdrs->x_op == XDR_FREE) {
    return true;
  }
  return false;
}

template<class T, int N>
struct xdr_helper {};

#define ALPS_DUMP_DO_TYPE(T,X) \
  template<int N> struct xdr_helper<T, N> { \
    static bool xdr_do_type(XDR * xdrs, T * v) { return X (xdrs, v); } \
  };
#define ALPS_DUMP_DO_TYPE_N(T,N,X) \
  template<> struct xdr_helper<T, N> { \
    static bool xdr_do_type(XDR * xdrs, T * v) { return X (xdrs, v); } \
  };
ALPS_DUMP_DO_TYPE(bool, alps::detail::xdr_bool)
ALPS_DUMP_DO_TYPE(char, xdr_char)
ALPS_DUMP_DO_TYPE(signed char, xdr_s_char)
ALPS_DUMP_DO_TYPE(unsigned char, xdr_u_char)
ALPS_DUMP_DO_TYPE(short, xdr_short)
ALPS_DUMP_DO_TYPE(unsigned short, xdr_u_short)
ALPS_DUMP_DO_TYPE(int, xdr_int)
ALPS_DUMP_DO_TYPE(unsigned int, xdr_u_int)
ALPS_DUMP_DO_TYPE_N(long, 4, xdr_long)
ALPS_DUMP_DO_TYPE_N(unsigned long, 4, xdr_u_long)

ALPS_DUMP_DO_TYPE_N(long, 8, xdr_long_8)
ALPS_DUMP_DO_TYPE_N(unsigned long, 8, xdr_u_long_8)
#ifdef BOOST_HAS_LONG_LONG
ALPS_DUMP_DO_TYPE(long long, alps::detail::xdr_hyper)
ALPS_DUMP_DO_TYPE(unsigned long long, alps::detail::xdr_u_hyper)
#endif
ALPS_DUMP_DO_TYPE(float, xdr_float)
ALPS_DUMP_DO_TYPE(double, xdr_double)
ALPS_DUMP_DO_TYPE(long double, xdr_long_double)
#undef ALPS_DUMP_DO_TYPE
#undef ALPS_DUMP_DO_TYPE_N

} // namespace detail

//-----------------------------------------------------------------------
// get and set the position in the stream
//-----------------------------------------------------------------------

uint32_t OXDRDump::getPosition() const
{
  return xdr_getpos((XDR*) &xdr_); // cast to non-const necessary
}

void OXDRDump::setPosition(uint32_t pos)
{
  if (!xdr_setpos(&xdr_,pos))
    boost::throw_exception(std::runtime_error("failed to reposition OXDRDump"));
}

#define ALPS_DUMP_DO_TYPE(T) \
void OXDRDump::write_simple(T x)  \
{ \
  if (!detail::xdr_helper<T, int(sizeof(T))>::xdr_do_type(&xdr_, const_cast<T*>(&x))) \
    boost::throw_exception(std::runtime_error("failed to write type "#T" to an OXDRDump"));\
} \
void OXDRDump::write_array(size_t n, const T* p)  \
{ \
  int l = n; \
  if (!xdr_vector(&xdr_, reinterpret_cast<char*>(const_cast<T*>(p)), l, int(sizeof(T)), (xdrproc_t) &detail::xdr_helper<T, int(sizeof(T))>::xdr_do_type)) \
    boost::throw_exception ( std::runtime_error("failed to write array of type "#T" to an OXDRDump")); \
} \
void IXDRDump::read_simple(T& x)\
{ \
  if (!detail::xdr_helper<T, int(sizeof(T))>::xdr_do_type(&xdr_, &x)) \
    boost::throw_exception(std::runtime_error("failed to read type "#T" from an IXDRDump")); \
} \
void IXDRDump::read_array(size_t n, T* p) \
{ \
  int l = n; \
  if (!xdr_vector(&xdr_, reinterpret_cast<char*>(p), l, int(sizeof(T)), (xdrproc_t) &detail::xdr_helper<T, int(sizeof(T))>::xdr_do_type)) \
    boost::throw_exception ( std::runtime_error("failed to read array of type "#T" from an IXDRDump")); \
}

ALPS_DUMP_DO_TYPE(bool)
ALPS_DUMP_DO_TYPE(char)
ALPS_DUMP_DO_TYPE(signed char)
ALPS_DUMP_DO_TYPE(unsigned char)
ALPS_DUMP_DO_TYPE(short)
ALPS_DUMP_DO_TYPE(unsigned short)
ALPS_DUMP_DO_TYPE(int)
ALPS_DUMP_DO_TYPE(unsigned int)
ALPS_DUMP_DO_TYPE(long)
ALPS_DUMP_DO_TYPE(unsigned long)
#ifdef BOOST_HAS_LONG_LONG
ALPS_DUMP_DO_TYPE(long long)
ALPS_DUMP_DO_TYPE(unsigned long long)
#endif
ALPS_DUMP_DO_TYPE(float)
ALPS_DUMP_DO_TYPE(double)
ALPS_DUMP_DO_TYPE(long double)

#undef ALPS_DUMP_DO_TYPE

void OXDRDump::write_string(size_t n,const char *p) 
{ 
  int l=n; 
  char* ncp = const_cast<char*>(p);
  if (!xdr_string(&xdr_,&ncp,l))
    boost::throw_exception ( std::runtime_error("failed to write a string to an OXDRDump"));
} 

void IXDRDump::read_string(size_t n, char *p) 
{ 
  int l=n; 
  if (!xdr_string(&xdr_,&p,l))
    boost::throw_exception ( std::runtime_error("failed to read a string from an IXDRDump"));
} 

//-----------------------------------------------------------------------
// get and set the position in the stream
//-----------------------------------------------------------------------

uint32_t IXDRDump::getPosition() const 
{
  return xdr_getpos((XDR*) &xdr_);
}

void IXDRDump::setPosition(uint32_t pos)
{
  if (!xdr_setpos(&xdr_,pos))
    boost::throw_exception( std::runtime_error("failed to reposition IXDRDump"));
}

//=======================================================================
// OXDRFileDump
// 
// implements a dump for writing into a file using the XDR format
//-----------------------------------------------------------------------

// reopen a file
void OXDRFileDump::open_file(const std::string& fn,bool append)
{
  file_ = std::fopen(fn.c_str(),(append ? "ab" : "wb"));
  if(file_)
      xdrstdio_create(&xdr_,file_,XDR_ENCODE);
  else  {
      // opening failed
      std::string text = "failed to open file \"";
      text += fn;
      text += "\" for writing";
      boost::throw_exception(std::runtime_error(text));
    }
}


// create a new dump file
OXDRFileDump::OXDRFileDump(const boost::filesystem::path& fn, bool append)
{
  open_file(fn.string(),append);
}

// destructor closes the stream and file
OXDRFileDump::~OXDRFileDump()
{
  xdr_destroy(&xdr_);
  if(file_)
    std::fclose(file_);
}

void OXDRFileDump::flush() 
{
  std::fflush(file_);
}


//=======================================================================
// IXDRFileDump
// 
// implements a dump for reading from a file using the XDR format
//-----------------------------------------------------------------------

// open a dump file
IXDRFileDump::IXDRFileDump(const boost::filesystem::path& p)
{
  open_file(p.string());
}

// open a file for reading at a specified position
void IXDRFileDump::open_file(const std::string& fn)
{
  valid_ = true;
  file_ = std::fopen(fn.c_str(),"rb");

  if(file_) // open succeeded
    xdrstdio_create(&xdr_,file_,XDR_DECODE);
  else   {
      // open failed
      std::string text = "failed to open file ";
      text += fn;
      text += " for reading";
      valid_=false;
#ifndef BOOST_NO_EXCEPTIONS
      boost::throw_exception (std::runtime_error(text));
#else
      std::cerr << "Osiris error: " << text << "\n";
#endif
    }
}


// destructor closes XDR stream and the file
IXDRFileDump::~IXDRFileDump()
{
  if(valid_) {
      xdr_destroy(&xdr_);
      if(file_)
        std::fclose(file_);
    }
}

} // namespace alps
