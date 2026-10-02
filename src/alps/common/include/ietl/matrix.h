/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: matrix.h,v 1.3 2003/03/17 10:52:47 prakash Exp $ */

namespace ietl {
  template <class T>
    class FortranMatrix {
    public:
    typedef std::size_t size_type;
    FortranMatrix(size_type n, size_type m):n_(n),m_(m){
      p = new T[m*n];
    }
    
    ~FortranMatrix() {delete[] p; }
    T* data() { return p;}
    const T* data() const { return p;}

    T operator()(size_type i, size_type j) const {
      return p[i + j*n_];
    }
    
    T& operator()(size_type i, size_type j) {
      return p[i + j*n_];
    }
    void resize(size_type n, size_type m) {
      m_ = m;
      n_ = n;
      delete[] p;
      p = new T[m_*n_];
    }
    
  size_type minor() {
    return n_;
  } 
   
    private:
  FortranMatrix(const FortranMatrix<T>&) {}
  void operator=(const FortranMatrix<T> &){}
  T* p;
  unsigned int n_;
  unsigned int m_;  
  };
} // end of namespace.
