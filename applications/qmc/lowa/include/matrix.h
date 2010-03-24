/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Lode Pollet <lpollet@physics.harvard.edu>,
*                            Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id: matrix.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- done
 * 2) Replacing raw pointers -- not done yet
 *
 */



#ifndef ALPS_APPLICATIONS_LOWA_MATRIX_H
#define ALPS_APPLICATIONS_LOWA_MATRIX_H

#include <iostream>
#include <cmath>
#include <cassert>


namespace alps {
namespace applications {
namespace lowa {


template <class S, class T>
class matrix
{
public:
  typedef S size_type;
  typedef S index_type;
  typedef T value_type;  


  inline size_type size() const {return _size;}
  inline size_type dim1() const {return _row;}
  inline size_type dim2() const {return _column;}

  inline value_type& operator() (index_type const row, index_type const column) const 
  { return _values[row*_column + column];}    //CHANGED TO FORTRAN STYLE!!!


private:
  bool _is_initialized;

  value_type* _values;
  size_type   _row;
  size_type   _column;
  size_type   _size;


public:
  // destructor
  ~matrix()  
  {  
    delete [] _values;  
  }


  // constructors
  matrix() 
    : _is_initialized(false)
  {}

  matrix(size_type row, size_type column)      
    : _is_initialized(true) 

    , _row(row)
    , _column(column)
    , _size(row*column)
  {
    _values = new value_type [_size];
  }

  matrix(matrix const & mymatrix) 
    : _is_initialized(mymatrix._is_initialized)

    , _row(mymatrix._row)
    , _column(mymatrix._column)
    , _size(mymatrix._size)

    , _values(new value_type [mymatrix._size])
  {
    std::copy(&mymatrix._values[0], &mymatrix._values[0]+_size, &_values[0]);
  }


  // assignment operator
  matrix& operator=(const matrix& mymatrix) {
    if (!_is_initialized)
    {
      _is_initialized(true);

      _row = mymatrix._row;
      _column = mymatrix._column;
      _size  = mymatrix._size;

      _values = new value_type [_size];
      std::copy(&mymatrix._values[0], &mymatrix._values[0]+_size, &_values[0]);

      return *this;
    }
    else
    {
      assert(_size == mymatrix._size);
      std::copy(&mymatrix._values[0], &mymatrix._values[0]+_size, &_values[0]);

      return *this;
    }
  }

  void init(size_type row, size_type column) 
  {  
    _is_initialized = true;

    _row = row; 
    _column = column;
    _size  = row*column; 

    _values = new value_type [_size];
  }

};


} // ending namespace lowa
} // ending namespace applications
} // ending namespace alps

#endif
