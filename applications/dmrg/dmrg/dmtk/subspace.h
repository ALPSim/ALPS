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

#ifndef __DMTK_SUBSPACE_H__
#define __DMTK_SUBSPACE_H__

#include <iostream>
#include <iosfwd>
#include "qn.h" 
#include "range.h" 

namespace dmtk
{

class SubSpace: public Range
{
  private:
    QN _qn;
  public:

    SubSpace(): Range(0,0), _qn(0) {};
    SubSpace(const QN &qn, int _begin, int _end): 
      Range(_begin,_end), _qn(qn)  {}
    SubSpace(const Range &s): Range(s) {}

    size_t dim() const { return std::slice::size(); }
    QN& qn() { return _qn; }
    QN qn() const { return _qn; }

    // Streams

    void write(std::ostream &s) const
    {
      _qn.write(s);
      int _begin = begin(), _end = end();
      s.write((const char *)&_begin, sizeof(int));
      s.write((const char *)&_end, sizeof(int));
    }

    void read(std::istream &s)
    {
      QN qn;
      int _begin, _end;

      qn.read(s);
      s.read((char *)&_begin, sizeof(int));
      s.read((char *)&_end, sizeof(int));

      *this = SubSpace(qn, _begin, _end);
    }

};

} // namespace dmtk

#endif // __DMTK_SUBSPACE_H__
