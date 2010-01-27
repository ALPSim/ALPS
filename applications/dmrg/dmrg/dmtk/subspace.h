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
