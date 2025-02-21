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

#ifndef __RANGE_H__
#define __RANGE_H__

namespace dmtk
{

class Range: public std::slice
{
  public:
    Range(size_t _o, size_t _d, size_t _s): std::slice(_o, (_d-_o+1)/_s, _s) {}
    Range(size_t _o, size_t _d): std::slice(_o, (_d-_o+1), 1) {}
    Range(const Range &r): std::slice(r) {}

    bool operator==(const Range &r)
     { 
       if(begin() == r.begin() && end() == r.end() && stride() == r.stride())
         return true;
       return false;
     }

    size_t begin() const { return std::slice::start(); }
    size_t end() const { return std::slice::start()+std::slice::size()-1; }
};

} // namespace dmtk

#endif // __RANGE_H__
