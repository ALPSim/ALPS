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
