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
