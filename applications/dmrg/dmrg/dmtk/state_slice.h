#ifndef __DMTK_STATE_SLICE_H__
#define __DMTK_STATE_SLICE_H__

#include <valarray>
#include "slice_iter.h" 
#include "gslice_iter.h" 

namespace dmtk
{

template <class T>
class state_slice
{
  public:
    state_slice(): v(0),s1(slice(0,0,0)),s2(slice(0,0,0)),s3(slice(0,0,0)),s4(slice(0,0,0)) {};
    state_slice(valarray<T>* vv, slice _s1, slice _s2, slice _s3, slice _s4):
      v(vv),s1(_s1),s2(_s2),s3(_s3),s4(_s4){};
    state_slice(const state_slice<T>& ss):
      v(ss.v),s1(ss.s1),s2(ss.s2),s3(ss.s3),s4(ss.s4){}

    state_slice& operator=(state_slice<T> ss)
      { 
        size_t n1 = std::min(ss.s1.size(), s1.size());
        size_t n2 = std::min(ss.s2.size(), s2.size());
        size_t n3 = std::min(ss.s3.size(), s3.size());
        size_t n4 = std::min(ss.s4.size(), s4.size());
        for(size_t i = 0; i < n1; i++)
          for(size_t j = 0; j < n2; j++)
            for(size_t k = 0; k < n3; k++)
              for(size_t l = 0; l < n4; l++)
                this->operator()(i,j,k,l) = ss(i,j,k,l);
        return *this; 
      }

    gslice_iter<T> operator()(size_t i1, size_t i2, slice r3, slice r4) const
      {
        slice _s3(s1.start()+i1*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+r3.start()*s3.stride(),
                  r3.size(), r3.stride()*s3.stride());
        slice _s4(s4.start()+r4.start()*s4.stride(),r4.size(),
                  r4.stride()*s4.stride());
        return gslice_iter<T>(v,_s3,_s4);
      }
    gslice_iter<T> operator()(size_t i1, slice r2, slice r3, size_t i4) const
      {
        slice _s2(s1.start()+i1*s1.stride()+
                  s2.start()+r2.start()*s2.stride()+
                  s4.start()+i4*s4.stride(),
                  r2.size(), r2.stride()*s2.stride());
        slice _s3(s3.start()+r3.start()*s3.stride(),r3.size(),
                  r3.stride()*s3.stride());
        return gslice_iter<T>(v,_s2,_s3);
      }
    gslice_iter<T> operator()(slice r1, slice r2, size_t i3, size_t i4) const
      {
        slice _s1(s1.start()+r1.start()*s1.stride()+
                  s3.start()+i3*s3.stride()+
                  s4.start()+i4*s4.stride(),
                  r1.size(), r1.stride()*s1.stride());
        slice _s2(s2.start()+r2.start()*s2.stride(),r2.size(),
                  r2.stride()*s2.stride());
        return gslice_iter<T>(v,_s1,_s2);
      }
    gslice_iter<T> operator()(size_t i1, slice r2, size_t i3, slice r4) const
      {
        slice _s2(s1.start()+i1*s1.stride()+
                  s2.start()+r2.start()*s2.stride()+
                  s3.start()+i3*s3.stride(),
                  r2.size(), r2.stride()*s2.stride());
        slice _s4(s4.start()+r4.start()*s4.stride(),r4.size(),
                  r4.stride()*s4.stride());
        return gslice_iter<T>(v,_s2,_s4);
      }
    gslice_iter<T> operator()(slice r1, size_t i2, size_t i3, slice r4) const
      {
        slice _s1(s1.start()+r1.start()*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+i3*s3.stride(),
                  r1.size(), r1.stride()*s1.stride());
        slice _s4(s4.start()+r4.start()*s4.stride(),r4.size(),
                  r4.stride()*s4.stride());
        return gslice_iter<T>(v,_s1,_s4);
      }
    gslice_iter<T> operator()(slice r1, size_t i2, slice r3, size_t i4) const
      {
        slice _s1(s1.start()+r1.start()*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s4.start()+i4*s4.stride(),
                  r1.size(), r1.stride()*s1.stride());
        slice _s3(s3.start()+r3.start()*s3.stride(),r3.size(),
                  r3.stride()*s3.stride());
        return gslice_iter<T>(v,_s1,_s3);
      }
    slice_iter<T> operator()(size_t i1, size_t i2, size_t i3, slice r4) const
      {
        slice _s4(s1.start()+i1*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+i3*s3.stride()+
                  s4.start()+r4.start()*s4.stride(),
                  r4.size(), r4.stride()*s4.stride());
        return slice_iter<T>(v,_s4);
      }
    slice_iter<T> operator()(size_t i1, size_t i2, slice r3, size_t i4) const
      {
        slice _s3(s1.start()+i1*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+r3.start()*s3.stride()+
                  s4.start()+i4*s4.stride(),
                  r3.size(), r3.stride()*s3.stride());
        return slice_iter<T>(v,_s3);
      }
    slice_iter<T> operator()(size_t i1, slice r2, size_t i3, size_t i4) const
      {
        slice _s2(s1.start()+i1*s1.stride()+
                  s2.start()+r2.start()*s3.stride()+
                  s3.start()+i3*s3.stride()+
                  s4.start()+i4*s4.stride(),
                  r2.size(), r2.stride()*s2.stride());
        return slice_iter<T>(v,_s2);
      }
    slice_iter<T> operator()(slice r1, size_t i2, size_t i3, size_t i4) const
      {
        slice _s1(s1.start()+r1.start()*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+i3*s3.stride()+
                  s4.start()+i4*s4.stride(),
                  r1.size(), r1.stride()*s1.stride());
        return slice_iter<T>(v,_s1);
      }

    size_t stride1() const { return s1.stride(); }
    size_t stride2() const { return s2.stride(); }
    size_t stride3() const { return s3.stride(); }
    size_t stride4() const { return s4.stride(); }

    size_t start1() const { return s1.start(); }
    size_t start2() const { return s2.start(); }
    size_t start3() const { return s3.start(); }
    size_t start4() const { return s4.start(); }

    size_t size1() const { return s1.size(); }
    size_t size2() const { return s2.size(); }
    size_t size3() const { return s3.size(); }
    size_t size4() const { return s4.size(); }

    size_t size() const { return size1() * size2() * size3() * size4(); }

    T& operator()(size_t i, size_t j, size_t k, size_t l)
      {return ref(i, j, k, l);}
    T& operator()(size_t i)
      {return (*v)(i);}

    size_t index(size_t i, size_t j, size_t k, size_t l) const
      { 
        size_t _index = s1.start() + i * s1.stride() +
                        s2.start() + j * s2.stride() +
                        s3.start() + k * s3.stride() +
                        s4.start() + l * s4.stride();
        return _index;
      }

  private:
    valarray<T>* v;
    slice s1;
    slice s2;
    slice s3;
    slice s4;

    T& ref(size_t i, size_t j, size_t k , size_t l) 
      { 
        size_t _index = s1.start() + i * s1.stride() +
                        s2.start() + j * s2.stride() +
                        s3.start() + k * s3.stride() +
                        s4.start() + l * s4.stride();
        return (*v)[_index]; 
      }
};

template <class T>
class cstate_slice
{
  public:
    cstate_slice(): v(0),s1(slice(0,0,0)),s2(slice(0,0,0)),s3(slice(0,0,0)),s4(slice(0,0,0)) {};
    cstate_slice(const valarray<T>* vv, slice _s1, slice _s2, slice _s3, slice _s4):
      v(vv),s1(_s1),s2(_s2),s3(_s3),s4(_s4){};
    cstate_slice(const cstate_slice<T>& ss):
      v(ss.v),s1(ss.s1),s2(ss.s2),s3(ss.s3),s4(ss.s4){}

    cgslice_iter<T> operator()(size_t i1, size_t i2, slice r3, slice r4) const
      {
        slice _s3(s1.start()+i1*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+r3.start()*s3.stride(),
                  r3.size(), r3.stride()*s3.stride());
        slice _s4(s4.start()+r4.start()*s4.stride(),r4.size(),
                  r4.stride()*s4.stride());
        return cgslice_iter<T>(v,_s3,_s4);
      }
    cgslice_iter<T> operator()(size_t i1, slice r2, slice r3, size_t i4) const
      {
        slice _s2(s1.start()+i1*s1.stride()+
                  s2.start()+r2.start()*s2.stride()+
                  s4.start()+i4*s4.stride(),
                  r2.size(), r2.stride()*s2.stride());
        slice _s3(s3.start()+r3.start()*s3.stride(),r3.size(),
                  r3.stride()*s3.stride());
        return cgslice_iter<T>(v,_s2,_s3);
      }
    cgslice_iter<T> operator()(slice r1, slice r2, size_t i3, size_t i4) const
      {
        slice _s1(s1.start()+r1.start()*s1.stride()+
                  s3.start()+i3*s3.stride()+
                  s4.start()+i4*s4.stride(),
                  r1.size(), r1.stride()*s1.stride());
        slice _s2(s2.start()+r2.start()*s2.stride(),r2.size(),
                  r2.stride()*s2.stride());
        return cgslice_iter<T>(v,_s1,_s2);
      }
    cgslice_iter<T> operator()(size_t i1, slice r2, size_t i3, slice r4) const
      {
        slice _s2(s1.start()+i1*s1.stride()+
                  s2.start()+r2.start()*s2.stride()+
                  s3.start()+i3*s3.stride(),
                  r2.size(), r2.stride()*s2.stride());
        slice _s4(s4.start()+r4.start()*s4.stride(),r4.size(),
                  r4.stride()*s4.stride());
        return cgslice_iter<T>(v,_s2,_s4);
      }
    cgslice_iter<T> operator()(slice r1, size_t i2, size_t i3, slice r4) const
      {
        slice _s1(s1.start()+r1.start()*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+i3*s3.stride(),
                  r1.size(), r1.stride()*s1.stride());
        slice _s4(s4.start()+r4.start()*s4.stride(),r4.size(),
                  r4.stride()*s4.stride());
        return cgslice_iter<T>(v,_s1,_s4);
      }
    cgslice_iter<T> operator()(slice r1, size_t i2, slice r3, size_t i4) const
      {
        slice _s1(s1.start()+r1.start()*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s4.start()+i4*s4.stride(),
                  r1.size(), r1.stride()*s1.stride());
        slice _s3(s3.start()+r3.start()*s3.stride(),r3.size(),
                  r3.stride()*s3.stride());
        return cgslice_iter<T>(v,_s1,_s3);
      }
    cslice_iter<T> operator()(size_t i1, size_t i2, size_t i3, slice r4) const
      {
        slice _s4(s1.start()+i1*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+i3*s3.stride()+
                  s4.start()+r4.start()*s4.stride(),
                  r4.size(), r4.stride()*s4.stride());
        return cslice_iter<T>(v,_s4);
      }
    cslice_iter<T> operator()(size_t i1, size_t i2, slice r3, size_t i4) const
      {
        slice _s3(s1.start()+i1*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+r3.start()*s3.stride()+
                  s4.start()+i4*s4.stride(),
                  r3.size(), r3.stride()*s3.stride());
        return cslice_iter<T>(v,_s3);
      }
    cslice_iter<T> operator()(size_t i1, slice r2, size_t i3, size_t i4) const
      {
        slice _s2(s1.start()+i1*s1.stride()+
                  s2.start()+r2.start()*s2.stride()+
                  s3.start()+i3*s3.stride()+
                  s4.start()+i4*s4.stride(),
                  r2.size(), r2.stride()*s2.stride());
        return cslice_iter<T>(v,_s2);
      }
    cslice_iter<T> operator()(slice r1, size_t i2, size_t i3, size_t i4) const
      {
        slice _s1(s1.start()+r1.start()*s1.stride()+
                  s2.start()+i2*s2.stride()+
                  s3.start()+i3*s3.stride()+
                  s4.start()+i4*s4.stride(),
                  r1.size(), r1.stride()*s1.stride());
        return cslice_iter<T>(v,_s1);
      }

    size_t start1() const { return s1.start(); }
    size_t start2() const { return s2.start(); }
    size_t start3() const { return s3.start(); }
    size_t start4() const { return s4.start(); }

    size_t stride1() const { return s1.stride(); }
    size_t stride2() const { return s2.stride(); }
    size_t stride3() const { return s3.stride(); }
    size_t stride4() const { return s4.stride(); }

    size_t size1() const { return s1.size(); }
    size_t size2() const { return s2.size(); }
    size_t size3() const { return s3.size(); }
    size_t size4() const { return s4.size(); }

    size_t size() const { return size1() + size2() + size3() + size4(); }

    T operator()(size_t i, size_t j, size_t k, size_t l) const
      {return ref(i, j, k, l);}
    T operator()(size_t i) const
      {return (*v)(i);}

    size_t index(size_t i, size_t j, size_t k, size_t l) const
      { 
        size_t _index = s1.start() + i * s1.stride() +
                        s2.start() + j * s2.stride() +
                        s3.start() + k * s3.stride() +
                        s4.start() + l * s4.stride();
        return _index;
      }

  private:
    const valarray<T>* v;
    slice s1;
    slice s2;
    slice s3;
    slice s4;

    T ref(size_t i, size_t j, size_t k, size_t l) const
      { 
        size_t _index = s1.start() + i * s1.stride() +
                        s2.start() + j * s2.stride() +
                        s3.start() + k * s3.stride() +
                        s4.start() + l * s4.stride();
        return (*v)[_index]; 
      }
};

} // namespace dmtk 

#endif // __DMTK_STATE_SLICE_H__


