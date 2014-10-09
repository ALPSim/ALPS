/*********************************************************
*  Emanuel Gull, 2006 - 2007
*  BLAS Compatibility headers
*  blas vector class.
*********************************************************/
#ifndef BLAS_SVECTOR
#define BLAS_SVECTOR

#include "./blasheader.h"
#include <cmath>
#include <iostream>

namespace blas{
  //forward declaration of matrix:
  class smatrix;
  
  //simple BLAS vector that uses BLAS calls for operations.
  class svector{
public:
    svector(int size){
      values_=size>0?new float[size]:0;
      size_=size;
    }
    svector(){
      size_=0;
      values_=0;
    }
    svector(int size, int initial_value){
      if(size<0) std::cerr<<"choose positive size!"<<std::endl; abort();
      values_=size>0?new float[size]:0;
      size_=size;
      for(int i=0;i<size_;++i){
        values_[i]=initial_value;
      }
    }
    ~svector(){
      if(size_>0 && values_!=0){
        //std::cout<<"deleting values for size: "<<size_<<" values: "<<values_<<std::endl;
        delete[] values_;
      }
    }
    svector& operator=(const svector &V){
      if(&V==this) return *this;
      if(size_!= V.size_){
        size_=V.size_;
        delete[] values_;
        values_=size_>0?new float[size_]:0;
      }
      if(size_>0)
        memcpy(values_, V.values_, size_*sizeof(float));
      return *this;
    }
    svector(const svector &V){
      values_=V.size_>0?new float[V.size_]:0;
      size_=V.size_;
      if(size_>0)
        memcpy(values_, V.values_, sizeof(float)*V.size_);
    }
    inline float &operator()(const int i){return *(values_+i);}
    inline const float &operator()(const int i) const {return *(values_+i);}
    inline int size()const{return size_;}
    //inner product
    inline float operator*(const svector &v2){
      int inc=1;
      return sdot_(&size_, values_,&inc,v2.values_,&inc);
    }
    inline float dot(const svector &v2){ return (*this)*v2; }
    //multiply vector by constant
    inline svector operator *=(float lambda){
      int inc=1;
      sscal_(&size_, &lambda, values_, &inc);
      return *this;
    }
    void clear(){
      memset(values_, 0, sizeof(float)*size_);
    }
    inline const float *values() const{return values_; }
    inline float *values() {return values_; }
    void exp( float c){
      //std::cout<<"exp-ing for c: "<<c<<std::endl;
#ifdef VECLIB
      //osx veclib exp
      for(int i=0;i<size_;++i){ values_[i]*=c; }
      vecLib::vvexp_(values_, values_, &size_); 
#else
#ifdef ACML
      //amd acml vector exp
      float scaled_values[size_];
      for(int i=0;i<size_;++i){ scaled_values[i]=values_[i]*c; }
      acml::vrda_exp(size_, scaled_values, values_);
#else
#ifdef MKL
      //intel MKL vector exp
      float scaled_values[size_];
      memset(scaled_values, 0,sizeof(float)*size_);
      int inc=1;
      saxpy_(&size_, &c, values_, &inc, scaled_values, &inc);
      mkl::vfExp(size_, scaled_values, values_);
#else
      //pedestrian method
      for(int i=0;i<size_;++i){
        values_[i]=std::exp(c*values_[i]);
      }
#endif
#endif
#endif
    }
    void exp(const float &c, const blas::svector &v){
      resize(v.size());
      for(int i=0;i<size_;++i){
        operator()(i)=std::exp(c*v(i));
      }
    }
    void resize(int new_size){
      if(new_size==size_) return;
      if(size_!=0){
        //std::cout<<"deleting values"<<std::endl;
        delete[] values_;
        //std::cout<<"done"<<std::endl;
      }
      size_=new_size;
      values_=size_>0?new float[size_]:0;
    }
private:
    int size_; //current size of matrix
    float *values_; //where the actual values are stored
    friend class blas::smatrix;
  };
  inline std::ostream &operator<<(std::ostream &os, const svector&M){
    os<<"[ ";
    for(int i=0;i<M.size();++i){
      os<<M(i)<<" ";
    }
    os<<"]"<<std::endl;
    return os;
  }
} //namespace

#endif 
