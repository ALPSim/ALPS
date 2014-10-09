/*********************************************************
*  Emanuel Gull, 2006 - 2007
*  BLAS Compatibility headers
*  blas rssvector class.
*********************************************************/
#ifndef BLAS_RSSVECTOR
#define BLAS_RSSVECTOR

#include "./blasheader.h"
#include <cmath>
#include <iostream>

namespace blas{
  //forward declaration of matrix:
  class smatrix;
  
  //simple BLAS rssvector that uses BLAS calls for operations.
  class rssvector{
public:
    rssvector(int size){
      values_=size>0?new float[size]:0;
      size_=size;
      memory_size_=size;
    }
    rssvector(){
      size_=0;
      memory_size_=0;
      values_=0;
    }
    rssvector(int size, float initial_value){
      if(size<0) {std::cerr<<"choose positive size!"<<std::endl; abort();}
      values_=size>0?new float[size]:0;
      size_=size;
      memory_size_=size;
      for(int i=0;i<size_;++i){
        values_[i]=initial_value;
      }
    }
    ~rssvector(){
      if(memory_size_>0 && values_!=0){
        delete[] values_;
      }
    }
    rssvector& operator=(const rssvector &V){
      if(&V==this) return *this;
      if(memory_size_< V.size_){
        size_=V.size_;
        memory_size_=V.size_;
        delete[] values_;
        values_=size_>0?new float[size_]:0;
      }
      if(size_>0)
        memcpy(values_, V.values_, size_*sizeof(float));
      return *this;
    }
    rssvector(const rssvector &V){
      values_=V.size_>0?new float[V.size_]:0;
      size_=V.size_;
      memory_size_=V.size_;
      if(size_>0)
        memcpy(values_, V.values_, sizeof(float)*V.size_);
    }
    inline float &operator()(const int i){return *(values_+i);}
    inline const float &operator()(const int i) const {return *(values_+i);}
    inline int size()const{return size_;}
    //inner product
    inline float operator*(const rssvector &v2){
      int inc=1;
      return sdot_(&size_, values_,&inc,v2.values_,&inc);
    }
    inline float dot(const rssvector &v2){ return (*this)*v2; }
    //multiply rssvector by constant
    inline rssvector operator *=(float lambda){
      int inc=1;
      sscal_(&size_, &lambda, values_, &inc);
      return *this;
    }
    void clear(){
      memset(values_, 0, sizeof(float)*size_);
    }
    void exp( float c){
      //std::cout<<"exp-ing for c: "<<c<<std::endl;
#ifdef VECLIB
      //osx veclib exp
      for(int i=0;i<size_;++i){ values_[i]*=c; }
      vecLib::vvexp_(values_, values_, &size_); 
#else
#ifdef ACML
      //amd acml rssvector exp
      float scaled_values[size_];
      for(int i=0;i<size_;++i){ scaled_values[i]=values_[i]*c; }
      acml::vrda_exp(size_, scaled_values, values_);
#else
#ifdef MKL
      //intel MKL rssvector exp
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
    void exp(const float &c, const blas::rssvector &v){
      resize(v.size());
      for(int i=0;i<size_;++i){
        operator()(i)=std::exp(c*v(i));
      }
    }
    void resize(int new_size){
      if(new_size==size_) return;
      if(new_size<memory_size_){
        size_=new_size;
        return;
      }
      if(size_!=0){
        delete[] values_;
      }
      size_=new_size;
      memory_size_=new_size;
      values_=size_>0?new float[size_]:0;
    }
private:
    int size_; //current size of matrix
    int memory_size_; //current size of matrix
    float *values_; //where the actual values are stored
    friend class blas::smatrix;
  };
  inline std::ostream &operator<<(std::ostream &os, const rssvector&M){
    os<<"[ ";
    for(int i=0;i<M.size();++i){
      os<<M(i)<<" ";
    }
    os<<"]"<<std::endl;
    return os;
  }
} //namespace

#endif 
