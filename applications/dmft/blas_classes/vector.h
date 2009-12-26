/*****************************************************************************
 *
 * ALPS DMFT Project - BLAS Compatibility headers
 *  Square Matrix Class
 *
 * Copyright (C) 2005 - 2009 by 
 *                              Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/
 
#ifndef BLAS_VECTOR
#define BLAS_VECTOR

#include "./blasheader.h"
#include <string.h> //memcpy
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace blas{
  //forward declaration of matrix:
  class matrix;
  
  //simple BLAS vector that uses BLAS calls for operations.
  class vector{
public:
    vector(int size){
      values_=size>0?new double[size]:0;
      size_=size;
    }
    vector(){
      size_=0;
      values_=0;
    }
    vector(int size, int initial_value){
      if(size<0) throw(std::invalid_argument("choose positive size!"));
      values_=size>0?new double[size]:0;
      size_=size;
      for(int i=0;i<size_;++i){
        values_[i]=initial_value;
      }
    }
    ~vector(){
      if(size_>0 && values_!=0){
        //std::cout<<"deleting values for size: "<<size_<<" values: "<<values_<<std::endl;
        delete[] values_;
      }
    }
    vector& operator=(const vector &V){
      if(&V==this) return *this;
      if(size_!= V.size_){
        size_=V.size_;
        delete[] values_;
        values_=size_>0?new double[size_]:0;
      }
      if(size_>0)
        memcpy(values_, V.values_, size_*sizeof(double));
      return *this;
    }
    vector(const vector &V){
      values_=V.size_>0?new double[V.size_]:0;
      size_=V.size_;
      if(size_>0)
        memcpy(values_, V.values_, sizeof(double)*V.size_);
    }
    inline double &operator()(const int i){return *(values_+i);}
    inline const double &operator()(const int i) const {return *(values_+i);}
    inline int size()const{return size_;}
    //inner product
    inline double operator*(const vector &v2){
      int inc=1;
      return ddot_(&size_, values_,&inc,v2.values_,&inc);
    }
    inline double dot(const vector &v2){ return (*this)*v2; }
    //multiply vector by constant
    inline vector operator *=(double lambda){
      int inc=1;
      dscal_(&size_, &lambda, values_, &inc);
      return *this;
    }
    void clear(){
      memset(values_, 0, sizeof(double)*size_);
    }
    inline const double *values() const{return values_; }
    inline double *values() {return values_; }
    void exp( double c){
      //std::cout<<"exp-ing for c: "<<c<<std::endl;
#ifdef VECLIB
      //osx veclib exp
      //for(int i=0;i<size_;++i){ values_[i]*=c; }
      int one=1;
      blas::dscal_(&size_,&c,values_,&one);
      vecLib::vvexp(values_, values_, &size_); 
#else
#ifdef ACML
      //amd acml vector exp
      double scaled_values[size_];
      daxpy_(&size_, &c, values_, &inc, scaled_values, &inc);
      acml::vrda_exp(size_, scaled_values, values_);
#else
#ifdef MKL
      //intel MKL vector exp
      double scaled_values[size_];
      memset(scaled_values, 0,sizeof(double)*size_);
      int inc=1;
      daxpy_(&size_, &c, values_, &inc, scaled_values, &inc);
      mkl::vdExp(size_, scaled_values, values_);
#else
      //pedestrian method
      for(int i=0;i<size_;++i){
        values_[i]=std::exp(c*values_[i]);
      }
#endif
#endif
#endif
    }
    void exp(const double &c, const blas::vector &v){
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
      values_=size_>0?new double[size_]:0;
    }
private:
    int size_; //current size of matrix
    double *values_; //where the actual values are stored
    friend class blas::matrix;
  };
  inline std::ostream &operator<<(std::ostream &os, const vector&M){
    os<<"[ ";
    for(int i=0;i<M.size();++i){
      os<<M(i)<<" ";
    }
    os<<"]"<<std::endl;
    return os;
  }
} //namespace

#endif 
