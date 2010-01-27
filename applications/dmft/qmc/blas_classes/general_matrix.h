/*****************************************************************************
 *
 * ALPS DMFT Project - BLAS Compatibility headers
 *  General Matrix Class with BLAS bindings.
 *
 * Copyright (C) 2005 - 2009 by 
 *                              Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

#ifndef GENBLAS_MATRIX
#define GENBLAS_MATRIX

#include "./blasheader.h"
#include "./vector.h"
#include "./resizeable_vector.h"
#include <cmath>
#include <cassert>
#include <complex>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sys/types.h>
#ifdef USE_MATRIX_DISPATCH //use dispatch to tiny matrix functions for small matrices
#undef __APPLE_CC__

#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#define MSIZE_P (1)(2)(4) //use tiny matrix dispatcher for these sizes.
#define MSIZE_Q (1)(2)(4)
#define MSIZE_R (1)(2)(4)

#define MLOOP(unused, product) \
BOOST_PP_CAT(\
BOOST_PP_CAT(\
BOOST_PP_CAT(\
BOOST_PP_CAT(\
BOOST_PP_CAT(\
void matmult_, BOOST_PP_SEQ_ELEM(0, product)\
)\
, _\
)\
, BOOST_PP_SEQ_ELEM(1, product)\
)\
, _\
)\
, BOOST_PP_SEQ_ELEM(2, product)\
)(const double * A, const double * B, double * C);
//that produces the function definitions:
BOOST_PP_SEQ_FOR_EACH_PRODUCT(MLOOP, (MSIZE_P)(MSIZE_Q)(MSIZE_R))
//macros that will be used in the switch statement.
#define CASE_MACRO_R(ignored, nr, np_and_nq) \
BOOST_PP_CAT( \
BOOST_PP_CAT( \
BOOST_PP_CAT( \
BOOST_PP_CAT( \
BOOST_PP_CAT( \
case BOOST_PP_SEQ_ELEM(nr, MSIZE_R): matmult_ , BOOST_PP_SEQ_ELEM(BOOST_PP_SEQ_ELEM(0, np_and_nq),MSIZE_P) \
), \
_ \
), \
BOOST_PP_SEQ_ELEM(BOOST_PP_SEQ_ELEM(1, np_and_nq),MSIZE_Q)\
), \
_\
),\
BOOST_PP_SEQ_ELEM(nr, MSIZE_R)(values_, M2.values_, Mres.values_); break; \
)
#define CASE_MACRO_Q(ignored, nq, np) \
case BOOST_PP_SEQ_ELEM(nq,MSIZE_Q): switch(size2_){ BOOST_PP_REPEAT(BOOST_PP_SEQ_SIZE(MSIZE_R), CASE_MACRO_R, (np)(nq)) default: use_blas=true; } break;
#define CASE_MACRO_P(ignored, np, ignored2) \
case BOOST_PP_SEQ_ELEM(np,MSIZE_P): switch(Mres.size2_){ BOOST_PP_REPEAT(BOOST_PP_SEQ_SIZE(MSIZE_Q), CASE_MACRO_Q, np) default: use_blas=true;} break;
#endif
namespace blas{
  //simple BLAS matrix that uses BLAS calls for rank one and matrix vector.
  //contains size 1 and size 2 that do not have to be the same.
  //std::ostream &operator<<(std::ostream &os, const general_matrix &M); //forward declaration
  template<typename T>class general_matrix{
        public:
    typedef T value_type;
    general_matrix(int size1, int size2){
      //std::cout<<"new general matrix, sizes : "<<size1<<" "<<size2<<std::endl;
      assert(size1>=0 && size2 >=0);
      if(size1*size2>0){
        values_=new T[size1*size2];
        //memset(values_, 0, sizeof(T)*size1*size2);
      }else
        values_=0;
      size1_=size1;
      size2_=size2;
      total_memory_size_=size1*size2;
    }
    general_matrix(int size1){
      //std::cout<<"new general matrix, size : "<<size1<<std::endl;
      assert(size1>=0);
      int size2=size1;
      if(size1>0){
        values_=new T[size1*size2];
        //memset(values_, 0, sizeof(T)*size1*size2);
      }else 
        values_=0;
      size1_=size1;
      size2_=size2;
      total_memory_size_=size1*size2;
    }
    ~general_matrix(){
      if(size1_*size2_>0) //delete of a 0 pointer *is* legal.
        delete[] values_;
    }
    general_matrix(){
      size1_=0;
      size2_=0;
      values_=0;
      total_memory_size_=0;
    }
    template <class M> general_matrix(const M& mat){
      assert(mat.size1()>=0 && mat.size2()>=0);
      size1_=mat.size1();
      size2_=mat.size2();
      total_memory_size_=size1_*size2_;
      if(size1_*size2_>0){
        values_ =new T[size1_*size2_];
        for(int i=0;i<size1_;++i){
          for(int j=0;j<size2_;++j){
            operator()(i,j)=mat(i,j);
          }
        }
      }
      else
        values_=0;
    }
    general_matrix& operator=(const general_matrix &M){
      //std::cout<<"using operator= for matrix: "<<M<<" this is: "<<this<<std::endl;
      if(this==&M) return *this;
      if(size1_!=M.size1_ || size2_!=M.size2_){
        size1_=M.size1_;
        size2_=M.size2_;
        total_memory_size_=size1*size2;
        delete [] values_;
        values_=new T[M.size1_*M.size2_];
      }
      assert(size1_>=0 && size2_ >=0);
      if(size1_*size2_>0){
        memcpy(values_, M.values_, size1_*size2_*sizeof(T));
      }else{
        delete[] values_;
        total_memory_size_=0;
        values_=0;
      }
      return *this;
    }
    general_matrix(const general_matrix&M){
      size1_=M.size1_;
      size2_=M.size2_;
      total_memory_size_=size1_*size2_;
      assert(size1_>=0 && size2_ >=0);
      if(size1_*size2_>0){
        values_=new T[size1_*size2_];
        memcpy(values_, M.values_, size1_*size2_*sizeof(T));
      }else{
        delete[] values_;
        values_=0;
      }
    }
    T* values() { return values_; }
    const T* values() const { return values_; }
    inline T&operator()(const unsigned i, const unsigned j){return *(values_+(i*size2_+j));}
    inline const T&operator()(const unsigned i, const unsigned j) const {return *(values_+(i*size2_+j));}
    //matrix size 
    inline const int &size1()const{return size1_;}
    //matrix size 
    inline const int &size2()const{return size2_;}
    inline void resize(int size1, int size2){
      assert(size1 >=0 && size2>=0);
      if(size1_!=size1 || size2_!=size2){
        size1_=size1; size2_=size2;
        if(size1_*size2_>0){
          if(size1_*size2_>total_memory_size_){
            delete[] values_;
            values_=new T[size1_*size2_];
            total_memory_size_=size1*size2;
            //memset(values_, 0, sizeof(T)*size1*size2);
          }else{
            size1_=size1;
            size2_=size2;
            //memset(values_, 0, sizeof(T)*size1*size2);
          }
        }
        else{
          delete[] values_;
          values_=0;
          total_memory_size_=0;
        }
      }
    }
    //Mres=this*M2
    inline void matrix_right_multiply(const general_matrix<T> &M2,  general_matrix<T> &Mres, T beta=0) const{ //check template specialization!
      //call the BLAS routine for matrix matrix multiplication:
      assert(Mres.size1()==size1_);
      assert(Mres.size2()==M2.size2());
      assert(size2_==M2.size1());
      std::cout<<"please use the optimized version!"<<std::endl;
      abort();
      general_matrix<T> M2t(M2);
      M2t.transpose();
      int s2=Mres.size2();
      for(int i=0;i<size1_*s2;++i){ Mres.values_[i]*=beta; }
      for(int i=0;i<size1_;++i){
        for(int k=0;k<size2_;++k){
          T visk=values_[i*size2_+k];
          T *vjsk=M2t.values_+k;
          T *visj=Mres.values_+i*s2;
          for(int j=0;j<s2;++j){
            *visj+=visk*(*vjsk);
            visj++;
            vjsk+=size2_;
          }
        }
      }
    }
    void clear(){
      if(size1_*size2_>0)
        memset(values_, 0, size1_*size2_*sizeof(T));
    }
    bool operator!=(const general_matrix &M)const {return ! operator==(M);}
    bool operator==(const general_matrix &M)const {
      if(size1_!=M.size1_) return false;
      if(size2_!=M.size2_) return false;
      for(int i=0;i<size1_;++i){
        for(int j=0;j<size2_;++j){
          if(std::fabs(M(i,j) - operator()(i,j))>1.e-15){
            return false;
          }
        }
      }
      return true;
    }
    //multiply matrix by value.
    general_matrix<T> &operator *=(double lambda){
      int inc=1;
      int total_size=size1_*size2_;
      dscal_(&total_size, &lambda, values_, &inc);
      return *this;
    }
    general_matrix<T> operator+(const general_matrix &M2) const{
      general_matrix M(*this);
      for(int i=0;i<size1_;++i){
        for(int j=0;j<size2_;++j){
          M(i,j)+=M2(i,j);
        }
      }
      return M;
    }
    general_matrix<T> &operator+=(const general_matrix &M2){
      for(int i=0;i<size1_;++i){
        for(int j=0;j<size2_;++j){
          operator()(i,j)+=M2(i,j);
        }
      }
      return *this;
    }
    general_matrix<T> operator-(const general_matrix &M2) const{
      //std::cout<<"in operator-"<<std::endl;
      general_matrix M(*this);
      for(int i=0;i<size1_;++i){
        for(int j=0;j<size2_;++j){
          M(i,j)-=M2(i,j);
        }
      }
      //std::cout<<"in operator- returning"<<std::endl;
      return M;
    }
    general_matrix<T> operator-() const{
      general_matrix<T> M2(*this);
      for(int i=0;i<size1_;++i){
        for(int j=0;j<size2_;++j){
          M2(i,j)=-operator()(i,j);
        }
      }
      return M2;
    }
    void swap(general_matrix &M2){
      std::swap(size1_, M2.size1_);
      std::swap(size2_, M2.size2_);
      std::swap(total_memory_size_, M2.total_memory_size_);
      std::swap(values_, M2.values_); //just exchange pointers
    }
    inline double max()const{
      int max_index;
      int total_size=size1_*size2_;
      int inc=1;
      if(total_size==0) return 0;
      if(total_size==1) return std::abs(values_[0]);
      else
        max_index=idamax_(&total_size, values_,&inc);
      return std::fabs(values_[max_index-1]);
    }
    //compute M=D*M, where D is a diagonal matrix represented by
    //the
    //vector. Mij*=Dii
    void left_multiply_diagonal_matrix(const blas::vector &diagonal_matrix){
      assert(size1_==diagonal_matrix.size()); 
      for(int i=0;i<size1_;++i){ //for each col: multiply col. with constant
        //int start=i*size2_;
        for(int j=0;j<size2_;++j){ //for each col: multiply col. with constant
          /*values_[start++]*/operator()(i,j)*=diagonal_matrix(i);
        }
      }     
    }       
    //same, but M=M*D. Slow because of much larger stride.. Mij *=Djj
    void right_multiply_diagonal_matrix(const blas::vector &diagonal_matrix){
      assert(size2_==diagonal_matrix.size()); //must have as many rows as diag has size.
      for(int i=0;i<size1_;++i){ 
        //int start=i*size2_;
        for(int j=0;j<size2_;++j){ 
          /*values_[start++]*/operator()(i,j)*=diagonal_matrix(j);
        }
      }     
    }
    inline void right_multiply(const vector &v1, vector &v2) const{ //perform v2[i]=M[ij]v1[j]
      assert(size1_==size2_);
      //call the BLAS routine for matrix vector multiplication:
      char trans='T';
      double alpha=1., beta=0.;       //no need to multiply a constant or add a vector
      int inc=1;
      dgemv_(&trans, &size1_, &size1_, &alpha, values_, &size1_, &(v1(0)), &inc, &beta, &(v2(0)), &inc);
    }
    inline void right_multiply(const rsvector &v1, rsvector &v2) const{ //perform v2[i]=M[ij]v1[j]
      assert(size1_==size2_);
      //call the BLAS routine for matrix vector multiplication:
      char trans='T';
      double alpha=1., beta=0.;       //no need to multiply a constant or add a vector
      int inc=1;
      dgemv_(&trans, &size1_, &size1_, &alpha, values_, &size1_, &(v1(0)), &inc, &beta, &(v2(0)), &inc);
    }

    inline void multiply_row(const int i, const T&val){ //blas specialization see below
      T *k=values_+i*size2_;
      for(int i=0;i<size1_;++i){
        (*k)*=val; k++;
      }
    }
    inline double trace() const{
      //std::cout<<"trace: sizes: "<<size1_<<" "<<size2_<<std::endl;
      assert(size1_==size2_);
      T trace=0;
      for(int i=0;i<size1_;++i){
        trace+=operator()(i,i);
      }
      return trace;
    }
    void set_to_identity(){
      assert(size1_==size2_);
      memset(values_, 0, size1_*size2_*sizeof(T));
      for(int i=0;i<size1_;++i){
        operator()(i,i)=1;
      }
    }
    void transpose(){
      if(values_==0 || size1_==0 || size2_==0) return;
      T*values=new T[size1_*size2_];
      total_memory_size_=size1_*size2_;
      for(int i=0;i<size1_;++i){
        for(int j=0;j<size2_;++j){
          values[j*size1_+i]=values_[i*size2_+j];
        }
      }
      std::swap(size1_, size2_);
      delete[] values_;
      values_=values;
    }
    general_matrix operator*(const general_matrix M2)const{
      general_matrix M(size1_, M2.size2_);
      matrix_right_multiply(M2, M);
      return M;
    }
    general_matrix &operator*=(const general_matrix M2){
      general_matrix M(size1_, M2.size2_);
      matrix_right_multiply(M2, M);
      swap(M);
      return *this;
    }
    general_matrix operator*(const blas::vector &v)const{
      general_matrix M(*this);
      right_multiply_diagonal_matrix(v);
      return M;
    }
    general_matrix &operator*=(const blas::vector &v){
      right_multiply_diagonal_matrix(v);
      return *this;
    }
    general_matrix &invert(){
      throw(std::logic_error(std::string("you linked the general case for invert. Please use the specializations.")));
    }
        private:
    int size1_; //current size of matrix
    int size2_; //current size of matrix
    int total_memory_size_; //total memory allocated for this matrix
    T *values_; //where the actual values are stored
  };
  
  //this is crap!! reimplement!
#ifndef USE_MATRIX_DISPATCH
  template<> inline void general_matrix<double>::matrix_right_multiply(const general_matrix<double> &M2,  general_matrix<double> &Mres, double beta) const{
    //call the BLAS routine for matrix matrix multiplication:
    //std::cerr<<"sizes are: "<<size1_<<" "<<Mres.size1()<<" "<<Mres.size2()<<std::endl;
    assert(Mres.size1()==size1_);
    assert(Mres.size2()==M2.size2());
    assert(size2_==M2.size1());
    if(size1_ < 4 || size2_ < 4 || M2.size2_ < 4){ //do it by hand! a test showed this to be the optimal size on my intel.
      int s2=Mres.size2();
      double *mres_ij=&Mres(0,0);
      for(int i=0;i<size1_;++i){
        for(int j=0;j<s2;++j){
          (*mres_ij)*=beta;
          const double *thism_ik=&operator()(i,0);
          const double *m2_kj=&M2(0,j);
          for(int k=0;k<size2_;++k){
            (*mres_ij)+=(*thism_ik)*(*m2_kj);
            m2_kj+=s2;
            thism_ik++;
          }
          mres_ij++;
        }
      }
    }
    else{
      double one_double=1.;
      char notrans='N';
      //std::cout<<"invoking dgemm_"<<std::endl;
      dgemm_(&notrans, &notrans, &M2.size2_, &size1_, &size2_, &one_double, M2.values_, &M2.size2_, values_,
             &size2_, &beta, Mres.values_, &M2.size2_);
    }
  }
#else
  template<> inline void general_matrix<double>::matrix_right_multiply(const general_matrix<double> &M2,  general_matrix<double> &Mres, double beta) const{
    //call the BLAS routine for matrix matrix multiplication:
    assert(Mres.size1()==size1_);
    assert(Mres.size2()==M2.size2());
    assert(size2_==M2.size1());
    double one_double=1.;
    char notrans='N';
    bool use_blas=false;
    switch(Mres.size1()){ //boost pp dispatcher for matrix multiplication - see also dispatcher.
        BOOST_PP_REPEAT(BOOST_PP_SEQ_SIZE(MSIZE_P), CASE_MACRO_P, )
      default: use_blas=true; break;
    }
    if(use_blas){
      dgemm_(&notrans, &notrans, &M2.size2_, &size1_, &size2_, &one_double, M2.values_, &M2.size2_, values_,
             &size2_, &beta, Mres.values_, &M2.size2_);
    }
  }
#endif //matrix dispatch
  template<> inline void general_matrix<std::complex<double> >::matrix_right_multiply(const general_matrix<std::complex<double> > &M2,  general_matrix<std::complex<double> >  &Mres, std::complex<double> beta) const{
    //call the BLAS routine for matrix matrix multiplication:
    assert(Mres.size1()==size1_);
    assert(Mres.size2()==M2.size2());
    assert(size2_==M2.size1());
    std::complex<double> one_double=1.;
    char notrans='N';
    if(size1_ < 4 && size2_ < 4 && M2.size2()< 4){ //do it by hand! a test showed this to be the optimal size on my intel.
      general_matrix M2t(M2);
      M2t.transpose();
      int s2=Mres.size2();
      for(int i=0;i<size1_*s2;++i){ Mres.values_[i]*=beta; }
      for(int i=0;i<size1_;++i){
        for(int k=0;k<size2_;++k){
          std::complex<double> visk=values_[i*size2_+k];
          std::complex<double> *vjsk=M2t.values_+k;
          std::complex<double> *visj=Mres.values_+i*s2;
          for(int j=0;j<s2;++j){
            *visj+=visk*(*vjsk);
            visj++;
            vjsk+=size2_;
          }
        }
      }
    }
    else{
      zgemm_(&notrans, &notrans, &M2.size2_, &size1_, &size2_, &one_double, M2.values_, &M2.size2_, values_,
             &size2_, &beta, Mres.values_, &M2.size2_);
    }
  }
  template<> inline double general_matrix<std::complex<double> >::trace() const{
    assert(size1_==size2_);
    std::complex<double> trace=0;
    for(int i=0;i<size1_;++i){
      trace+=operator()(i,i);
    }
    if(trace.imag() > 1.e-10) std::cerr<<"careful, trace is complex: "<<trace<<std::endl;
    return trace.real();
  }
  template<> inline double general_matrix<std::complex<double> >::max()const{
    int max_index;
    int total_size=size1_*size2_;
    int inc=1;
    if(total_size==0) return 0;
    if(total_size==1) return std::abs(values_[0]);
    else
      max_index=izamax_(&total_size, values_,&inc);
    return std::abs(values_[max_index-1]);
  }
  template<typename T> inline std::ostream &operator<<(std::ostream &os, const general_matrix<T> &M){
    os<<"[ ";
    for(int i=0;i<M.size1();++i){
      os<<"[ ";
      for(int j=0;j<M.size2();++j){
        os<<M(i,j)<<" ";
      }
      os<<" ]"<<std::endl;
    }
    os<<"]"<<std::endl;
    return os;
  }
  template<> inline std::ostream &operator<<(std::ostream &os, const general_matrix<std::complex<double> > &M){
    for(int i=0;i<M.size1();++i){
      for(int j=0;j<M.size2();++j){
        os<<M(i,j).real()<<"+i*"<<M(i,j).imag()<<" ";
      }
      os<<std::endl;
    }
    os<<std::endl;
    return os;
  }
  template<> inline std::ostream &operator<<(std::ostream &os, const general_matrix<double > &M){
    for(int i=0;i<M.size1();++i){
      for(int j=0;j<M.size2();++j){
        os<<M(i,j)<<" ";
      }
      os<<std::endl;
    }
    os<<std::endl;
    return os;
  }
  template<> inline general_matrix<std::complex<double> >& general_matrix<std::complex<double> >::invert()
  {
    general_matrix<std::complex<double> > B(size1_, size1_);
    int *ipiv=new int[size1_];
    int info;
    B.set_to_identity();
    lapack::zgesv_(&size1_, &size1_, values_, &size1_, ipiv, &(B(0,0)), &size1_, &info);
    delete[] ipiv;
    if(info){ throw(std::logic_error("in dgesv: info was not zero.")); }
    swap(B);
    return *this;
  }
  template<> inline general_matrix<double >& general_matrix<double >::invert()
  {
    general_matrix<double> B(size1_, size1_);
    int *ipiv=new int[size1_];
    int info;
    B.set_to_identity();
    lapack::dgesv_(&size1_, &size1_, values_, &size1_, ipiv, &(B(0,0)), &size1_, &info);
    delete[] ipiv;
    if(info){ throw(std::logic_error("in dgesv: info was not zero.")); }
    swap(B);
    //std::cout<<"B: "<<B<<" this: "<<*this<<std::endl;
    return *this;
  }

  template<typename T> inline general_matrix<T> operator*(const T &lambda, const general_matrix<T> &M){
    general_matrix<T> M2(M);
    M2*=lambda;
    return M2;
  }
  template<typename T> inline general_matrix<T> operator*(const general_matrix<T> &M, const T &lambda){
    general_matrix<T> M2(M);
    M2*=lambda;
    return M2;
  }
  template<typename T> inline general_matrix<T>operator*(const blas::vector &v, const general_matrix<T> &M){
    general_matrix<T> M2(M);
    M2.left_multiply_diagonal_matrix(v);
    return M2;
  }
  typedef general_matrix<double> double_matrix;
  typedef general_matrix<std::complex< double> > complex_double_matrix;
} //namespace

#endif
