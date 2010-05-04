/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

#ifndef DMFT_QMC_WEAK_COUPLING_RESIZABLE_MATRIX_H
#define DMFT_QMC_WEAK_COUPLING_RESIZABLE_MATRIX_H

#include<fstream>
#include "operator.hpp"
//this implements all functinons necessary for an inverse 'M' - matrix.

//BLAS calls
extern "C" void dgemv_(const void *trans, const void *size1, const void *size2, const void *alpha, const void *values_, const void *memory_size,
                       const void *v1, const void *inc, const void *beta, void *v2, const void *incb);
extern "C" void dgemm_(const void *transa, const void *transb, const void *m, const void *n, const void *k, 
		       const void *alpha, const void *A, const void *lda,
		       const void *B, const void *ldb, const void *beta, const void *C, const void *ldc);
extern "C" void dger_(const void *m, const void* n, const void* alpha,const void*  x, const void* incx, const void*y, const void*  incy,
                      const void* a,const void* lda) ;
extern "C" int idamax_(const unsigned int *n,const double*  x,const unsigned int* incx);
extern "C" int dnrm2_ (const unsigned int *n,const double*  x,const unsigned int* incx);
extern "C" double ddot_(const int *size,const double * v1,const int *inc,const double *v2,const int *inc2);
double inner_prod(const double* v1, const double *v2, const int size);
extern "C" void dscal_(const unsigned int *size, const double *alpha, double *v, const int *inc);
void scale(const double alpha, double *v, const unsigned int size);


inline double vector_max(const double *v, const unsigned int noperators)
{
  unsigned int inc=1;
  unsigned int max_index=idamax_(&noperators, v ,&inc);
  return fabs(v[max_index-1]); //fortran indexing
}


class resizeable_matrix
{
public:

  resizeable_matrix(const unsigned int ms=1)
  {
    values_=new double[ms*ms];
    size_=0;
    memory_size_=ms;
  }

  ~resizeable_matrix()
  {
    delete[] values_;
  }
  
  const resizeable_matrix& operator=(const resizeable_matrix &M)
  {
    memory_size_=M.memory_size_;
    size_=M.size_;
    max_size_=M.max_size_;
    memcpy(values_, M.values_, memory_size_*memory_size_*sizeof(double));
    creators_=M.creators_;
    annihilators_=M.annihilators_;
    alpha_=M.alpha_;
    return *this;
  }
  
  resizeable_matrix(const resizeable_matrix &M)
  {
    values_=new double[M.memory_size_*M.memory_size_];
    operator=(M);
  }
  
  inline double &operator()(const unsigned int i, const unsigned int j){return *(values_+(i*memory_size_+j));}
  inline const double &operator()(const unsigned int i, const unsigned int j) const {return *(values_+(i*memory_size_+j));}

  void resize(const unsigned int new_size)
  {
    if(new_size<=size_){ //down is easy
      if((int)new_size < (int)(size_-5) && (int) new_size > 10){
        double *new_values=new double[new_size*new_size];
        for(unsigned int i=0;i<size_;++i){
          memcpy(new_values+i*new_size, values_+i*memory_size_, sizeof(double)*new_size); //for each row: copy the entire row.
        }
        delete[] values_;       //free memory
        values_=new_values;    //let the matrix point to the new memory location.
        size_=new_size;
        memory_size_=new_size;
      }
      size_=new_size;
      return;
    } else if(new_size<= memory_size_){ //up is easy as long as we don't have to allocate new memory
      size_=new_size;
    } else{ //get new memory
      double *new_values=new double[new_size*new_size];
      for(unsigned int i=0;i<size_;++i){
        memcpy(new_values+i*new_size, values_+i*memory_size_, sizeof(double)*size_); //for each row: copy the entire row.
      }
      delete[] values_;	//free memory
      values_=new_values;	//let the matrix point to the new memory location.
      size_=new_size;
      memory_size_=new_size;
    }
  }
  
  inline unsigned int size()const{return size_;}
  
  inline void add_outer_prod(const double *v1, const double *v2)
  {
  //smart people call this 'rank one update'.
    double alpha=1.;
    unsigned int inc=1;
    dger_(&size_, &size_, &alpha, v2, &inc, v1, &inc, values_, &memory_size_); 
  }

  inline void right_multiply(const double *v1, double *v2)
  { //perform v2[i]=M[ij]v1[j]
    //call the BLAS routine for matrix vector multiplication:
    char trans='T';
    double alpha=1., beta=0.;	//no need to multiply a constant or add a vector
    int inc=1;
    dgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1, &inc, &beta, v2, &inc);
  }

  inline void matrix_right_multiply(const double *M1, double *M2, const unsigned int columns)
  { //perform M2(i,j)=sum_k M(i,k)*M1(k,j)
    //call the BLAS routine for matrix matrix multiplication:     
    char transa='N';
    char transb='N';
    double alpha=1.;
    double beta=0.;
    dgemm_(&transa, &transb, &columns, &size_, &size_, &alpha, M1, &columns, values_, &memory_size_, &beta, M2, &columns);
  }
  
  void left_multiply(const double *v1, double *v2)
  { //perform v2[i]=v1[j]M[ji]
    //call the BLAS routine for matrix vector multiplication:
    char trans='N';
    double alpha=1., beta=0.;       //no need to multiply a constant or add a vector
    int inc=1;
    dgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1, &inc, &beta, v2, &inc);
  }
  
  bool check()
  {
    bool pass=true;
    for(unsigned int i=0;i<size_;++i){
      for(unsigned int j=0;j<size_;++j){
        if(operator()(i,j)>2){
          std::cout<<"M is really large: "<<i<<" "<<j<<" "<<operator()(i,j)<<" "<<size_<<std::endl;
          pass=false;
        }
      }
    }
    return pass;
  }
  
  double max() const
  {
    double *rowmax=new double[size_];
    unsigned int rowmax_index, max_index, inc=1;
    for(unsigned int i=0;i<size_;++i){
      rowmax_index=idamax_(&size_, values_+i*memory_size_,&inc);
      rowmax[i]=*(values_+i*memory_size_+rowmax_index-1); //fortran convention: start counting from one
    }
    max_index=idamax_(&size_, rowmax ,&inc);
    double m=fabs(rowmax[max_index-1]);
    delete[] rowmax;
    return m;

  }
  
  // return functions
  void clear(){ size_=0;}
  inline unsigned int memory_size() const {return memory_size_;}
  
  std::vector<creator> &creators(){ return creators_;}
  const std::vector<creator> &creators() const{ return creators_;}
  std::vector<annihilator> &annihilators(){ return annihilators_;}
  const std::vector<annihilator> &annihilators()const{ return annihilators_;}
  std::vector<double> &alpha(){ return alpha_;}
  const std::vector<double> &alpha() const{ return alpha_;}
  
  
private:
  std::vector<creator> creators_;         //an array of creation operators c_dagger corresponding to the row of the matrix
  std::vector<annihilator> annihilators_; //an array of to annihilation operators c corresponding to the column of the matrix
  std::vector<double> alpha_;             //an array of doubles corresponding to the alphas of Rubtsov for the c, cdaggers at the same index.
  
  unsigned int memory_size_; //size of matrix for which we have memory allocated
  unsigned int size_; //current size of matrix
  unsigned int max_size_; //max size of matrix
  double *values_; //where the actual values are stored
};

std::ostream & operator<<(std::ostream &os, const resizeable_matrix &M);

#endif
