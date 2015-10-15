/*********************************************************
*  Emanuel Gull, 2006 - 2008
*  BLAS Compatibility headers
*  Square Matrix Class - Single Precision
*********************************************************/
#ifndef BLAS_SMATRIX
#define BLAS_SMATRIX

#include "./blasheader.h"
#include "./svector.h"
#include "./resizeable_svector.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sys/types.h> //on osx for uint


namespace blas{
//simple BLAS matrix that uses BLAS calls for rank one and matrix vector.
std::ostream &operator<<(std::ostream &os, const smatrix &M); //forward declaration
class smatrix{
	public:
	smatrix(int size){
		if(size>0)
		  values_=new float[size*size];
		else values_=0;
                size_=size;
		memory_size_=size;
	}
        smatrix(){ 
          values_=0;
          size_=0;
          memory_size_=0;
        }
	~smatrix(){
              if(values_!=0)
		delete[] values_;
	}
	smatrix(int size1, int size2){
	  assert(size1==size2); //square matrix only
	  if(size1>0)
	    values_=new float[size1*size2];
	  else
            values_=0;
	  size_=size1;
	  memory_size_=size1;
	}
        const smatrix& operator=(const smatrix &M){
          if(this==&M) return *this;
		size_=M.size_;
                if(memory_size_!=M.memory_size_){
		  memory_size_=M.memory_size_;
                  if(values_!=0)
		    delete[] values_;
                  if(memory_size_>0){
                    values_=new float[M.memory_size_*M.memory_size_];
                  }else values_=0;
                }
                if(values_!=0)
		  memcpy(values_, M.values_, memory_size_*memory_size_*sizeof(float));
		return *this;
	}
	smatrix(const smatrix&M){
	  size_=M.size_;
	  memory_size_=M.memory_size_;
          if(memory_size_>0){
            values_=new float[M.memory_size_*M.memory_size_];
	    memcpy(values_, M.values_, memory_size_*memory_size_*sizeof(float));
          } else values_=0;
	}
	inline float &operator()(const uint i, const uint j){return *(values_+(i*memory_size_+j));}
	inline const float &operator()(const uint i, const uint j) const {return *(values_+(i*memory_size_+j));}
	//matrix size
	inline int size()const{return size_;}
	//matrix size 
	inline int size1()const{return size_;}
	inline int size2()const{return size_;}
	inline void add_outer_product(const blas::svector &v1, const blas::svector &v2, float alpha=1.){
          add_outer_product(&v1(0), &v2(0), alpha);
        }
	inline void add_outer_product(const blas::rssvector &v1, const blas::rssvector &v2, float alpha=1.){
          add_outer_product(&v1(0), &v2(0), alpha);
        }
	inline void add_outer_product(const float *v1, const float *v2, float alpha=1.){
	  int inc=1;
          if(size_>1){
	    sger_(&size_, &size_, &alpha,v2, &inc, v1, &inc, values_, &memory_size_); 
          }else if(size_==1){
            values_[0]+=alpha*v1[0]*v2[0];
          }else
            return;
	}
        inline float sum(){ 
          float sum=0;
          for(int i=0;i<size_;++i){
            for(int j=0;j<size_;++j){
              sum+=operator()(i,j);
            }
          }
          return sum;
        }
        inline void insert_row_column_last(blas::svector &row, blas::svector &col, float Mkk){
          resize(size_+1);
          int one=1;
          int oldsize=size_-1;
          scopy_(&oldsize, &(col(0)), &one, &(values_[oldsize             ]), &memory_size_); //copy in row (careful: col. major)
          scopy_(&oldsize, &(row(0)), &one, &(values_[oldsize*memory_size_]), &one         );   //copy in column
          operator()(oldsize, oldsize)=Mkk;
        }
        inline void getrow(int k, float *row) const{
          int one=1;
          scopy_(&size_, &(values_[k*memory_size_]), &one, row, &one);
        }
        inline void getcol(int k, float *col) const{
          int one=1;
          scopy_(&size_, &(values_[k]), &memory_size_, col, &one);
        }
        inline void setrow(int k, const float *row){
          int one=1;
          scopy_(&size_, row, &one, &(values_[k*memory_size_]), &one);
        }
        inline void setcol(int k, const float *col){
          int one=1;
          scopy_(&size_, col, &one, &(values_[k]), &memory_size_);
        }
        inline void insert_row_column(blas::svector &row, blas::svector &col, float Mkk, int k){
          resize(size_+1);
          float *new_values=new float[size_*size_];
          //insert a row and a column at position k. could be done more
          //efficiently, I suppose.
          for(int j=size_-2;j>=0;--j){
            int j_new=j<k?j:j+1;
            for(int i=size_-2;i>=0;--i){
              int i_new=i<k?i:i+1;
              new_values[i_new*size_+j_new]=operator()(i, j);
            }
            new_values[k*size_+j_new]=row(j);
            new_values[j_new*size_+k]=col(j);
          }
          new_values[k*size_+k]=Mkk;
          delete[] values_;
          values_=new_values;
          memory_size_=size_;
        }
        //delete column k.
        inline void remove_row_column(int c){
          int new_size=size_-1;
          blas::smatrix new_matrix(new_size);
          for(int i=0, k=0;i<new_size;++i,++k){
            if(k==c) ++k;
            for(int j=0,l=0;j<new_size;++j,++l){
              if(c==l) ++l;
              new_matrix(i,j)=operator()(k, l);
            }
          }
          this->swap(new_matrix);
        }
        //delete last column
        inline void remove_row_column_last(){
          size_--;
        }
        //swap two columns:
        inline void swap_row_column(int c1, int c2){
          if(c1==c2) return;
          int one=1;
          sswap_(&size_, &(values_[c1]), &memory_size_, &(values_[c2]), &memory_size_);
          sswap_(&size_, &(values_[c1*memory_size_]), &one, &(values_[c2*memory_size_]), &one);
        }
        template<class mat> smatrix convert_from(const mat &M){ //convert dense or sparse matrices (or any other with  (i,j) and size1() to blas matrix.
	  assert(M.size1()==M.size2());
          resize_nocopy(M.size1());
          for(int i=0;i<(int)M.size1();++i){
            for(int j=0;j<(int)M.size1();++j){
              operator()(i,j)=M(i,j);
            }
          }
          return *this;
        }
        inline void right_multiply(const svector &v1, svector &v2) const{ //perform v2[i]=M[ij]v1[j]
		//call the BLAS routine for matrix vector multiplication:
		char trans='T';
		float alpha=1., beta=0.;	//no need to multiply a constant or add a vector
		int inc=1;
		sgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
	}
        inline void right_multiply(const rssvector &v1, rssvector &v2) const{ //perform v2[i]=M[ij]v1[j]
		//call the BLAS routine for matrix vector multiplication:
		char trans='T';
		float alpha=1., beta=0.;	//no need to multiply a constant or add a vector
		int inc=1;
		sgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);
	}
	inline void left_multiply(const svector &v1, svector &v2) const{ //perform v2[i]=v1[j]M[ji]
                //call the BLAS routine for matrix vector multiplication:
                char trans='N';
                float alpha=1., beta=0.;       //no need to multiply a constant or add a vector
                int inc=1;
                sgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);

	}
	inline void left_multiply(const rssvector &v1, rssvector &v2) const{ //perform v2[i]=v1[j]M[ji]
                //call the BLAS routine for matrix vector multiplication:
                char trans='N';
                float alpha=1., beta=0.;       //no need to multiply a constant or add a vector
                int inc=1;
                sgemv_(&trans, &size_, &size_, &alpha, values_, &memory_size_, v1.values_, &inc, &beta, v2.values_, &inc);

	}
        //Mres=this*M2
	inline void matrix_right_multiply(const smatrix &M2,  smatrix &Mres) const{
	  //call the BLAS routine for matrix matrix multiplication:
	  float one_float=1.;
	  float zero_float=0.;
	  char notrans='N';
	  sgemm_(&notrans, &notrans, &size_, &size_, &size_, &one_float, M2.values_, &memory_size_, values_, &memory_size_, &zero_float, Mres.values_, &Mres.memory_size_);
	}
	inline float trace()const{
	  float tr=0;
	  for(int i=0;i<size_;++i) tr+=operator()(i,i);
	  return tr;
	}
        inline float determinant() const{
          //the simple ones...
          if(size_==0) return 1;
          if(size_==1) return values_[0];
          if(size_==2) return values_[0]*values_[3]-values_[1]*values_[2];
          int info=0;
          int *ipiv = new int[size_];
          assert(size_==memory_size_); //otherwise think about plugging ing memory size to lda
          
          blas::smatrix det_matrix(*this); //will be changed...
          blas::smatrix identity(*this);
          identity.set_to_identity();
          //LU factorization
          lapack::sgesv_(&size_, &size_, det_matrix.values_, &memory_size_, ipiv, identity.values_, &memory_size_,&info);
          if(info < 0) {
              std::cout << "LAPACK ERROR IN APPLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
              std::cout << "INFO:" << info << std::endl;
          } else if(info > 0){
            //check dgesv page: that means that we haven an exactly singular
            //matrix and the det is therefore =0: 
            return 0.;
          }
          //compute the determinant:
          float det=1.;
          //CAREFUL when pivoting: fortran uses array indexing starting
          //from one. But we're in C here -> 'one off' error
          for(int i=0;i<size_;++i){
            if(ipiv[i]-1!=i){
              det*=-det_matrix(i,i);
            }
            else{
              det*=det_matrix(i,i);
            }
          }
          delete[] ipiv;
          return det;
        }
	inline void transpose(){
	  for(int i=0;i<size_;++i){
	    for(int j=0;j<i;++j){
	      std::swap(operator()(i,j), operator()(j,i));
	    }
	  }
	}
	void set_to_identity(){
	  for(int i=0;i<size_;++i){
	    for(int j=0;j<size_;++j){
	      operator()(i,j)=0;
	    }
	    operator()(i,i)=1;
	  }
	}
	void clear(){
          for(int i=0;i<size_;++i){
            for(int j=0;j<size_;++j){
              operator()(i,j)=0;
            }
	  }
	}
	bool operator!=(const smatrix &M)const {return ! operator==(M);}
	bool operator==(const smatrix &M)const {
	  if(size_!=M.size_) return false;
	  for(int i=0;i<size_;++i){
	    for(int j=0;j<size_;++j){
	      if(std::fabs(M(i,j) - operator()(i,j))>1.e-5){
	        return false;
	      }
	    }
	  }
	  return true;
	}
	//multiply matrix by value.
	smatrix &operator *=(float lambda){
	  //int inc=1;
	  //int total_size=memory_size_*memory_size_;
          //dscal_(&total_size, &lambda, values_, &inc);
          for(int i=0;i<size_;++i){
            for(int j=0;j<size_;++j){
              operator()(i,j)*=lambda;
            }
          }
	  return *this;
	}
        smatrix operator+(const smatrix &M2) const{
          smatrix Msum(*this);
          Msum+=M2;
          return Msum;
        }
        smatrix &operator+=(const smatrix &M2) {
		assert(size_==M2.size_);
          for(int i=0;i<size_;++i){
            for(int j=0;j<size_;++j){
              operator()(i,j)+=M2(i,j);
            }
          }
          return *this;
        }
        smatrix operator-(const smatrix &M2) const{
          smatrix Msum(*this);
          for(int i=0;i<size_;++i){
            for(int j=0;j<size_;++j){
              Msum(i,j)=operator()(i,j)-M2(i,j);
            }
          }
          return Msum;
        }
	void resize(int size1, int size2){
	  assert(size1==size2);
	  resize(size1);
	}
	void resize_nocopy(int new_size){
          if(size_!=0)
	    delete[] values_;
	  values_=new float[new_size*new_size];
	  memset(values_, 0, sizeof(float)*new_size*new_size);
	  size_=new_size;
	  memory_size_=new_size;
	}
	void resize(int new_size){
          if(new_size==size_) return;
          if(new_size<size_){ //down is easy
            if((int)new_size < (int)(memory_size_-30) && (int) new_size > 10){
              float *new_values_=new float[new_size*new_size];
              for(int i=0;i<new_size;++i){
                memcpy(new_values_+i*new_size, values_+i*memory_size_, sizeof(float)*new_size); //for each row: copy the entire row.
              }
              delete[] values_;       //free memory
              values_=new_values_;    //let the matrix point to the new memory location.  
              size_=new_size;
              memory_size_=new_size;
            }
            size_=new_size;
            return;
          } else if(new_size<= memory_size_){ //up is easy as long as we don't have to allocate new memory
            size_=new_size;
          } else{ //get new memory */
          float *new_values_=new float[new_size*new_size];
            for(int i=0;i<size_;++i){
              memcpy(new_values_+i*new_size, values_+i*memory_size_, sizeof(float)*size_); //for each row: copy the entire row.
	    }
            delete[] values_;       //free memory
            values_=new_values_;    //let the matrix point to the new memory location.
            size_=new_size;
            memory_size_=new_size;
          }
	}
        float max() const{
          float rowmax[size_];
          int rowmax_index, max_index, inc=1;
          for(int i=0;i<size_;++i){
            rowmax_index=blas::isamax_(&size_, values_+i*memory_size_,&inc);
            rowmax[i]=*(values_+i*memory_size_+rowmax_index-1); //fortran convention: start counting from one
          }
          max_index=blas::isamax_(&size_, rowmax ,&inc);
          return std::abs(rowmax[max_index-1]);
        }
	void swap(smatrix &M2){
	  std::swap(size_, M2.size_);
	  std::swap(memory_size_, M2.memory_size_);
	  std::swap(values_, M2.values_); //just exchange pointers
	}
        void eigenvalues_eigenvectors_symmetric( svector &eigenvalues, smatrix &eigenvectors) const{
          //perform ssyev call (LAPACK)
          eigenvectors=*this;
          assert(memory_size_==size_);
          int lwork=-1;
          int info;
          float work_size;
          char jobs='V';
          char uplo='L';
          //get optimal size for work
          lapack::ssyev_(&jobs, &uplo, &size_, eigenvectors.values_, &size_, &(eigenvalues(0)),&work_size, &lwork, &info);
          lwork=(int)work_size;
          float work[lwork];
          lapack::ssyev_(&jobs, &uplo, &size_, eigenvectors.values_, &size_, &eigenvalues(0),work, &lwork, &info); 
          /*matrix eigenvectors_trans(eigenvectors);
          eigenvectors_trans.transpose();
          std::cout<<*this<<std::endl;
          matrix diag(size_);
          diag.clear();
          for(int i=0;i<size_;++i){
            diag(i,i)=eigenvalues(i);
          }
          std::cout<<eigenvectors<<std::endl;
          std::cout<<diag<<std::endl;
          std::cout<<(eigenvectors_trans*diag*eigenvectors)<<std::endl;
          std::cout<<(eigenvectors*(*this)*eigenvectors_trans)<<std::endl;*/
        }
        //compute M=D*M, where D is a diagonal matrix represented by the
        //vector.
        void multiply_diagonal_matrix(const blas::svector &diagonal_matrix){
          int inc=1;
          for(int i=0;i<size_;++i){
            sscal_(&size_, &diagonal_matrix(i), values_+i*memory_size_, &inc);
          }
        }
        //same, but M=M*D. Slow because of much larger stride.
        void right_multiply_diagonal_matrix(const blas::svector &diagonal_matrix){
          for(int i=0;i<size_;++i){
            sscal_(&size_, &diagonal_matrix(i), values_+i, &memory_size_);
          }
        }
	private:
	int size_; //current size of matrix
	int memory_size_; //current size of matrix
	float *values_; //where the actual values are stored
};

inline svector operator*(const smatrix &M, const svector &v1){
	assert(v1.size()==M.size());
	svector vres(v1.size());
        if(M.size()>0)
	  M.right_multiply(v1, vres);
	return vres;
}
inline svector operator*(const svector &v1, const smatrix &M){
	assert(v1.size()==M.size());
	svector vres(v1.size());
        if(M.size()>0)
          M.left_multiply(v1, vres);
	return vres;
}
inline smatrix operator*(const smatrix &M, const float lambda){
        smatrix M2=M;
        M2*=lambda;
        return M2;
}
inline smatrix operator*(const float lambda, const smatrix &M){
        smatrix M2=M;
        M2*=lambda;
        return M2;
}
inline smatrix operator*(const smatrix M, const smatrix M2){
          smatrix Mres(M.size());
          M.matrix_right_multiply(M2, Mres);
          return Mres;
}
inline std::ostream &operator<<(std::ostream &os, const smatrix &M){
    os<<"[ ";
  for(int i=0;i<M.size();++i){
    os<<"[ ";
    for(int j=0;j<M.size();++j){
      os<<M(i,j)<<" ";
    }
    os<<" ]"<<std::endl;
  }
  os<<"]"<<std::endl;
  return os;
}
} //namespace

#endif
