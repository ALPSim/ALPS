/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef __DMTK_LANCZOS_CC__
#define __DMTK_LANCZOS_CC__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include "range.h"
#include "vector.h"
#include "matrix.h"

using namespace std;

namespace dmtk
{

#define SIGN(a,b) ((b) < 0 ? -fabs(a) : fabs(a))
#define SWAP(a,b) { itemp=(a);(a)=(b);(b)=itemp; }

template<class T>
class AuxIndex
{
  public:
    T v;
    size_t index;

    AuxIndex():v(T()), index(0) {}
    AuxIndex(const T& _v, const size_t &i): v(_v), index(i) {}
    AuxIndex & operator=(const AuxIndex& o) { v = o.v; index = o.index; return *this; } 

    bool operator>(const AuxIndex &o) const
      {
        if(v > o.v) return true;
        if(v == o.v && index > o.index) return true;
        return false;
      }
    bool operator<(const AuxIndex &o) const
      {
        if(v < o.v) return true;
        if(v == o.v && index < o.index) return true;
        return false;
      }
    bool operator>=(const AuxIndex &o) const
      {
        return (v >= o.v && index >= o.index); 
      }
    bool operator<=(const AuxIndex &o) const
      {
        return(v <= o.v && index <= o.index);
      }
    bool operator==(const AuxIndex &o) const
      {
        return (v == o.v);
      }
    bool operator!=(const AuxIndex &o) const
      {
        return (v != o.v);
      }
};


template <class T>
int sort2(size_t n, Vector<T>& arr, bool ascending = true)
{
  size_t j,l;
  T itemp;
  int nswaps = 0;
  bool changed = true;
  T a;

  if(n == 0){
    cout << "*** WARNING: sort width (n == 0) \n";
    return 0;
  }

  Vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 

  while(changed){
    changed = false;
    for(l = 0; l < n-1; l++){
      if(aux[l] > aux[l+1]){
        SWAP(aux[l],aux[l+1]); 
        nswaps++;
        changed = true;
      }
    }
  }

  if(!ascending && n > 1)
    for(j = 0; j <= n/2-1; j++)
      if(aux[j] != aux[n-j-1]) {
        SWAP(aux[j],aux[n-j-1]);
        nswaps++;
      }

  for(int i = 0; i < n; i++) arr[i] = aux[i].v;  
  return nswaps;
}

template <class T, class A>
int indexx2(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true)
{
  size_t j,l,itemp;
  bool changed = true;
  int nswaps = 0;

  Vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 

  for (j=0;j<n;j++) indx[j]=j;
  while(changed){
    changed = false;
    for(l = 0; l < n-1; l++){
      if(aux[indx[l]] > aux[indx[l+1]]){
        SWAP(indx[l],indx[l+1]); 
        nswaps++;
        changed = true;
      }
    }
  }

  if(!ascending && n > 1)
    for(j = 0; j <= n/2-1; j++){
      SWAP(indx[j],indx[n-j-1]);
      nswaps++;
    }

  return nswaps;
}

template <class T>
void sort(size_t n, Vector<T>& arr, bool ascending = true)
{
  if (n < 2) return;

  T itemp;
                                                                              
  std::vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 

  std::sort(aux.begin(), aux.end());

  ra = &aux[0];
  for(int i = 0; i < n; i++, ra++) arr[i] = ra->v;  

  if(!ascending && n > 1)
    for(int j = 0; j <= n/2-1; j++)
      if(arr[j] != arr[n-j-1]) SWAP(arr[j],arr[n-j-1]);
}

template<class T, class A>
void indexx(size_t n, const A& arr, Vector<size_t>& indx, bool ascending = true)
{
  if (n == 1) { indx[0] = 0; return; }

  size_t itemp;
                                                                              
  std::vector<AuxIndex<T> > aux(n);
  AuxIndex<T> *ra = &aux[0];
  for(int i = 0; i < n; i++) { ra->v = arr[i]; ra->index = i; ra++; } 
  std::sort(aux.begin(), aux.end());

  ra = &aux[0];
  for(int i = 0; i < n; i++, ra++) indx[i] = ra->index;  

  if(!ascending && n > 1)
    for(int j = 0; j <= n/2-1; j++)
      SWAP(indx[j],indx[n-j-1]);
}

void tqli(dmtk::Vector<double>& d, dmtk::Vector<double>& e, 
          int n, dmtk::Matrix<double>& z, bool calc_vectors = false)
{
  char calc = calc_vectors ? 'V' : 'N';
  int lwork = calc_vectors ? 1 + 4*n + n*n : 1;
  Vector<double> work(lwork);
  int liwork = calc_vectors ? 3 + 5*n : 1;
  Vector<int> iwork(liwork);
  int info;

  FORTRAN_ID(dstevd)(calc,n,d.array(),e.array(),z.array(),n,work.array(),lwork,iwork.array(),liwork,info);
  if(info != 0) cout << "ERROR in dstevd\n";
}

//////////////////////////////////////////////////////////////////////

template<class T, class A, class B>
void lanczos(A& m,
             Vector<B>& gs, const B& seed, 
             Vector<double>& ener, int nvectors,
             dmtk::Vector<double>& a, dmtk::Vector<double>& b, 
             int &maxiter, double tol,
             bool use_seed, bool calc_gs, const char* vectors_file, bool force_maxiter = false)
{  
  B x1(gs(0)), x2(gs(0)), aux(gs(0));
  dmtk::Vector<double> d(101), e(101), xc(101);
  dmtk::Matrix<double> z(101,101);
  uint j, col, iter = 0;
  // uint i;
  double eini, e0;
  int control_max = maxiter;
  maxiter = 0;

  if(control_max == 0) { gs[0] = T(1); maxiter = 1; return; }
#ifndef USE_LANCZOS_VECTORS
  ofstream outputfile(vectors_file,std::ios::out|std::ios::binary);
  if(!outputfile) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";
#endif
 
  e0 = 10000.0f;
  x1 = x2 = aux = T(0);

  while(true) // Two iterations: 1) calculate energy 2) build gs
  { 
    a = 0.0f;
    b = 0.0f;
    if(use_seed){
      x1 = seed;
    } else { 
#ifdef USE_LANCZOS_BASIC
      x1 = 1.;
//      x1 = 0;
//      x1[0] = 1.;
#else
      x1.randomize();
#endif
    } 

    iter = 0;

    b(0) = sqrt(real(dot_product(x1,x1)));
    x1 = x1 / T(b(0));
    x2 = T(0);
    b(0) = 1.;


    uint nmax = std::min(99, (int)gs(0).size());
    for(iter = 1; iter <= nmax; iter++){ // Lanczos iteration
      eini = e0;

      if(b(iter - 1) != 0.){
        aux = x1;
        x1 = -T(b(iter-1)) * x2;
        x2 = aux / T(b(iter-1));
      }

//---------------------------------------------------------------
//******* REORTHOGONALIZATION
/*
      outputfile.close();
      ifstream inputfile(vectors_file,std::ios::in|std::ios::binary);
      if(!inputfile) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";

      B aux1 = x2, aux2 = x2;
      T overlap;
      if(iter >= 1){ 
         for(int i = 1; i < iter; i++){
           aux2.read(inputfile);
           overlap = product(aux1,aux2);
           x2 -= overlap*aux2;
         }
         x2 /= product(x2,x2);

//         inputfile.close();
//         ifstream inputfile(vectors_file,std::ios::in|std::ios::binary);
//         for(int i = 1; i < iter; i++){
//           aux2.read(inputfile);
//           overlap = product(x2,aux2);
//           cout << i << " " << overlap << endl;
//         }

      }

      inputfile.close();
      ofstream outputfile(vectors_file,std::ios::out|std::ios::binary|std::ios::app);
      if(!outputfile) cerr << "*** ERROR 1: Lanczos could not open vectors.dat\n";
*/
//---------------------------------------------------------------
      aux = product(m,x2);
      x1 = x1 + aux;
      T xa = dot_product(x1, x2);
      a(iter) = real(xa);
      x1 = x1 - x2*T(a(iter));
      b(iter) = sqrt(real(dot_product(x1, x1)));
   
//      cout << "Iter=" << iter << " " << xa << " " << a(iter) << " " << b(iter) << endl;

#ifndef USE_LANCZOS_VECTORS
      if(calc_gs) x2.write(outputfile);
#endif

      if(maxiter > 0) {  // We are building the ground state;

#ifdef USE_LANCZOS_VECTORS
        gs = gs + xc(iter-1) * x2;
        if(iter == maxiter) return; 
#endif
      } else{  // we are calculating energy

        if(iter >= 2){
          d(Range(0,iter-1)) = a(Range(1,iter));
          e(Range(0,iter-1)) = b(Range(1,iter));
          // call tqli without eigenvectors
          tqli(d, e, iter, z, false);
  
          e0 = 10000.0f;
          for(j = 0; j < iter; j++){
            if(d(j) < e0) {
              e0 = d(j);
              col = j;
            }
          }

#ifdef DMTK_DEBUG
          cout << setprecision(12) << "Iter = " << iter << "  Ener = " << e0 << endl;
#endif // DMTK_DEBUG
          if((force_maxiter && iter >= control_max) || (!force_maxiter && (iter >= gs(0).size()-1 || iter == 99 || fabs(b(iter)) < tol || fabs(eini-e0) <= tol || iter >= control_max)))
          { 
             // converged
             if(fabs(b(iter)) < tol) cout << "b(iter) < tol " << b(iter) << endl;
             cout << setprecision(12) << "E0 = " << e0 << endl;
             maxiter = iter;
             if(!calc_gs) return; // We return with ground states energy
             z = I<double>(); //identity;
             d(Range(0,iter-1)) = a(Range(1,iter));
             e(Range(0,iter-1)) = b(Range(1,iter));
             // call tqli with eigenvectors
             tqli(d, e, iter, z, true);
             xc = z(col, Range(0,iter));

#ifdef USE_LANCZOS_VECTORS
             break; // Exit Lanczos iteration. Re-start for ground state.
#endif
             outputfile.flush();
             outputfile.close();
         
             Vector<size_t> indx(maxiter); 
//             d(Range(0,maxiter-1)) = d(Range(1,maxiter));
             indexx<double,Vector<double> >(maxiter, d, indx, true);

#ifdef DMTK_DEBUG
             for(int n = 0; n < maxiter; n++)
               cout << "LEVEL " << n << " " << " " << d(indx(n)) << endl;
#endif // DMTK_DEBUG

             for(int n = 0; n < nvectors; n++){
               ifstream inputfile(vectors_file,std::ios::in|std::ios::binary);
               if(!inputfile) cerr << "*** ERROR 2: Lanczos could not open vectors.dat\n";
               gs(n) = T(0);
               ener(n) = d(indx(n));
//               cout << n << " " << ener(n) << endl;
               for(int i = 0; i < maxiter; i++){
                 x2.read(inputfile);
                 gs(n) += T(z(indx(n),i)) * x2;
               }
               inputfile.close();

             }

             return;
          }
        } // diagonalization of tridiagonal matrix

      } 
    } // Lanczos iteration
  } // main iteration

}

template<class T, class A, class B>
void verif_hamiltonian(A& m, const B& gs)
{
  Matrix<T> h(gs.size(),gs.size());
  B aux(gs);

  cout << "***** VERIFYING HAMILTONIAN *****\n";   
  for(int col = 0; col < gs.size(); col++){
    aux = T(0);
    aux(col) = T(1);
    h.column(col) = product(m,aux);
    for(int row = 0; row < gs.size(); row++){
      cout << col << " " << row << " " << h.column(col)(row) << endl;
    }
  }
  for(int col = 0; col < gs.size(); col++){
    for(int row = col+1; row < gs.size(); row++){
      if(fabs(h(row,col)-h(col,row)) > 1.e-5) 
        cout << "ERROR: " << col << " " << row << " " << h(col,row) << " " << h(row,col) << endl;
    }
  }
}
 
#undef SWAP
#undef SIGN

} // namespace dmtk

#endif // __DMTK_LANCZOS_CC__
