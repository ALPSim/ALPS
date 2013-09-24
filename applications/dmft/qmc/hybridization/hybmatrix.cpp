/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
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

#include "hybmatrix.hpp"

//this was changed with update 51!
//nomenclature: the c        is at the segment end  , so the time for c        is new_segment->t_end_
//              the c_dagger is at the segment start, so the time for c_dagger is new_segment->t_start_

//the direct matrix F has the structure F_ij=Delta(\tau_i_end - \tau_j_start) // (check)
//that means that the creation operators enter the rows (the i-s), and the annihilation operators form columns (the j-s) //unchanged (check)

//Careful with the ordering: the ORDERED matrix has creation operators and annihilation operators ordered by start and end time,
//in case of a wraparound segment the wrapping c is the last entry.
//this matrix can always be turned into an ordered matrix by running rebuild_ordered_matrix.
//in general the matrix is not ordered (that would be inefficient), but permutation_sign keeps track of the shift of rows
//and columns.


//compute the hybridization weight change when an operator pair is inserted
double hybmatrix::hyb_weight_change_insert(const segment &new_segment, int orbital, const hybfun &Delta){
  Q.resize(size());
  R.resize(size());
  PinvQ.resize(size());
  //column  Delta_i,last
  for(hyb_map_t::const_iterator it=c_index_map_.begin();it!=c_index_map_.end();++it){
    Q[it->second]=Delta.interpolate(it->first-new_segment.t_start_, orbital); //this is the new column Q
  }
  
  //row Delta_last,i
  for(hyb_map_t::const_iterator it=cdagger_index_map_.begin();it != cdagger_index_map_.end();++it){
    R[it->second]=Delta.interpolate(new_segment.t_end_-it->first, orbital);  //this is the new row R
  }
  S=Delta.interpolate(new_segment.t_end_-new_segment.t_start_, orbital);  //this is the entry S
  //std::cout<<clblue<<"last entry S is: "<<S<<cblack<<std::endl;
  
  S_tilde_inv=S;
  fortran_int_t s=size();
  if(s>0){
    right_multiply(Q, PinvQ); //dgemv
    fortran_int_t inc=1;
    S_tilde_inv-=FORTRAN_ID(ddot)(&s, &(R[0]),&inc,&(PinvQ[0]),&inc);
  }
  //a -1 from the anticommutator from the wraparound segment
  if(new_segment.t_end_<new_segment.t_start_){
    weight_ratio_=-S_tilde_inv;
  }else{
    weight_ratio_=S_tilde_inv;
  }
  return weight_ratio_;
}

//actually insert an operator pair and change the configuration
void hybmatrix::insert_segment(const segment &new_segment, int orbital){
  //std::cout<<clred<<"upon entering insert segment: "<<*this<<cblack<<std::endl;
  //consistency_check();
  //enlarge the M-matrix by one
  int last=size();
  resize(size()+1);
  
  //last element
  operator()(last, last)=1./S_tilde_inv;
  
   fortran_int_t sm1=size()-1;
  if(sm1>0){ //this is exactly the content of the loops above, in dger/dgemv blas calls.
    char trans='T', notrans='N';
    double alpha=-1./S_tilde_inv, beta=0.;
    fortran_int_t inc=1;
    fortran_int_t ms=memory_size();
    FORTRAN_ID(dgemv)(&  trans, &sm1, &sm1, &alpha, &(operator()(0,0)), &ms, &(Q[0]), &inc, &beta, &(operator()(0,last)), &ms);
    FORTRAN_ID(dgemv)(&notrans, &sm1, &sm1, &alpha, &(operator()(0,0)), &ms, &(R[0]), &inc, &beta, &(operator()(last,0)), &inc);
    alpha=S_tilde_inv;
    FORTRAN_ID(dger)(&sm1, &sm1, &alpha,&(operator()(last,0)), &inc, &(operator()(0,last)), &ms, &(operator()(0,0)), &ms);
  }
  
  // add the new segment times:
  cdagger_index_map_.insert(std::make_pair(new_segment.t_start_, last));
  c_index_map_      .insert(std::make_pair(new_segment.t_end_  , last));
  
  //keep track of the wraparound sign
  if(new_segment.t_start_>new_segment.t_end_){
    permutation_sign_*=-1.;
  }
  //std::cout<<clred<<*this<<cblack<<std::endl;
  //consistency_check();
  
}
//compute the hybridization weight change when an operator pair is removed
double hybmatrix::hyb_weight_change_remove(const segment &new_segment, int orbital, const hybfun &Delta){
  //std::cout<<clgreen<<"proposing to remove the segment: "<<new_segment<<cblack<<std::endl;
  int k1=cdagger_index_map_[new_segment.t_start_];
  int k2=c_index_map_[new_segment.t_end_];
  
  S_tilde = operator()(k1,k2);
  weight_ratio_=1./S_tilde;
  
  // take care of sign changes due to wraparound segments
  if(new_segment.t_start_>new_segment.t_end_){
    weight_ratio_ *=-1;
  }
  
  //std::cout<<"returning removal weight ratio: "<<weight_ratio_<<" S_tilde is: "<<S_tilde<<std::endl;
  //std::cout<<S_tilde<<" weight ratio: "<<weight_ratio_<<std::endl;
  return weight_ratio_;
}

//actually remove an operator pair and change the configuration
void hybmatrix::remove_segment(const segment &new_segment, int orbital){
  //std::cout<<clblue<<"upon entering remove segment: "<<*this<<cblack<<std::endl;
  //consistency_check();
  
  //find row and column indices
  std::size_t thisrow=cdagger_index_map_[new_segment.t_start_];
  std::size_t thiscolumn=c_index_map_[new_segment.t_end_];
  //if(thisrow != thiscolumn) throw std::logic_error("row and column got screwed up!");
  std::size_t last=size()-1;
  
  //swap row and column of thisrow and the last row. Det picks up a minus sign for each interchange.
  double row_column_sign=1.;
  if(thisrow != last){
    swap_row(thisrow, last);
    row_column_sign*=-1.;
  }
  if(thiscolumn != last){
    swap_column(thiscolumn, last);
    row_column_sign*=-1.;
  }
  
  //std::cout<<"row column sign is: "<<row_column_sign<<std::endl;
  //std::cout<<"row: "<<thisrow<<" column: "<<thiscolumn<<" size: "<<size()<<std::endl;
  permutation_sign_*=row_column_sign;
  
  //perform rank one update
  /*for (std::size_t i=0; i<last; i++) {
   for (std::size_t j=0; j<last; j++) {
   operator()(i,j) -=operator()(i,last)*operator()(last,j)/operator()(last,last);
   }
   }*/
  fortran_int_t sm1=size()-1;
  if(sm1>0){ //this is exactly the content of the loops above, in dger/dgemv blas calls.
    double alpha=-1./operator()(last,last);
    fortran_int_t inc=1;
    fortran_int_t ms=memory_size();
    FORTRAN_ID(dger)(&sm1, &sm1, &alpha,&(operator()(last,0)), &inc, &(operator()(0,last)), &ms, &(operator()(0,0)), &ms);
  }
  
  
  //adjust index of operator that pointed to last, let it point to thisrow instead
  for(hyb_map_t::iterator it=c_index_map_.begin();it!=c_index_map_.end();++it){
    if(it->second==last){ it->second=thiscolumn;  break; }
  }
  for(hyb_map_t::iterator it=cdagger_index_map_.begin();it!=cdagger_index_map_.end();++it){
    if(it->second==last){ it->second=thisrow;  break; }
  }
  if(new_segment.t_start_>new_segment.t_end_){
    permutation_sign_*=-1.;
    //std::cout<<"additional permutation sign flip: t_start: "<<new_segment.t_start_<<" t_end: "<<new_segment.t_end_<<std::endl;
  }
  
  //shrink matrix size by one
  resize(size()-1);
  cdagger_index_map_.erase(new_segment.t_start_);
  c_index_map_      .erase(new_segment.t_end_);
  
  //std::cout<<clblue<<*this<<cblack<<std::endl;
  //consistency_check();
  //debug
  /*determinant_/=S_tilde;
   std::cout<<*((blas_matrix*)this)<<std::endl;
   std::cout<<clred<<"incremental determinant: "<<determinant_<<" actual determinant: "<<determinant()<<" prev determinant: "<<determinant_old_<<cblack<<std::endl;
   determinant_old_=determinant_;*/
}
std::ostream &operator<<(std::ostream &os, const hybmatrix &hyb_mat){
  os<<"hyb matrix size: "<<hyb_mat.size()<<" permutation sign: "<<hyb_mat.permutation_sign_<<std::endl;
  os<<"c map: ";
  for(hyb_map_t::const_iterator it=hyb_mat.c_index_map_.begin(); it!=hyb_mat.c_index_map_.end();++it){
    os<<"( "<<it->first<<" , "<<it->second<<" ) ";
  }
  os<<std::endl;
  os<<"cdagger map: ";
  for(hyb_map_t::const_iterator it=hyb_mat.cdagger_index_map_.begin(); it!=hyb_mat.cdagger_index_map_.end();++it){
    os<<"( "<<it->first<<" , "<<it->second<<" ) ";
  }
  std::cout<<std::endl;
  return os;
}
void hybmatrix::rebuild_hyb_matrix(int orbital, const hybfun &Delta){
  blas_matrix bup(*this);
  //build the matrix inverse:
  for(hyb_map_t::const_iterator it_start=c_index_map_.begin();it_start != c_index_map_.end();++it_start){
    for(hyb_map_t::const_iterator it_end=cdagger_index_map_.begin();it_end != cdagger_index_map_.end();++it_end){
      operator()(it_start->second,it_end->second)=Delta.interpolate(it_start->first-it_end->first, orbital);
    }
  }
  //...then invert it.
  invert();
}
void hybmatrix::rebuild_ordered_hyb_matrix(int orbital, const hybfun &Delta){
  if(size()<2) return;
  //std::cout<<"on entry rebuild orderd: full weight: "<<full_weight()<<" permutation sign: "<<permutation_sign_<<std::endl;
  //std::cout<<*this<<std::endl;
  //std::cout<<clblue<<*(blas_matrix*)this<<cblack<<std::endl;
  //order the times properly
  int k=0;
  hyb_map_t::iterator it_bup;
  for(hyb_map_t::iterator it_end=cdagger_index_map_.begin();it_end != cdagger_index_map_.end();){
    it_bup=it_end++;
    std::pair<double,int> new_entry=*it_bup;
    new_entry.second=k++;
    cdagger_index_map_.erase(it_bup);
    cdagger_index_map_.insert(new_entry);
  }
  k=0;
  for(hyb_map_t::iterator it_start=c_index_map_.begin();it_start != c_index_map_.end();){
    it_bup=it_start++;
    std::pair<double,int> new_entry=*it_bup;
    new_entry.second=((it_bup==c_index_map_.begin()) && (it_bup->first<cdagger_index_map_.begin()->first))?c_index_map_.size()-1:k++;
    c_index_map_.erase(it_bup);
    c_index_map_.insert(new_entry);
  }  //if we have an overlapping segment we need a permutation sign of -1, otherwise it is 1 in the ordered case.
  if(size()==0){
    permutation_sign_=1.;
  }else{
    if(c_index_map_.begin()->first<cdagger_index_map_.begin()->first){
      permutation_sign_=-1.;
    }else{
      permutation_sign_=1.;
    }
  }
  //then rebuild the hybridization matrix
  rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<*this<<std::endl;
  //std::cout<<clred<<*(blas_matrix*)this<<cblack<<std::endl;
  //std::cout<<"on exit rebuild orderd: full weight: "<<full_weight()<<" permutation sign: "<<permutation_sign_<<std::endl;
}
double hybmatrix::full_weight() const{
  //std::cout<<clcyan<<"det: "<<determinant()<<" ps: "<<permutation_sign_<<cblack<<std::endl;
  return determinant()*permutation_sign_;
}
void hybmatrix::measure_G(std::vector<double> &G, std::vector<double> &F, const std::map<double,double> &F_prefactor, double sign) const{
  double N_div_beta=(G.size()-1)/beta_;
  static std::vector<double> cdagger_times(size()); cdagger_times.resize(size());
  static std::vector<double> c_times(size()); c_times.resize(size());
  for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) {
    c_times[it->second] = it->first;
  }
  for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) {
    cdagger_times[it->second] = it->first;
  }
  //we measure G(tau-tau'):=-<T c(tau) c^dagger(tau')>
  for (int i = 0; i < size(); i++) {
    double f_pref=(F_prefactor.find(c_times[i]))->second;
    for (int j = 0; j < size(); j++) {
      double argument = c_times[i] - cdagger_times[j];
      double bubble_sign = sign;
      if (argument < 0) {
        bubble_sign *=-1.;
        argument += beta_;
      }
      int index = (int) (argument * N_div_beta + 0.5);
      double g = operator() (j, i) * bubble_sign;
      //NOTE:  - corresponds to -<T c(tau) c^dag(tau')>
      G[index] -= g; //changed this to-; check consistency with ALPS DMFT loop!
      F[index] -= g*f_pref;
    }
  }
}
void hybmatrix::consistency_check() const{
  for(hyb_map_t::const_iterator it1=c_index_map_.begin(); it1!= c_index_map_.end();++it1){
    for(hyb_map_t::const_iterator it2=c_index_map_.begin(); it2!= c_index_map_.end();++it2){
      if(it1->first != it2->first && it1->second==it2->second){
        std::cout<<clcyan<<"problem; inconsistent c map."<<cblack<<std::endl;
        std::cout<<*this;
        throw std::logic_error("...");
      }
    }
  }
  for(hyb_map_t::const_iterator it1=cdagger_index_map_.begin(); it1!= cdagger_index_map_.end();++it1){
    for(hyb_map_t::const_iterator it2=cdagger_index_map_.begin(); it2!= cdagger_index_map_.end();++it2){
      if(it1->first != it2->first && it1->second==it2->second){
        std::cout<<clcyan<<"problem; inconsistent c map."<<cblack<<std::endl;
        std::cout<<*this;
        throw std::logic_error("...");
      }
    }
  }
}


void hybmatrix::measure_Gl(std::vector<double> &Gl, std::vector<double> &Fl , const std::map<double,double> &F_prefactor, double sign) const{
  static std::vector<double> cdagger_times(size()); cdagger_times.resize(size());
  static std::vector<double> c_times(size()); c_times.resize(size());
  static std::vector<std::complex<double> > cdagger_exp(size()); cdagger_exp.resize(size());
  static std::vector<std::complex<double> > c_exp(size()); c_exp.resize(size());
  int N_l=Gl.size();

  //create map of creator and annihilator times
  for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) {
    c_times[it->second] = it->first;
  }
  for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) {
    cdagger_times[it->second] = it->first;
  }

  //measures the Legendre coefficients of G(tau-tau'):=-<T c(tau) c^dagger(tau')>
  for (int i = 0; i < size(); i++) {
    double f_pref=(F_prefactor.find(c_times[i]))->second;
    for (int j = 0; j < size(); j++) {
      double M_ji = operator() (j, i) * sign;
      double argument = c_times[i] - cdagger_times[j];
      double bubble_sign = sign;
      if (argument < 0) {
        bubble_sign *=-1.;
        argument += beta_;
      }
      double x=2.0*argument/beta_-1.0;
      double pl_2=1; double pl_1=x; double legendre_p;
      for(int l=0; l<N_l; l++){
        if(l==0) legendre_p=1;
        else if(l==1) legendre_p=x;
        else{ 
          legendre_p=((2*l-1)*x*pl_1-(l-1)*pl_2)/static_cast<double>(l);//l
          pl_2=pl_1; //l-2
          pl_1=legendre_p; //l-1
        }
        double gl = M_ji*legendre_p*bubble_sign;
        double fl = gl*f_pref;
        Gl[l]-=gl/beta_;
        Fl[l]-=fl/beta_;
      }
    }
  }
}


