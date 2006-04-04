//Data structure for repulsion for density density bands.

#include <alps/parameters.h>
#include "types.h"

typedef class U_matrix{
  public:
    U_matrix(const alps::Parameters &parms){
      uint n_site=parms.value_or_default("SITES", 1);
      uint n_flavor=parms.value_or_default("FLAVORS", 2);
      assert(n_site==1); //we do NOT know how to build U if nsite != 0.
      ns_=n_site;
      nf_=n_flavor;
      val_=new double[nf_*nf_];
      if(parms.defined("U_MATRIX")){
        std::ifstream u_file(std::string(parms["U_MATRIX"]).c_str());
        int i;
        int j;
        double U_ij;
        while(u_file>>i>>j>>U_ij){
          operator()(i,j)=U_ij;
        }
      } else if(n_flavor==1){
       //special case: only 1 orbital
        operator()(0,0)=(double)(parms["U"]);
      }else if (n_flavor ==2){
        //your ordinary two site problem
        assert(parms.defined("U"));
        double U=(double)(parms["U"]);
        operator()(0,0)=0; operator()(1,1)=0;
        operator()(0,1)=U; operator()(1,0)=U;
      } else{
        assert(parms.defined("U") && parms.defined("J"));     
        double U=(double)(parms["U"]);
        double J=(double)(parms["J"]);
        assemble(n_site, n_flavor, U, J);
      }
    }
    void assemble(site_t n_site, spin_t n_flavor, double U, double J){
      //this implements the U matrix for the special case of n_flavor/2 degenerate bands
      //convention according to PRB 72, 155106 (2005)
      assert(n_site==1);
      assert(n_flavor%2==0);
      double Uprime=U-2*J;
      for(spin_t flavor_i=0;flavor_i<n_flavor;flavor_i+=2){
        operator()(flavor_i  , flavor_i  )=0; //Pauli
        operator()(flavor_i+1, flavor_i+1)=0; //Pauli
        operator()(flavor_i  ,flavor_i+1)=U; //Hubbard repulsion same band
        operator()(flavor_i+1,flavor_i  )=U; //Hubbard repulsion same band
        for(spin_t flavor_j=0;flavor_j<n_flavor;flavor_j+=2){
          if(flavor_j==flavor_i) continue;
          operator()(flavor_i  ,flavor_j  )= Uprime; //Hubbard repulsion interband 
          operator()(flavor_i+1,flavor_j+1)= Uprime; //Hubbard repulsion interband 
          operator()(flavor_i  ,flavor_j+1)= Uprime-J; //Hubbard repulsion interband Hund 
          operator()(flavor_i+1,flavor_j  )= Uprime-J; //Hubbard repulsion interband Hund
          operator()(flavor_j  ,flavor_i  )= Uprime; //Hubbard repulsion interband 
          operator()(flavor_j+1,flavor_i+1)= Uprime; //Hubbard repulsion interband 
          operator()(flavor_j  ,flavor_i+1)= Uprime-J; //Hubbard repulsion interband Hund
          operator()(flavor_j+1,flavor_i  )= Uprime-J; //Hubbard repulsion interband Hund
        }
      }
      ns_=n_site;
      nf_=n_flavor;
    } 
    ~U_matrix(){
      delete[] val_;
    }
    double &operator()(spin_t flavor_i, spin_t flavor_j){
      return val_[flavor_i*nf_+flavor_j];
    }
    const double &operator() (spin_t flavor_i, spin_t flavor_j)const {
      return val_[flavor_i*nf_+flavor_j];
    }
    spin_t nf()const {return nf_;}
    spin_t ns()const {return ns_;}
  private:
    double *val_;
    site_t ns_;
    spin_t nf_;

} U_matrix;

std::ostream &operator<<(std::ostream &os, const U_matrix &U);
