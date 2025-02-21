/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2012 by Vito Scarola <scarola@vt.edu>
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/


#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/upper.hpp>
#include <alps/config.h>
#include <alps/alea.h>

/* Author: Vito Scarola Nov. 2012 scarola@vt.edu */

/*
 This code diagonalizes the full Hilbert space of the Lowest Landau level (LLL)
 Coulomb Hamiltonian on the sphere.   For details of the Hamiltonian see:
 G. Fano, F. Ortolani, and E. Colombo PRB 34,2670 (1986)
 J.K. Jain, Composite Fermions, Cambridge University Press 2007 Chapter 3
*/

/*********** Input Parameters ************/
int N = 5; /* Number of particles, N>=8 takes a while*/
double q = 3.0*(N-1)/2.0; /* Flux 2q=m(N-1) yields the nu=1/m  ground states approximated by the Laughlin wavefunction */
int maxstatestoprint=50; /* Maximum number of states to print out */
/*****************************************/

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;
typedef boost::numeric::ublas::vector<double> dvector_type;
typedef boost::numeric::ublas::vector<int> ivector_type;
typedef ublas::matrix<double, ublas::column_major> dmatrix_type;
typedef ublas::matrix<int, ublas::column_major> imatrix_type;

int iq = static_cast<int> (2*q);
int mmax=iq;
dvector_type factorial(4*mmax+1);

class Simulation
{
public:
    Simulation(std::string output_file)
        : filename_(output_file)
    {;}

    void run() {

        int nmax=0; /* index for the total number of basis states*/
        int Lz= (int) (N*q); /* Total z-component of momentum */
        int maxprime=56; /* use at most the first 56 primes in defining matrix*/
        int i1,i2,lin,kin,kfin,lfin,im1,im2;
        int iunequal, iunequal2,ites;
        int it11,it22,it1,it2,it,it12;
        int mi,mj,mk,ml;
        double vv;
        ivector_type ipr(maxprime);
        ivector_type m1(N+1),m2(N+1),mdt(N+1);

        /* Use N! as initial guess for matrix dimension */
        double sumlogx=0.0;
        for(int i=1; i<=N; i++) sumlogx += log(i);
        int nguess= (int) exp(sumlogx);
        

        imatrix_type lst(nguess+1,N+1);
        imatrix_type lstm(nguess+1,N+1);

        std::cout << " Computing basis states" << std::endl;
        /* Code to generate LLL basis states and compute nmax*/
        #include "states_lll.c"

        std::cout << nmax << " " << "states found" << std::endl;

        dmatrix_type HamiltonianMatrix(nmax,nmax);
        dvector_type Ltot(nmax);
        dvector_type EigenvaluesVector(nmax);
        dvector_type Gap(nmax);
        dvector_type wavevector(nmax);
        double vint[mmax+1][mmax+1][mmax+1];
        define_factorial();
        /*
        for(int l1=1; l1<=nmax; l1++){
            std::cout << std::endl;
            for(int l2=1; l2<=N; l2++)
                std::cout << lst(l1,l2) << " " ;
        }
         */


        /* Vector to define the first "maxprime" prime numbers*/
        #include "vector_of_primes.c"

        /********** Store matrix elements for lookup *******/
        for(int i=0; i<=mmax; i++) {
            for(int j=0; j<=mmax; j++) {
                for(int k=0; k<=i+j; k++) {
                    if(k <= mmax) {
                        vint[i][j][k]=vc(i,j,k,i+j-k);
                    }
                }
            }
        }

        for(int i=0; i<nmax; i++)
            for(int j=0; j<nmax; j++)
                HamiltonianMatrix(i,j)=0.;

        std::cout << " Computing Matrix Elements " << std::endl;

        /**** Generate Hamiltonian Matrix  (Eq. 23 in Fano et al. PRB 43, 2670 1986)  **********/
        for(int i1=1; i1<=nmax; i1++) {
            for(int i2=1; i2<=i1; i2++) {
                vv=0.0;
                for(int k=1; k<=N; k++) {
                    m1(k)=lst(i1,k);
                    m2(k)=lst(i2,k);
                }
                iunequal=0;
                iunequal2=0;
                for(int im1=1; im1<=N; im1++) {
                    ites=0;
                    for(int im2=1; im2<=N; im2++) {
                        if(m1(im1)==m2(im2))ites=1;
                    }
                    if(ites==0)iunequal=iunequal+1;
                }
               
                if((iunequal <= 2) && (iunequal2 <= 2)) {
                    for(int i=1; i<=N; i++) {
                        for(int j=1; j<=N; j++) {
                            kin=1;
                            kfin=N;
                            lin=1;
                            lfin=N;
                            for(int k=kin; k<=kfin; k++) {
                                for(int l=lin; l<=lfin; l++) {
                                    if((i != j) && (k != l)) {
                                        mi=m1(i);
                                        mj=m1(j);
                                        mk=m2(k);
                                        ml=m2(l);
                                       
                                        if((mi+mj) == (mk+ml)) {
                                            it11=1;
                                            it12=1;
                                            for(it1=1; it1<=N; it1++) {
                                                if(m1(it1)>maxprime || m2(it1)>maxprime) {
                                                    std::cout<<"Not enough prime numbers";
                                                    exit(EXIT_FAILURE);
                                                }
                                                if((it1 != i) && (it1 !=j)) {
                                                    it11=it11*ipr(m1(it1));
                                                }
                                                if ((it1 != k) && (it1 != l)) {
                                                    it12=it12*ipr(m2(it1));
                                                }
                                            }
                                            if(it11==it12) {
                                                double ifac=1.0;
                                                if (i>j) {
                                                    ifac=ifac*pow(-1.0,(i+j));
                                                }
                                                if(i<=j) {
                                                    ifac=ifac*pow(-1.0,(i+j+1));
                                                }
                                                if (k>l) {
                                                    ifac=ifac*pow(-1.0,(k+l));
                                                }
                                                if(k<=l) {
                                                    ifac=ifac*pow(-1.0,(k+l+1));
                                                }

                                                /*ifac accounts for fermion antisymmetry*/
                                                vv=vv+ifac*vint[mi][mj][mk];
                                              
                                            }/*end if(it11==it12) */
                                        }/*end  if((mi+mj) == (mk+ml))*/
                                    }/*end if (i != j) && (k != l)*/
                                }
                            }
                        }
                    }
                }
                // Matrix is symmetric so store upper triangular
                HamiltonianMatrix(i2-1,i1-1)=vv;
            } /* end i2 loop*/
        }/* end i1 loop*/
        std::cout << " Computing Matrix Elements Finished "<< std::endl;

        std::cout << " Diagonalizing Matrix "<< std::endl;
        bindings::lapack::syev('V', bindings::upper(HamiltonianMatrix), EigenvaluesVector);
        std::cout << " Diagonalization finished "<< std::endl;


        /**** Compute angular momentum  using < L^+L^- >  *****/
        std::cout << " Computing states with M=-1" << std::endl;
        /* Generate basis with M=-1*/
        Lz=(int) (N*q-1);
        int nmaxm=0;
        #include "states_lll_mz_minus_one.c"
        std::cout << " Finished Computing states with M=-1" << std::endl;

        std::cout << " Computing Total L of each eigenvector" << std::endl;
        dvector_type am(nmaxm+1);
        int nprint=nmax;
        if(maxstatestoprint<nmax) nprint=maxstatestoprint;

        for(int i=1; i<=nprint; i++) {
            for(int im=1; im<=nmaxm; im++) am(im)=0.0;
            double angmom=0.0;
            for(int j=1; j<=nmax; j++) {
                for(int k=1; k<=N; k++) mdt(k) = lst(j,k);
                    for(int k=1; k<=N; k++) {
                    mdt(k) = lst(j,k)-1;
                    for(int im=1; im<=nmaxm; im++) {
                        int test=0;
                        for(int k2=1; k2<=N; k2++) if(lstm(im,k2)!=mdt(k2))test = 1;
                        if(test==0) {
                            am(im) += HamiltonianMatrix(j-1,i-1)*sqrt(lst(j,k)*(iq-lst(j,k)+1) );
                            im = nmaxm+1;
                        }
                    }/*end im loop*/
                    mdt(k)=lst(j,k);
                } /* end k loop */
            } /* end j loop*/

            for(int im=1; im<=nmaxm; im++) angmom += am(im)*am(im);
            Ltot(i-1)= (-1+sqrt(1+4.0*angmom))/2.0;

        } /* end i loop */
        std::cout << " Finished Computing Total L of each eigenvector " << std::endl;
        std::cout << std::endl;


        std::cout << " Parameters" <<  std::endl;
        std::cout << " Number of Particles N=" << N << std::endl;
        std::cout << " Twice the flux iq=" << iq << std::endl;
        std::cout << " Number of basis states nmax=" << nmax << std::endl;
        std::cout << std::endl;

        double backgroundenergy = N*N/2.0/sqrt(q);
        std::cout << "L " << " Ground state energy per particle "  << std::endl;
        std::cout << (int) Ltot(0) << "      "  << (EigenvaluesVector(0)/sqrt(q) - backgroundenergy)/N  << std::endl;
        std::cout << std::endl;
        std::cout << "L " << " Excited State Energy per particle "  << std::endl;
        for(int i=2; i<=nprint; i++) std::cout << (int) Ltot(i-1) << "      "  << (EigenvaluesVector(i-1)/sqrt(q) - backgroundenergy)/N  << std::endl;
        std::cout << std::endl;

        for(int i=1; i<=nprint; i++) {
            wavevector(i-1)=Ltot(i-1)/sqrt(q);
            Gap(i-1)=(EigenvaluesVector(i-1) - EigenvaluesVector(0))/sqrt(q);
        }
        std::cout << "Wavevector " << " Energy Gap "  << std::endl;
        for(int i=2; i<=nprint; i++) std::cout << wavevector(i-1) << "      " << Gap(i-1)<< std::endl;


        //save the observables
        save(filename_);
    }
    void define_factorial() {
        factorial(0)=1;
        factorial(1)=1;
        for(int i=2; i<4*mmax+1; i++)
            factorial(i) = factorial(i-1)*i;
    }
    /* Coulumb matrix elements on the sphere*/
    double vc(int m1,int m2,int m3,int m4) {
        double answer=0.,v,v1,v2;
        double rm1= m1-q;
        double rm2= m2-q;
        double rm3= m3-q;
        double rm4= m4-q;

        for(int j=0; j<=iq; j++) {
            v1 = factorial(2*iq-2*j)*factorial(2*iq+2*j+2)/(factorial(iq-j)*factorial(iq+j+1)*factorial(iq-j)*factorial(iq+j+1));
            v2 = factorial(2*iq+2)/(factorial(iq+1)*factorial(iq+1));
            v = v1/v2/v2;
            answer += cgcoeff(q,q,rm1,rm2,j)*cgcoeff(q,q,rm3,rm4,j)*v;
        }
        return answer;
    }
    /* a function to define the Clebsh-Gordon coefficients*/
    double cgcoeff(double j1,double j2,double m1,double m2,double j3)
    {
        int mh, me, J3,v,i;
        double SUM;
        double Cl;
        const int twoje = (int)(2*j1);
        const int twojh = (int)(2*j2);
        const int low = (int)(fabs(j1-j2) +j1+j2);

        
        dvector_type Lfact(3*twojh +3*twoje +10);;

        
        J3= (int)(j3+j1+j2);
        me= (int)(m1+j1);
        mh= (int)(m2+j2);

        /************* Definition of log of factorial *************/
        Lfact(0)= 0.0;
        for(i=1; i<=3*twojh+3*twoje+10; i++) {
            Lfact(i)= Lfact(i-1)+log((double)i);
        }
        /*********** Clebsh-Gordon Coefficients *******************/
        if( (j3>=fabs(j1-j2)) && (j3<=j1+j2) && (m1>=-j1) && (m1<=j1) &&
                (m2>=-j2) && (m2<=j2) ) {

            SUM = 0.0;
            for(v=0; v<= J3 -twojh +mh; v++) {
                if( (twojh-mh+v >= 0) && (J3-twojh-v >= 0)
                        && (J3+mh+me-twoje-twojh-v >= 0) && (v+twojh-me-mh >= 0) )

                    SUM += pow(-1.0,v+me)
                           * exp (Lfact(J3-twojh+mh-v)
                                  + Lfact(twojh-mh+v)
                                  - Lfact(v)
                                  - Lfact(J3-twojh-v)
                                  - Lfact(J3+mh+me-twoje-twojh-v)
                                  - Lfact(v+twojh-me-mh) );
            }
            if((J3+mh+me-twoje-twojh >= 0) && (J3-me-mh >= 0)) {
                Cl= SUM*pow(-1.0,J3)*sqrt((2*J3-twojh-twoje+1)
                                      *exp(Lfact(J3-twoje)
                                           + Lfact(J3-twojh)
                                           + Lfact(twoje+twojh-J3)
                                           + Lfact(J3+me+mh-twoje-twojh)
                                           + Lfact(J3-me-mh)
                                           - Lfact(J3+1)
                                           - Lfact(twojh-mh)
                                           - Lfact(mh)
                                           - Lfact(twoje-me)
                                           - Lfact(me) ) );

            } else {
                Cl= 0.0;
            }
        } else {
            Cl= 0.0;
        }
        return Cl;
    }
    /* save some things to hdf5 file*/
    void save(std::string const & filename) {
        alps::hdf5::archive ar(filename,"w");
        ar << alps::make_pvp("/parameters/N", N);
        ar << alps::make_pvp("/parameters/q", q);
        /*need to include write of output vectors EigenvaluesVector() and L() to hdf5 here*/ 
    }


private:
    std::string filename_;

};


int main()
{
    std::stringstream output_name;
    output_name << "output.q_"<<q<<"N_"<<N<<".h5" ;
    Simulation sim(output_name.str());
    sim.run();
    return 0;
}



