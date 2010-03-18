/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id: band_structure_calculations_BHM.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- mostly done, except the ultimate use of STL library
 * 2) Replacing raw pointers -- not done yet
 *
 */


#ifndef ALPS_APPLICATIONS_BAND_STRUCTURE_CALCULATIONS_BHM_H
#define ALPS_APPLICATIONS_BAND_STRUCTURE_CALCULATIONS_BHM_H


#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>


namespace alps {
namespace applications {


extern "C" void dsteqr_(char*,int*,double*,double*,double*,int*,double*,int*); // (compz,n,d,e,z,ldz,work,info) --diagonalize a tridiagonal matrix

template<class S,class T>
class band_structure_calculations_BHM {
public:
  // _typedefs
  typedef S site_type;
  typedef S index_type;
  typedef T coordinate_type;
  typedef T parm_type;
  typedef T obs_type;


  // get functions
  inline site_type site0()  const   {  return _site0;  }

#ifdef SITE_DEPENDENT_T_U 
  inline parm_type tx0() const { return _tx_plus[_site0]; }
  inline parm_type ty0() const { return _ty_plus[_site0]; }
  inline parm_type tz0() const { return _tz_plus[_site0]; }
  inline parm_type U0()  const { return _U[_site0]; }
#else
  inline parm_type tx0() const { return _tx_raw; }
  inline parm_type ty0() const { return _ty_raw; }
  inline parm_type tz0() const { return _tz_raw; }
  inline parm_type U0()  const { return _U; }
#endif
  inline parm_type V0()  const { return _V[_site0]; }

#ifdef SITE_DEPENDENT_T_U
  inline parm_type tx_plus(site_type index)  const  { return _tx_plus[index];  }
  inline parm_type tx_minus(site_type index) const  { return _tx_minus[index]; }
  inline parm_type ty_plus(site_type index)  const  { return _ty_plus[index];  }
  inline parm_type ty_minus(site_type index) const  { return _ty_minus[index]; }
  inline parm_type tz_plus(site_type index)  const  { return _tz_plus[index];  }
  inline parm_type tz_minus(site_type index) const  { return _tz_minus[index]; }
  inline parm_type U(site_type index)        const  { return _U[index];  }
#else
  inline parm_type tx_plus(site_type index)  const  { return _tx_raw;  }
  inline parm_type tx_minus(site_type index) const  { return _tx_raw; }
  inline parm_type ty_plus(site_type index)  const  { return _ty_raw;  }
  inline parm_type ty_minus(site_type index) const  { return _ty_raw; }
  inline parm_type tz_plus(site_type index)  const  { return _tz_raw;  }
  inline parm_type tz_minus(site_type index) const  { return _tz_raw; }
  inline parm_type U(site_type index)        const  { return _U;  }
#endif
  inline parm_type V(site_type index)        const  { return _V[index];  }


private:
  // physical constants
  static const double pi   = 3.141592654;
  static const double hbar = 1.05457148;
  static const double amu  = 1.66053886;
  static const double kB   = 1.3806503;
  static const double a0   = 0.052917720859;

  double pi_sq, twice_pi;
  double _xEr2nK, _yEr2nK, _zEr2nK;  


  // variables
  int _label;     // important upon parallelization

  site_type _Nx, _Nxy, _N;
  site_type _Nx_minus, _Nx_half;
  bool _is_Nx_even;

  parm_type _xV0, _yV0, _zV0;
  parm_type _xVT, _yVT, _zVT;
  parm_type _xlambda, _ylambda, _zlambda;
  parm_type _xwaist, _ywaist, _zwaist;
  parm_type _as;
  parm_type _mass;
  parm_type g;

  site_type _site0;

  double* x_bloch;
  double* x_wannier;
  double* x;

  std::complex<obs_type>* bloch;
  obs_type* global_bloch_real;
  obs_type* global_bloch_imag;
  obs_type* wannier;

  obs_type* energy_q;

#ifdef SITE_DEPENDENT_T_U
  parm_type *_tx_raw, *_ty_raw, *_tz_raw;
  parm_type *_tx_plus, *_tx_minus, *_ty_plus, *_ty_minus, *_tz_plus, *_tz_minus;
  parm_type* _U;
#else
  parm_type _tx_raw, _ty_raw, _tz_raw;
  parm_type _U;
#endif
  parm_type* _V;


public:
  // destructors
  ~band_structure_calculations_BHM()
  {
    delete [] x;
    delete [] x_bloch, bloch, global_bloch_real, global_bloch_imag;
    delete [] x_wannier, wannier;

    delete [] energy_q;
#ifdef SITE_DEPENDENT_T_U
    delete [] _tx_raw, _ty_raw, _tz_raw;
    delete [] _tx_plus, _tx_minus, _ty_plus, _ty_minus, _tz_plus, _tz_minus;
    delete [] _U;
#endif
    delete [] _V;
  }


  // constructors
  band_structure_calculations_BHM() 
  {
    // hardcore lapack initialization
    vect = 'I'; 
    n    = MATRIX_SIZE; 
    ldz  = MATRIX_SIZE;

    // hardcore const number initialization
    pi_sq    = pi*pi;
    twice_pi = 2.*pi;

    twice_n_bloch = 2 * n_bloch;
    half_n_bloch  = n_bloch / 2;
    h_x    = 1./n_bloch;
    h_x_sq = h_x * h_x;
    h_kx   = pi/n_bloch;

    n_wannier      = N_wannier * n_bloch;
    half_n_wannier = n_wannier / 2;
  }

  band_structure_calculations_BHM(int label, site_type Nx, parm_type xV0, parm_type yV0, parm_type zV0, parm_type as, parm_type mass, parm_type xlambda, parm_type ylambda, parm_type zlambda, parm_type xwaist, parm_type ywaist, parm_type zwaist, parm_type xVT, parm_type yVT, parm_type zVT) 
  {
    // hardcore lapack initialization
    vect = 'I';
    n    = MATRIX_SIZE;
    ldz  = MATRIX_SIZE;
  
    pi_sq    = pi*pi;
    twice_pi = 2.*pi;
  
    twice_n_bloch = 2 * n_bloch;
    half_n_bloch  = n_bloch / 2;
    h_x    = 1./n_bloch;
    h_x_sq = h_x * h_x;
    h_kx   = pi/n_bloch;

    n_wannier      = N_wannier * n_bloch;
    half_n_wannier = n_wannier / 2;

    init(label,Nx,xV0,yV0,zV0,as,mass,xlambda,ylambda,zlambda,xwaist,ywaist,zwaist,xVT,yVT,zVT);  
  }
  
  // init functions
  void init(int label, site_type Nx, parm_type xV0, parm_type yV0, parm_type zV0, parm_type as, parm_type mass, parm_type xlambda, parm_type ylambda, parm_type zlambda, parm_type xwaist, parm_type ywaist, parm_type zwaist, parm_type xVT, parm_type yVT, parm_type zVT); 


  // print functions
#ifdef PRINT_BLOCH_BSCBHM
  void print_bloch(int);
#endif
#ifdef PRINT_WANNIER_BSCBHM
  void print_wannier(void);
#endif
#ifdef PRINT_BAND_PARAMETERS_BSCBHM
  void print_band_parameters(void);
#endif


  // other functions
  inline obs_type evaluate_W()  {  return std::abs(3. * (energy_q[_Nx_minus] - energy_q[0])); }
  inline site_type site_no(site_type const xindex, site_type const yindex, site_type const zindex)   { return (zindex * _Nxy) + (yindex * _Nx) + xindex; }    // d
  inline obs_type normalize_by_pi_sq(obs_type const value)    {  return (value / pi_sq);  }


private:
  // binning details for bloch functions and wannier functions
  static const int n_bloch       = 20;
  int twice_n_bloch, half_n_bloch;
  double h_x, h_x_sq, h_kx;
  static const int N_wannier      = 7;
  int n_wannier, half_n_wannier;
  int n_global_bloch, half_n_global_bloch;


  // lapack definitions
  static const int MATRIX_SIZE = 11;
  char vect;                            // starting from tridiagonal matrix : 'I'
  int n;                                // matrix size
  double d[MATRIX_SIZE];                // before: diagonal; after: eigenvalues
  double e[MATRIX_SIZE-1];              // before: off-diagonal
  double z[MATRIX_SIZE][MATRIX_SIZE];   // after: eigenvector
  double work[2*MATRIX_SIZE-2];         // work array -- must be at least 2*n-2
  int ldz;
  int info;


  // file management
  std::ofstream outFile;
  std::ifstream inFile;

  bool load_band_parameters();


  // other private member functions
  void solve_bloch_eigenstate_and_store_into_global_bloch(site_type,parm_type,parm_type);
  void setup_normalized_wannier();
  parm_type calculate_t();
  parm_type calculate_Uraw();
  inline double calculate_U(double Ux_raw, double Uy_raw, double Uz_raw)     { return (g * Ux_raw * Uy_raw * Uz_raw); }


  void perform_calculations();
};


template <class S,class T>
void band_structure_calculations_BHM<S,T>::init(int label, site_type Nx, parm_type xV0, parm_type yV0, parm_type zV0, parm_type as, parm_type mass, parm_type xlambda, parm_type ylambda, parm_type zlambda, parm_type xwaist, parm_type ywaist, parm_type zwaist, parm_type xVT, parm_type yVT, parm_type zVT)
{
  // variables
  _label = label;

  _Nx = Nx;  _Nxy = Nx*Nx;  _N = Nx*Nx*Nx;  
  _Nx_half = Nx/2;  _Nx_minus = Nx-1;
  _is_Nx_even = (((Nx%2) == 0) ? 1 : 0);    
  _site0 = site_no(_Nx_half,_Nx_half,_Nx_half);
 
  _xV0 = xV0;  _yV0 = yV0;  _zV0 = zV0;
  _xVT = xVT;  _yVT = yVT;  _zVT = zVT;   // if VT is in Er
  /*
   * _xVT = ((0.5*pi_sq*mass*amu*xVT*xVT*xlambda*xlambda)/(kB))*(1e-13);   // if VT is in Hz
   * _yVT = ((0.5*pi_sq*mass*amu*yVT*yVT*ylambda*ylambda)/(kB))*(1e-13);
   * _zVT = ((0.5*pi_sq*mass*amu*zVT*zVT*zlambda*zlambda)/(kB))*(1e-13);
   */
  _xlambda = xlambda;  _ylambda = ylambda;  _zlambda = zlambda;
  _xwaist = (2000*xwaist)/xlambda;   _ywaist = (2000*ywaist)/ylambda;  _zwaist = (2000*zwaist)/zlambda;
  _as = as;
  _mass = mass;
  g = ((32 * pi * a0 * hbar * hbar * _as) / (amu * kB * mass * xlambda * ylambda * zlambda)) * (1e9);         // in nK

  _xEr2nK = ((2 * pi_sq * hbar * hbar) / (_mass * amu * kB * _xlambda * _xlambda)) * (1e9);
  _yEr2nK = ((2 * pi_sq * hbar * hbar) / (_mass * amu * kB * _ylambda * _ylambda)) * (1e9);
  _zEr2nK = ((2 * pi_sq * hbar * hbar) / (_mass * amu * kB * _zlambda * _zlambda)) * (1e9);


  // further binning definitions
  n_global_bloch         = n_bloch * (_Nx_half + 1);
  half_n_global_bloch    = n_global_bloch / 2;


  // allocating memory for arrays  *** further cleanups necessary
  x                                = new double [_Nx];

  x_bloch                          = new double [n_bloch];
  bloch                            = new std::complex<obs_type> [n_bloch];
  global_bloch_real                = new obs_type [n_global_bloch];
  global_bloch_imag                = new obs_type [n_global_bloch];

  x_wannier                        = new double [n_wannier];
  wannier                          = new obs_type [half_n_wannier + 1];
  
  energy_q                         = new obs_type [_Nx];

#ifdef SITE_DEPENDENT_T_U
  _tx_raw = new parm_type [_N];
  _ty_raw = new parm_type [_N];
  _tz_raw = new parm_type [_N];

  _tx_plus  = new parm_type [_N];
  _tx_minus = new parm_type [_N];
  _ty_plus  = new parm_type [_N];
  _ty_minus = new parm_type [_N];
  _tz_plus  = new parm_type [_N];
  _tz_minus = new parm_type [_N];
  _U        = new parm_type [_N];  
#endif  
  _V        = new parm_type [_N];


  // initialization of arrays
  for (int index=0; index < n_bloch; ++index)     {  x_bloch[index]   = (index - half_n_bloch)   * h_x;  }
  for (int index=0; index < n_wannier; ++index)   {  x_wannier[index] = (index - half_n_wannier) * h_x;  }
  for (int index=0; index < _Nx; ++index)         {  x[index] = (_is_Nx_even ? ((index - _Nx_half) + 0.5) : (index - _Nx_half));  }


  // calculate actual band structure calculations
  perform_calculations();
}



template <class S,class T>
void band_structure_calculations_BHM<S,T>::solve_bloch_eigenstate_and_store_into_global_bloch(site_type q_no,parm_type q,parm_type V0) {

   // setting vectors
   parm_type dummy1 = 0.25 * V0 * pi_sq;
   for (int counter1=0; counter1 < MATRIX_SIZE; ++counter1) {
      parm_type dummy = ((counter1 - MATRIX_SIZE/2)*twice_pi) + q;
      d[counter1] = dummy * dummy + 2 * dummy1;
   }
   for (int counter1=0; counter1 < MATRIX_SIZE-1; ++counter1) {
      e[counter1] = -dummy1;
   }

   dsteqr_(&vect,&n,d,e,*z,&ldz,work,&info);

   // Store unnormalized Bloch function
   for (int counter1=0; counter1 < n_bloch; ++counter1) {
      std::complex<double> sum(0,0);
      for (int counter2=0; counter2 < MATRIX_SIZE; ++counter2) {
         obs_type dummy_x = (counter1 - half_n_bloch) * h_x;
         int    dummy_l = counter2 - MATRIX_SIZE/2;
         sum  +=  (z[0][counter2]) * std::polar(1.,twice_pi*dummy_l*dummy_x);
      }
      bloch[counter1] = sum;
   }

   obs_type chosen_phase_argument = std::arg(bloch[half_n_bloch]);
   std::complex<obs_type> chosen_phase = std::polar(1.,-chosen_phase_argument);
   for (int counter4=0 ; counter4 < n_bloch; ++counter4)   bloch[counter4] = chosen_phase * (bloch[counter4]);

   // Normalizing Bloch function (Trapezium rule is sufficient)    // FULLY DEBUGGED!!!
   obs_type total_sum = 0;
   for(int counter5=0; counter5 < half_n_bloch; ++counter5) 
      total_sum += norm(bloch[counter5]);
   total_sum *= 2;
   total_sum += norm(bloch[half_n_bloch]);
   total_sum *= h_x;
   obs_type sqrt_total_sum = std::sqrt(total_sum);
   for (int counter5=0; counter5 < n_bloch; ++counter5)   
      bloch[counter5] /= sqrt_total_sum;  

   // Storing rescaled energy (useful for the calculation of t)
   energy_q[q_no] = normalize_by_pi_sq(d[0]);

   // Storing into Global Bloch Array
   int start_counter   = q_no * n_bloch;
   for(int counter7=0; counter7 < n_bloch; ++counter7)   {
     int dummy_counter = start_counter + counter7;
     global_bloch_real[dummy_counter] = bloch[counter7].real();
     global_bloch_imag[dummy_counter] = bloch[counter7].imag();
   }
}


template <class S,class T>
void band_structure_calculations_BHM<S,T>::setup_normalized_wannier(void) {   

   // Determining Unnormalized Wannier Function
   if (_is_Nx_even) {
      for (int j=0; j <= half_n_wannier; ++j) {
         obs_type sum = 0;
         int bloch_counter = j % n_bloch;
         sum += ((global_bloch_real[bloch_counter] * std::cos(-pi*x_wannier[j])) - (global_bloch_imag[bloch_counter] * std::sin(-pi*x_wannier[j])));
         sum += ((global_bloch_real[_Nx_half * n_bloch + bloch_counter]));
         for (int q_no=1; q_no < _Nx_half; ++q_no) {
            parm_type q = -pi + ((2*pi / _Nx) * q_no);
            int dummy_counter = q_no * n_bloch + bloch_counter;
            sum += 2 * ((global_bloch_real[dummy_counter] * std::cos(q*x_wannier[j])) - (global_bloch_imag[dummy_counter] * std::sin(q*x_wannier[j])));
         }
         wannier[j] = sum;
      }
   }
   else {
      for (int j=0; j <= half_n_wannier; ++j) {
         obs_type sum = 0;
         int bloch_counter = j % n_bloch;
       sum += ((global_bloch_real[_Nx_half * n_bloch + bloch_counter]));
       for (int q_no=0; q_no < _Nx_half; ++q_no) {
          parm_type q = -pi + (pi/_Nx) + ((2*pi / _Nx) * q_no);
          int dummy_counter = q_no * n_bloch + bloch_counter;
          sum += 2 * ((global_bloch_real[dummy_counter] * std::cos(q*x_wannier[j])) - (global_bloch_imag[dummy_counter] * std::sin(q*x_wannier[j])));
       }
       wannier[j] = sum;
    }
 }

 // Normalizing ustd::sing Trapzium integration
 obs_type total_sum = 0;

 for(int counter=0; counter < half_n_wannier; ++counter) 
    total_sum += wannier[counter] * wannier[counter];
 total_sum += 0.5 * wannier[half_n_wannier] * wannier[half_n_wannier]; 
 total_sum *= 2 * h_x;
 obs_type sqrt_total_sum = std::sqrt(total_sum);

 for (int counter=0; counter <= half_n_wannier; counter++)
    wannier[counter] /= sqrt_total_sum;
}


template <class S,class T>
T band_structure_calculations_BHM<S,T>::calculate_t() {
 if (_is_Nx_even) {
    parm_type sum=0;
    for (site_type q_no=1; q_no < _Nx_half; ++q_no) {
       parm_type q = -pi + ((2 * pi) / _Nx) * q_no;
       sum += energy_q[q_no] * std::cos(q);
    }
    sum *= 2;
    sum += -energy_q[0];
    sum += energy_q[_Nx_half];

    return -(sum / _Nx);
 }
 else {
    parm_type sum=0;
    for (site_type q_no=0; q_no < _Nx_half; ++q_no) {
       parm_type q = -pi + (pi / _Nx) + (((2 * pi) / _Nx) * q_no);
       sum += energy_q[q_no] * std::cos(q);
    }
    sum *= 2;
    sum += energy_q[_Nx_half];

    return -(sum / _Nx);
 }
}   


template <class S,class T>
T band_structure_calculations_BHM<S,T>::calculate_Uraw(void) {
 parm_type sum=0;
 for (int index=1; index < half_n_wannier; ++index) {
    sum += std::pow(wannier[index],4);
 }
 sum *= 2;
 sum += std::pow(wannier[0],4);
 sum += std::pow(wannier[half_n_wannier],4);
 return (sum * h_x);
}


#ifdef PRINT_BLOCH_BSCBHM
template <class S,class T>
void band_structure_calculations_BHM<S,T>::print_bloch(int q_no) {
 std::ostringstream ss, ss2;
 ss  << _label;  
 ss2 << static_cast<double>(q_no);
 std::string filename1a  = "bloch_re_";
 std::string filename1b  = "bloch_im_";
 std::string filename2(ss.str());
 std::string filename3(ss2.str());
 filename1a  += filename2;  filename1a += "_"; filename1a += filename3; filename1a  += ".o";
 filename1b  += filename2;  filename1b += "_"; filename1b += filename3; filename1b  += ".o";

 int start_counter = q_no * n_bloch;

 outFile.open(filename1a.c_str(),std::ios::out);
 for (int index = 0; index < n_bloch; ++index) {
    outFile << x_bloch[index] << "    " << global_bloch_real[start_counter + index] << std::endl;
 }
 outFile.close();

 outFile.open(filename1b.c_str(),std::ios::out);
 for (int i = 0; i < n_bloch; i++) {
    outFile << x_bloch[i] << "    " << global_bloch_imag[start_counter + i] << std::endl;
 }
 outFile.close();
}
#endif


#ifdef PRINT_WANNIER_BSCBHM
template <class S,class T>
void band_structure_calculations_BHM<S,T>::print_wannier(void) {
 std::ostringstream ss;
 ss  << _label;
 std::string filename1  = "wannier_";
 std::string filename2(ss.str());
 filename1  += filename2;
 filename1 += ".o"; 

 outFile.open(filename1.c_str(),std::ios::out);
 for (int index = 0; index <= half_n_wannier; ++index) {
    outFile << x_wannier[index] << "    " << wannier[index] << std::endl;
 }
 for (int index = half_n_wannier+1; index < n_wannier; ++index) {
    outFile << x_wannier[index] << "    " << wannier[n_wannier-index] << std::endl;
 }

 outFile.close();
}
#endif


#ifdef PRINT_BAND_PARAMETERS_BSCBHM
template <class S,class T>
void band_structure_calculations_BHM<S,T>::print_band_parameters() {
std::ostringstream ss;
ss << _label;
std::string label_str(ss.str());

#define IMPLEMENT_FILE_OUTPUT1_BSCBHM(PHYSICAL_QUANTITY,FILENAME,LABEL) \
std::string LABEL = FILENAME; \
LABEL += label_str;  LABEL += ".o"; \
\
outFile.open(LABEL.c_str(),std::ios::out);  \
for (site_type index=0; index < _N; ++index)  {  outFile << PHYSICAL_QUANTITY[index] << "\n";  } \
outFile.close();
 
IMPLEMENT_FILE_OUTPUT1_BSCBHM(_tx_plus,"tx+_",filename1a)
IMPLEMENT_FILE_OUTPUT1_BSCBHM(_ty_plus,"ty+_",filename1b)
IMPLEMENT_FILE_OUTPUT1_BSCBHM(_tz_plus,"tz+_",filename1c)
IMPLEMENT_FILE_OUTPUT1_BSCBHM(_tx_minus,"tx-_",filename1d)
IMPLEMENT_FILE_OUTPUT1_BSCBHM(_ty_minus,"ty-_",filename1e)
IMPLEMENT_FILE_OUTPUT1_BSCBHM(_tz_minus,"tz-_",filename1f)
IMPLEMENT_FILE_OUTPUT1_BSCBHM(_U,"U_",filename2)
IMPLEMENT_FILE_OUTPUT1_BSCBHM(_V,"V_",filename3)

#define IMPLEMENT_FILE_OUTPUT2_BSCBHM(PHYSICAL_QUANTITY,FILENAME,LABEL) \
std::string LABEL = FILENAME; \
LABEL += label_str;  LABEL += ".o"; \
\
outFile.open(LABEL.c_str(),std::ios::out);  \
for (site_type index=_Nx_half; index < _Nx; ++index)  {  outFile << x[index] << "\t" << PHYSICAL_QUANTITY << "\n";  } \
outFile.close();

IMPLEMENT_FILE_OUTPUT2_BSCBHM(U0()/tx0(),"Utx_",filename2a)
IMPLEMENT_FILE_OUTPUT2_BSCBHM(U0()/ty0(),"Uty_",filename2b)
IMPLEMENT_FILE_OUTPUT2_BSCBHM(U0()/tz0(),"Utz_",filename2c)
}
#endif


template <class S,class T>
bool band_structure_calculations_BHM<S,T>::load_band_parameters(void) {
#ifdef SITE_DEPENDENT_T_U
std::ostringstream ss;
ss << _label;
std::string label_str(ss.str());

#define IMPLEMENT_FILE_INPUT_BSCBHM(PHYSICAL_QUANTITY,FILENAME,LABEL) \
std::string LABEL = FILENAME; \
LABEL += label_str;  LABEL += ".o"; \
\
inFile.open(LABEL.c_str(),std::ios::in);  \
if (!inFile.good())  {  return false;  }  \
for (site_type index=0; index < _N; ++index)  {  inFile >> PHYSICAL_QUANTITY[index];  } \
inFile.close();

IMPLEMENT_FILE_INPUT_BSCBHM(_tx_plus,"tx+_",filename1a)
IMPLEMENT_FILE_INPUT_BSCBHM(_ty_plus,"ty+_",filename1b)
IMPLEMENT_FILE_INPUT_BSCBHM(_tz_plus,"tz+_",filename1c)
IMPLEMENT_FILE_INPUT_BSCBHM(_tx_minus,"tx-_",filename1d)
IMPLEMENT_FILE_INPUT_BSCBHM(_ty_minus,"ty-_",filename1e)
IMPLEMENT_FILE_INPUT_BSCBHM(_tz_minus,"tz-_",filename1f)
IMPLEMENT_FILE_INPUT_BSCBHM(_U,"U_",filename2)
IMPLEMENT_FILE_INPUT_BSCBHM(_V,"V_",filename3)

return true;
#else
return false;
#endif
}


template <class S,class T>
void band_structure_calculations_BHM<S,T>::perform_calculations() {

#ifndef DEBUGMODE
  if (!load_band_parameters()) {
#endif

#ifndef DEBUGMODE
#ifdef SITE_DEPENDENT_T_U
  std::cout << "Starting to obtain site- and direction- dependent t_ij and U_i values...\n\n";
#else
  std::cout << "Starting to obtain direction- dependent only t and U values...\n\n";
#endif
#endif
#ifdef SITE_DEPENDENT_T_U
  for (site_type k_no = _Nx_half; k_no < _Nx; ++k_no) {
    site_type s      = k_no * _Nxy;
    site_type s_conj = (_Nx_minus - k_no) * _Nxy;
    coordinate_type z_coord_sq  = x[k_no] * x[k_no];

    std::cout << "Completed : " << ((static_cast<double>(k_no - _Nx_half) * 100)/(_Nx_half)) << " percent..." << std::endl;

    for (site_type j_no = _Nx_half; j_no < _Nx; ++j_no) {
      site_type q      = j_no * _Nx;
      site_type q_conj = (_Nx_minus - j_no) * _Nx;
      coordinate_type y_coord_sq = x[j_no] * x[j_no];

      for (site_type i_no= _Nx_half; i_no < _Nx; ++i_no) {
        site_type i_no_conj = _Nx_minus - i_no;
        coordinate_type x_coord_sq = x[i_no] * x[i_no];

        site_type p = s + q + i_no;
        coordinate_type r_coord_sq = x_coord_sq + y_coord_sq + z_coord_sq;

        parm_type Vx_correction_factor = exp(-2*r_coord_sq/(_xwaist*_xwaist));
        parm_type Vy_correction_factor = exp(-2*r_coord_sq/(_ywaist*_ywaist));
        parm_type Vz_correction_factor = exp(-2*r_coord_sq/(_zwaist*_zwaist));
#else
  parm_type Vx_correction_factor = 1.;
  parm_type Vy_correction_factor = 1.;
  parm_type Vz_correction_factor = 1.;
#endif
        parm_type tx, ty, tz, Ux, Uy, Uz, Uon, ep;     // Ux , Uy , Uz are raw U !!!

        if (_is_Nx_even) 
        {
          for (site_type i=0; i <= _Nx_half; ++i)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -pi+((2*pi/_Nx)*i),_xV0*Vx_correction_factor);
          tx = calculate_t() * _xEr2nK;
          setup_normalized_wannier();
          Ux = calculate_Uraw();

          for (site_type i=0; i <= _Nx_half; ++i)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -pi+((2*pi/_Nx)*i),_yV0*Vy_correction_factor);
          ty = calculate_t() * _yEr2nK;
          setup_normalized_wannier();
          Uy = calculate_Uraw();

          for (site_type i=0; i <= _Nx_half; ++i)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -pi+((2*pi/_Nx)*i),_zV0*Vz_correction_factor);
          tz = calculate_t() * _zEr2nK;
          setup_normalized_wannier();
          Uz = calculate_Uraw();

          Uon = calculate_U(Ux,Uy,Uz);
        }
        else  // _is_Nx_even is false 
        { 
          for (site_type i=0; i <= _Nx_half; ++i)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -pi+(pi/_Nx)+((2*pi/_Nx)*i),_xV0*Vx_correction_factor);
          tx = calculate_t() * _xEr2nK;
          setup_normalized_wannier();
          Ux = calculate_Uraw();

          for (site_type i=0; i <= _Nx_half; ++i)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -pi+(pi/_Nx)+((2*pi/_Nx)*i),_yV0*Vy_correction_factor);
          ty = calculate_t() * _yEr2nK;
          setup_normalized_wannier();
          Uy = calculate_Uraw();

          for (site_type i=0; i <= _Nx_half; ++i)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -pi+(pi/_Nx)+((2*pi/_Nx)*i),_zV0*Vz_correction_factor);
          tz = calculate_t() * _zEr2nK;
          setup_normalized_wannier();
          Uz = calculate_Uraw();

          Uon = calculate_U(Ux,Uy,Uz);
        }

#ifdef SITE_DEPENDENT_T_U
        ep = _xVT * x_coord_sq + _yVT * y_coord_sq + _zVT * z_coord_sq;

        if ((i_no == (_Nx-1)) || (j_no == (_Nx-1)) || (k_no == (_Nx-1)))  { ep = 1000000000000.; } 

                                            _tx_raw[p] = tx;  _ty_raw[p] = ty;  _tz_raw[p] = tz;  _U[p] = Uon;  _V[p] = ep;
        p = s + q + i_no_conj;              _tx_raw[p] = tx;  _ty_raw[p] = ty;  _tz_raw[p] = tz;  _U[p] = Uon;  _V[p] = ep;
        p = s + q_conj + i_no;              _tx_raw[p] = tx;  _ty_raw[p] = ty;  _tz_raw[p] = tz;  _U[p] = Uon;  _V[p] = ep;
        p = s + q_conj + i_no_conj;         _tx_raw[p] = tx;  _ty_raw[p] = ty;  _tz_raw[p] = tz;  _U[p] = Uon;  _V[p] = ep; 
        p = s_conj + q + i_no;              _tx_raw[p] = tx;  _ty_raw[p] = ty;  _tz_raw[p] = tz;  _U[p] = Uon;  _V[p] = ep;
        p = s_conj + q + i_no_conj;         _tx_raw[p] = tx;  _ty_raw[p] = ty;  _tz_raw[p] = tz;  _U[p] = Uon;  _V[p] = ep;
        p = s_conj + q_conj + i_no;         _tx_raw[p] = tx;  _ty_raw[p] = ty;  _tz_raw[p] = tz;  _U[p] = Uon;  _V[p] = ep;
        p = s_conj + q_conj + i_no_conj;    _tx_raw[p] = tx;  _ty_raw[p] = ty;  _tz_raw[p] = tz;  _U[p] = Uon;  _V[p] = ep;
#else
  _tx_raw = tx;  _ty_raw = ty;  _tz_raw = tz;  _U = Uon; 

  for (site_type k_no = 0; k_no < _Nx; ++k_no) {
    site_type s      = k_no * _Nxy;
    coordinate_type z_coord_sq  = x[k_no] * x[k_no];
    for (site_type j_no = 0; j_no < _Nx; ++j_no) {
      site_type q      = j_no * _Nx;
      coordinate_type y_coord_sq = x[j_no] * x[j_no];
      for (site_type i_no= 0; i_no < _Nx; ++i_no) {
        coordinate_type x_coord_sq = x[i_no] * x[i_no];
        site_type p = s + q + i_no;
        ep = _xVT * x_coord_sq + _yVT * y_coord_sq + _zVT * z_coord_sq;
        if ((i_no == 0) || (j_no == 0) || (k_no == 0) || (i_no == (_Nx-1)) || (j_no == (_Nx-1)) || (k_no == (_Nx-1)))  { ep = 100000000.; }
        _V[p]= ep;
      }
    }
  }

  std::cout << "V0 (x)  = " << _xV0 << "  ;  lambda (z) = " << _xlambda << "  ;  t (x) = " << tx << std::endl;
  std::cout << "V0 (y)  = " << _yV0 << "  ;  lambda (z) = " << _ylambda << "  ;  t (y) = " << ty << std::endl;
  std::cout << "V0 (z)  = " << _zV0 << "  ;  lambda (z) = " << _zlambda << "  ;  t (z) = " << tz << std::endl;
  std::cout << "W       = " << evaluate_W() << std::endl;
  std::cout << "U       = " << Uon << std::endl;
  std::cout << "U/t (x) = " << (U0() / tx0()) << std::endl;
  std::cout << "U/t (y) = " << (U0() / ty0()) << std::endl;
  std::cout << "U/t (z) = " << (U0() / tz0()) << std::endl;
  std::cout << std::endl;

#endif
#ifdef SITE_DEPENDENT_T_U
      }
    }
  }
#endif


#ifndef DEBUGMODE
  std::cout << "Completed totally!\n";
#endif

#ifndef DEBUGMODE
#ifdef SITE_DEPENDENT_T_U
  for (site_type k_no = 0; k_no < _Nx; ++k_no) {
    site_type s      = k_no * _Nxy;
    for (site_type j_no = 0; j_no < _Nx; ++j_no) {
      site_type q      = j_no * _Nx;
      for (site_type i_no = 0; i_no < _Nx; ++i_no) {
        site_type p = s + q + i_no;
        _tx_plus[p] = (i_no == _Nx_minus) ? 0 : 0.5*(_tx_raw[p] + _tx_raw[p+1]); 
        _ty_plus[p] = (j_no == _Nx_minus) ? 0 : 0.5*(_ty_raw[p] + _ty_raw[p+_Nx]);
        _tz_plus[p] = (k_no == _Nx_minus) ? 0 : 0.5*(_tz_raw[p] + _tz_raw[p+_Nxy]); 
        _tx_minus[p] = (i_no == 0) ? 0 : 0.5*(_tx_raw[p] + _tx_raw[p-1]);
        _ty_minus[p] = (j_no == 0) ? 0 : 0.5*(_ty_raw[p] + _ty_raw[p-_Nx]);
        _tz_minus[p] = (k_no == 0) ? 0 : 0.5*(_tz_raw[p] + _tz_raw[p-_Nxy]);
      }
    }
  }
#endif
#endif

#ifdef PRINT_BLOCH_BSCBHM
  for (site_type index=0; index < _Nx; ++index)   print_bloch(index);
#endif
#ifdef PRINT_WANNIER_BSCHM
  print_wannier();
#endif
#ifdef PRINT_BAND_PARAMETERS
  print_band_parameters(); 
#endif


#ifndef DEBUGMODE
  }
  else {
  std::cout << "Obtained BH parameters from file!" << std::endl;
#ifdef PRINT_BAND_PARAMETERS
  print_band_parameters();
#endif
  }
#endif
}


} // ending namespace applications
} // ending namespace alps


#endif
