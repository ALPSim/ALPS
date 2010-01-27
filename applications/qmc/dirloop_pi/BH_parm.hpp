/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2007-2010 by Ping Nang (Tama) Ma <pingnang@phys.ethz.ch>
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

// Version 1.0
// No precondition, postcondition, invariance concept has been used here. 
//
// Band Structure Calculation for single-band boson Hubbard model
// Option 1: Site- and Direction-dependent nearest neighbour hopping t_ij and onsite interaction strength U_i
// Option 2: Direction-dependent nearest neighbour hopping t and onsite interaction strength U  
// 
// The parameters will be printed to file.
// 
// Written by:  
// Ping Nang MA
// HIT 31.3 Institut. f. Theoretische Physik
// Wolfgang-Pauli-Strasse 27,
// ETH Honggerberg
// 8093 Zurich
//
// For any problem about this code, please direct an email to <pingnang@itp.phys.ethz.ch>
//


//#define DEBUGMODE
//#define SITE_DEPENDENT_T_U                       // please comment it if you are conidering site independent t and U


#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>


// dsteqr : computes all eigenvalues and eigenvectors of a symmetric or Hermitian matrix reduced to tridiagonal form (QR-factorization)
// argument : (compz,n,d,e,z,ldz,work,info)
extern "C" void dsteqr_(char*,int*,double*,double*,double*,int*,double*,int*);

//std::ofstream out_potential("out_V0",std::ios::app);
std::ofstream out_bloch;
std::ofstream out_wannier;
std::ofstream out_BH_parm;
std::ifstream in_BH_parm;

#define PI          3.141592654
#define hbar        1.05457148
#define k_B         1.3806503
#define bohr_radius 0.052917720859
#define amu         1.66053886

static const int MATRIX_SIZE = 11;               // *** modify this, odd is a must!!!
static const int n_bloch = 20;
static const int N_wannier = 7;

const int MATRIX_SIZE_minus = MATRIX_SIZE - 1;
const int half_MATRIX_SIZE  = MATRIX_SIZE / 2;
const int twice_MATRIX_SIZE = 2 * MATRIX_SIZE;

// LAPACK parameter
char vect = 'I';                    // 'vect' = 'compz'
int n = MATRIX_SIZE;
double d[MATRIX_SIZE];
double e[MATRIX_SIZE_minus];
double z[MATRIX_SIZE][MATRIX_SIZE];
double work[twice_MATRIX_SIZE];
int ldz = MATRIX_SIZE;                 // 'ldz' = 'ldq'
int info;


class BH_PARM {
public:
  BH_PARM(int label_, int N_x_, double V0_0_x_, double V0_0_y_, double V0_0_z_, double a_s_, double mass_, double lambda_x_, double lambda_y_, double lambda_z_, double waist_x_, double waist_y_, double waist_z_, double VT_x_, double VT_y_, double VT_z_) : 
    label(label_),
    N_x(N_x_),
    V0_0_x(V0_0_x_),
    V0_0_y(V0_0_y_),
    V0_0_z(V0_0_z_),
    a_s(a_s_),
    mass(mass_),
    lambda_x(lambda_x_),
    lambda_y(lambda_y_),
    lambda_z(lambda_z_),
    waist_x((2000*waist_x_)/lambda_x_),
    waist_y((2000*waist_y_)/lambda_y_),
    waist_z((2000*waist_z_)/lambda_z_),
    VT_x(VT_x_),
    VT_y(VT_y_),
    VT_z(VT_z_)
  {
    N_x_even = (((N_x%2) == 0) ? 1 : 0);     // N_x is even -> 1 , or , N_x is odd -> 0

#ifdef DEBUGMODE
    std::cout << "Nx = " << N_x << " ";
    std::cout << (N_x_even ? "(even)" : "(odd)") << std::endl;
#endif

    N_x_sq   = N_x * N_x;
    N_x_cube = N_x * N_x * N_x;
    half_N_x = N_x / 2;
    N_x_minus = N_x - 1;

    half_n_bloch  = n_bloch / 2;
    twice_n_bloch = 2 * n_bloch;

    n_wannier          = N_wannier * n_bloch;
    half_n_wannier     = n_wannier / 2;

    n_global_bloch         = n_bloch * (half_N_x + 1);
    half_n_global_bloch    = n_global_bloch / 2;

    h_x  = 1. / n_bloch;
    h_kx = PI / n_bloch;

    h_x_sq = h_x * h_x;

    g = ((32 * PI * bohr_radius * hbar * hbar * a_s) / (amu * k_B * mass * lambda_x * lambda_y * lambda_z)) * (1e9);         // in nK

    Er_to_nK_x = ((2 * PI * PI * hbar * hbar) / (mass * amu * k_B * lambda_x * lambda_x)) * (1e9);
    Er_to_nK_y = ((2 * PI * PI * hbar * hbar) / (mass * amu * k_B * lambda_y * lambda_y)) * (1e9);
    Er_to_nK_z = ((2 * PI * PI * hbar * hbar) / (mass * amu * k_B * lambda_z * lambda_z)) * (1e9);

    // // this is necessary if VT is in Hz
    // VT_x = ((0.5 * PI * PI * mass * amu * VT_x * VT_x * lambda_x * lambda_x)/(k_B)) * (1e-13);
    // VT_y = ((0.5 * PI * PI * mass * amu * VT_y * VT_y * lambda_y * lambda_y)/(k_B)) * (1e-13);
    // VT_z = ((0.5 * PI * PI * mass * amu * VT_z * VT_z * lambda_z * lambda_z)/(k_B)) * (1e-13);  

    PI_sq    = PI * PI;
    twice_PI = 2 * PI;

    x_bloch                          = new double [n_bloch];
    x_wannier                        = new double [n_wannier];
    x                                = new double [N_x];
    
    bloch                            = new std::complex<double> [n_bloch];
    energy_q                         = new double [N_x];
    global_bloch_real                = new double [n_global_bloch];
    global_bloch_imag                = new double [n_global_bloch];
    wannier                          = new double [half_n_wannier + 1];

    t_x_raw = new double [N_x_cube];
    t_y_raw = new double [N_x_cube];
    t_z_raw = new double [N_x_cube];

    t_x_plus_  = new double [N_x_cube];
    t_x_minus_ = new double [N_x_cube];
    t_y_plus_  = new double [N_x_cube];
    t_y_minus_ = new double [N_x_cube];
    t_z_plus_  = new double [N_x_cube];
    t_z_minus_ = new double [N_x_cube];
    U_         = new double [N_x_cube];    
    epsilon_   = new double [N_x_cube];

    initialize();
    perform_simulation();
  }
/*
  ~BH_PARM() {
     delete [] x_bloch;
     delete [] x_wannier;
     delete [] x;

     delete [] bloch;
     delete [] global_bloch_real;
     delete [] global_bloch_imag;
     delete [] wannier;

     delete [] energy_q;

     delete [] t_x_raw;
     delete [] t_y_raw;
     delete [] t_z_raw;

     delete [] t_x_plus_;
     delete [] t_x_minus_;
     delete [] t_y_plus_;
     delete [] t_y_minus_;
     delete [] t_z_plus_;
     delete [] t_z_minus_;
     delete [] U_;

     delete [] epsilon_;
  }
*/
  inline void t_x_plus(double* s)    { for (int i=0; i < N_x_cube; i++)    s[i] = t_x_plus_[i];  }
  inline void t_x_minus(double* s)   { for (int i=0; i < N_x_cube; i++)    s[i] = t_x_minus_[i]; }
  inline void t_y_plus(double* s)    { for (int i=0; i < N_x_cube; i++)    s[i] = t_y_plus_[i];  }
  inline void t_y_minus(double* s)   { for (int i=0; i < N_x_cube; i++)    s[i] = t_y_minus_[i]; }
  inline void t_z_plus(double* s)    { for (int i=0; i < N_x_cube; i++)    s[i] = t_z_plus_[i];  }
  inline void t_z_minus(double* s)   { for (int i=0; i < N_x_cube; i++)    s[i] = t_z_minus_[i]; }
  inline void U(double* s)           { for (int i=0; i < N_x_cube; i++)    s[i] = U_[i];         }
  inline void epsilon(double* s)     { for (int i=0; i < N_x_cube; i++)    s[i] = epsilon_[i];   }

  uint32_t centre_site;

private:
  inline uint32_t site_no(int x_no_in, int y_no_in, int z_no_in)                { return (z_no_in * N_x_sq) + (y_no_in * N_x) + x_no_in; }

  void initialize(void);
  void solve_bloch_eigenstate_and_store_into_global_bloch(int,double,double);
  void setup_normalized_wannier(void);

  void perform_simulation(void);

  inline double obtain_energy_eigenvalue(double E_in)    {  return (E_in / PI_sq);  }
  
  double calculate_t(void);
  double calculate_U_raw(void);
  inline double calculate_U(double U_x_raw, double U_y_raw, double U_z_raw)     { return (g * U_x_raw * U_y_raw * U_z_raw); }

  void print_bloch(int);
  void print_wannier(void);
  void print_BH_parm(void);

  bool obtain_BH_parm_from_file(void);

  int twice_n_bloch;
  int half_n_bloch;

  int n_wannier;
  int half_n_wannier;

  int n_global_bloch;
  int half_n_global_bloch;

  double h_x;
  double h_kx;

  double h_x_sq;

  double Er_to_nK_x;
  double Er_to_nK_y;
  double Er_to_nK_z;

  double PI_sq;
  double twice_PI;

  int label;

  bool N_x_even;

  int N_x;
  int N_x_sq;
  int N_x_minus;
  int half_N_x;
  uint32_t N_x_cube;
  double V0_0_x;
  double V0_0_y;
  double V0_0_z;
  double a_s;
  double mass;
  double lambda_x;
  double lambda_y;
  double lambda_z;
  double waist_x;
  double waist_y;
  double waist_z;
  double VT_x;
  double VT_y;
  double VT_z;

  double g;

  double* x_bloch;
  double* x_wannier;
  double* x;

  std::complex<double>* bloch;
  double* global_bloch_real;
  double* global_bloch_imag;
  double* wannier;

  double* energy_q;

  double* t_x_raw;
  double* t_y_raw;
  double* t_z_raw;

  double* t_x_plus_;
  double* t_x_minus_;
  double* t_y_plus_;
  double* t_y_minus_;
  double* t_z_plus_;
  double* t_z_minus_;
  double* U_;
  double* epsilon_;
};


void BH_PARM::initialize(void) {
   for (int i=0; i < n_bloch; i++)            x_bloch[i]   = (i - half_n_bloch) * h_x;
   for (int i=0; i < n_wannier; i++)          x_wannier[i] = (i - half_n_wannier) * h_x;  

   if (N_x_even) {
      for (int i=0; i < N_x; i++)             x[i]         = (i - half_N_x) + 0.5;
   }
   else {
      for (int i=0; i < N_x; i++)             x[i]         = (i - half_N_x);
   }
}


void BH_PARM::solve_bloch_eigenstate_and_store_into_global_bloch(int q_no,double q,double V0_0_in) {

   // setting vectors
   double dummy1 = 0.25 * V0_0_in * PI_sq;
   for (int counter1=0; counter1 < MATRIX_SIZE; counter1++) {
      double dummy = ((counter1 - half_MATRIX_SIZE)*twice_PI) + q;
      d[counter1] = dummy * dummy + 2 * dummy1;
   }
   for (int counter1=0; counter1 < MATRIX_SIZE_minus; counter1++) {
      e[counter1] = -dummy1;
   }

   dsteqr_(&vect,&n,d,e,*z,&ldz,work,&info);

   // Store unnormalized Bloch function
   for (int counter1=0; counter1 < n_bloch; counter1++) {
      std::complex<double> sum(0,0);
      for (int counter2=0; counter2 < MATRIX_SIZE; counter2++) {
         double dummy_x = (counter1 - half_n_bloch) * h_x;
         int    dummy_l = counter2 - half_MATRIX_SIZE;
         sum  +=  (z[0][counter2]) * std::polar(1.,twice_PI*dummy_l*dummy_x);
      }
      bloch[counter1] = sum;
   }

   double chosen_phase_argument = std::arg(bloch[half_n_bloch]);
   std::complex<double> chosen_phase = std::polar(1.,-chosen_phase_argument);
   for (int counter4=0 ; counter4 < n_bloch; counter4++)   bloch[counter4] = chosen_phase * (bloch[counter4]);

   // Normalizing Bloch function (Trapezium rule is sufficient)    // FULLY DEBUGGED!!!
   double total_sum = 0;
   for(int counter5=0; counter5 < half_n_bloch; counter5++) 
      total_sum += norm(bloch[counter5]);
   total_sum *= 2;
   total_sum += norm(bloch[half_n_bloch]);
   total_sum *= h_x;
   double sqrt_total_sum = sqrt(total_sum);
   for (int counter5=0; counter5 < n_bloch; counter5++)   
      bloch[counter5] /= sqrt_total_sum;  

   // Storing rescaled energy (useful for the calculation of t)
   energy_q[q_no] = obtain_energy_eigenvalue(d[0]);

   // Storing into Global Bloch Array
   int start_counter   = q_no * n_bloch;
   for(int counter7=0; counter7 < n_bloch; counter7++)   {
     int dummy_counter = start_counter + counter7;
     global_bloch_real[dummy_counter] = bloch[counter7].real();
     global_bloch_imag[dummy_counter] = bloch[counter7].imag();
   }
}


void BH_PARM::setup_normalized_wannier(void) {

   // Determining Unnormalized Wannier Function
   if (N_x_even) {
      for (int j=0; j <= half_n_wannier; j++) {
         double sum = 0;
         int bloch_counter = j % n_bloch;
         sum += ((global_bloch_real[bloch_counter] * cos(-PI*x_wannier[j])) - (global_bloch_imag[bloch_counter] * sin(-PI*x_wannier[j])));
         sum += ((global_bloch_real[half_N_x * n_bloch + bloch_counter]));
         for (int q_no=1; q_no < half_N_x; q_no++) {
            double q = -PI + ((2*PI / N_x) * q_no);
            int dummy_counter = q_no * n_bloch + bloch_counter;
            sum += 2 * ((global_bloch_real[dummy_counter] * cos(q*x_wannier[j])) - (global_bloch_imag[dummy_counter] * sin(q*x_wannier[j])));
         }
         wannier[j] = sum;
      }
   }
   else {
      for (int j=0; j <= half_n_wannier; j++) {
         double sum = 0;
         int bloch_counter = j % n_bloch;
         sum += ((global_bloch_real[half_N_x * n_bloch + bloch_counter]));
         for (int q_no=0; q_no < half_N_x; q_no++) {
            double q = -PI + (PI/N_x) + ((2*PI / N_x) * q_no);
            int dummy_counter = q_no * n_bloch + bloch_counter;
            sum += 2 * ((global_bloch_real[dummy_counter] * cos(q*x_wannier[j])) - (global_bloch_imag[dummy_counter] * sin(q*x_wannier[j])));
         }
         wannier[j] = sum;
      }
   }

   // Normalizing using Simpson's Rule
   double total_sum = 0;

   for(int counter=0; counter < half_n_wannier; counter++) 
      total_sum += wannier[counter] * wannier[counter];
   total_sum += 0.5 * wannier[half_n_wannier] * wannier[half_n_wannier]; 
   total_sum *= 2 * h_x;
   double sqrt_total_sum = sqrt(total_sum);

   for (int counter=0; counter <= half_n_wannier; counter++)
      wannier[counter] /= sqrt_total_sum;
}


double BH_PARM::calculate_t(void) {
   if (N_x_even) {
      double sum=0;
      for (int q_no=1; q_no < half_N_x; q_no++) {
         double q = -PI + ((2 * PI) / N_x) * q_no;
         sum += energy_q[q_no] * cos(q);
      }
      sum *= 2;
      sum += -energy_q[0];
      sum += energy_q[half_N_x];

      return -(sum / N_x);
   }
   else {
      double sum=0;
      for (int q_no=0; q_no < half_N_x; q_no++) {
         double q = -PI + (PI / N_x) + (((2 * PI) / N_x) * q_no);
         sum += energy_q[q_no] * cos(q);
      }
      sum *= 2;
      sum += energy_q[half_N_x];

      return -(sum / N_x);
   }
}   


double BH_PARM::calculate_U_raw(void) {
   double sum=0;
   for (int i=1; i < half_n_wannier; i++) {
      sum += pow(wannier[i],4);
   }
   sum *= 2;
   sum += pow(wannier[0],4);
   sum += pow(wannier[half_n_wannier],4);
   return (sum * h_x);
}


void BH_PARM::print_bloch(int q_no) {
   std::ostringstream ss, ss2;
   ss  << label;  
   ss2 << static_cast<double>(q_no);
   std::string filename1a  = "bloch_re_";
   std::string filename1b  = "bloch_im_";
   std::string filename2(ss.str());
   std::string filename3(ss2.str());
   filename1a  += filename2;  filename1a += "_"; filename1a += filename3; filename1a  += ".o";
   filename1b  += filename2;  filename1b += "_"; filename1b += filename3; filename1b  += ".o";

   int start_counter = q_no * n_bloch;

   out_bloch.open(filename1a.c_str(),std::ios::out);
   for (int i = 0; i < n_bloch; i++) {
      out_bloch << x_bloch[i] << "    " << global_bloch_real[start_counter + i] << std::endl;
   }
   out_bloch.close();

   out_bloch.open(filename1b.c_str(),std::ios::out);
   for (int i = 0; i < n_bloch; i++) {
      out_bloch << x_bloch[i] << "    " << global_bloch_imag[start_counter + i] << std::endl;
   }
   out_bloch.close();
}


void BH_PARM::print_wannier(void) {
   std::ostringstream ss;
   ss  << label;
   std::string filename1  = "wannier_";
   std::string filename2(ss.str());
   filename1  += filename2;
   filename1 += ".o"; 

   out_wannier.open(filename1.c_str(),std::ios::out);
   for (int i = 0; i <= half_n_wannier; i++) {
      out_wannier << x_wannier[i] << "    " << wannier[i] << std::endl;
   }
   for (int i = half_n_wannier+1; i < n_wannier; i++) {
      out_wannier << x_wannier[i] << "    " << wannier[n_wannier-i] << std::endl;
   }

   out_wannier.close();
}


void BH_PARM::print_BH_parm(void) {
   std::ostringstream ss;
   ss << label;

   std::string filename5a  = "t_dx_plus_";
   std::string filename5b  = "t_dy_plus_";
   std::string filename5c  = "t_dz_plus_";
   std::string filename5d  = "t_dx_minus_";
   std::string filename5e  = "t_dy_minus_";
   std::string filename5f  = "t_dz_minus_";

   std::string filename8a  = "U_t_x_";
   std::string filename8b  = "U_t_y_";
   std::string filename8c  = "U_t_z_";

   //std::string filename   = "wannier_k0_";
   std::string filename2  = "U_";
   std::string filename3  = "epsilon_";

   std::string filename6(ss.str());

   filename5a += filename6; filename5a  += ".o";
   filename5b += filename6; filename5b  += ".o";
   filename5c += filename6; filename5c  += ".o";
   filename5d += filename6; filename5d  += ".o";
   filename5e += filename6; filename5e  += ".o";
   filename5f += filename6; filename5f  += ".o";
   filename2 += filename6; filename2 += ".o";
   filename3 += filename6; filename3 += ".o";
   //filename += filename6; filename  += ".o";

   filename8a += filename6; filename8a  += ".o";
   filename8b += filename6; filename8b  += ".o";
   filename8c += filename6; filename8c  += ".o";

   // print U
   out_BH_parm.open(filename2.c_str(),std::ios::out);
   for (int i = 0; i < N_x_cube; i++)     out_BH_parm << U_[i] << std::endl;
   out_BH_parm.close();

   // print epsilon
   out_BH_parm.open(filename3.c_str(),std::ios::out);
   for (uint32_t i = 0; i < N_x_cube; i++)     out_BH_parm << epsilon_[i] << std::endl;
   out_BH_parm.close();

   // print t
   out_BH_parm.open(filename5a.c_str(),std::ios::out);
   for (uint32_t i = 0; i < N_x_cube; i++)     out_BH_parm << t_x_plus_[i] << std::endl;
   out_BH_parm.close();

   out_BH_parm.open(filename5b.c_str(),std::ios::out);
   for (uint32_t i = 0; i < N_x_cube; i++)     out_BH_parm << t_y_plus_[i] << std::endl;
   out_BH_parm.close();

   out_BH_parm.open(filename5c.c_str(),std::ios::out);
   for (uint32_t i = 0; i < N_x_cube; i++)     out_BH_parm << t_z_plus_[i] << std::endl;
   out_BH_parm.close();

   out_BH_parm.open(filename5d.c_str(),std::ios::out);
   for (uint32_t i = 0; i < N_x_cube; i++)     out_BH_parm << t_x_minus_[i] << std::endl;
   out_BH_parm.close();

   out_BH_parm.open(filename5e.c_str(),std::ios::out);
   for (uint32_t i = 0; i < N_x_cube; i++)     out_BH_parm << t_y_minus_[i] << std::endl;
   out_BH_parm.close();

   out_BH_parm.open(filename5f.c_str(),std::ios::out);
   for (uint32_t i = 0; i < N_x_cube; i++)     out_BH_parm << t_z_minus_[i] << std::endl;
   out_BH_parm.close();

   out_BH_parm.open(filename8a.c_str(),std::ios::out);
   for (int i=half_N_x; i < N_x; i++) {
      uint32_t current_site = half_N_x * N_x_sq + half_N_x * N_x + i;
      out_BH_parm << x[i] << "    " << U_[current_site] / t_x_minus_[current_site] << std::endl;
   }
   out_BH_parm.close();

   out_BH_parm.open(filename8b.c_str(),std::ios::out);
   for (int i=half_N_x; i < N_x; i++) {
      uint32_t current_site = half_N_x * N_x_sq + i * N_x + half_N_x;
      out_BH_parm << x[i] << "    " << U_[current_site] / t_y_minus_[current_site] << std::endl;
   }
   out_BH_parm.close();

   out_BH_parm.open(filename8c.c_str(),std::ios::out);
   for (int i=half_N_x; i < N_x; i++) {
      uint32_t current_site = i * N_x_sq + half_N_x * N_x + half_N_x;
      out_BH_parm << x[i] << "    " << U_[current_site] / t_z_minus_[current_site] << std::endl;
   }
   out_BH_parm.close();
}


bool BH_PARM::obtain_BH_parm_from_file(void) {
   std::ostringstream ss;
   ss << label;

   std::string filename5a  = "t_dx_plus_";
   std::string filename5b  = "t_dy_plus_";
   std::string filename5c  = "t_dz_plus_";
   std::string filename5d  = "t_dx_minus_";
   std::string filename5e  = "t_dy_minus_";
   std::string filename5f  = "t_dz_minus_";

   //std::string filename   = "wannier_k0_";
   std::string filename2  = "U_";
   std::string filename3  = "epsilon_";

   std::string filename6(ss.str());

   filename5a += filename6; filename5a  += ".o";
   filename5b += filename6; filename5b  += ".o";
   filename5c += filename6; filename5c  += ".o";
   filename5d += filename6; filename5d  += ".o";
   filename5e += filename6; filename5e  += ".o";
   filename5f += filename6; filename5f  += ".o";
   filename2 += filename6; filename2 += ".o";
   filename3 += filename6; filename3 += ".o";
   //filename += filename6; filename  += ".o";

   
   // default: if "U_<label>.o" exists, then read from file...
   in_BH_parm.open(filename2.c_str(),std::ios::in);   
   if (!in_BH_parm.good()) {
      return (false);
   }
   for (int i=0; i < N_x_cube; i++)       in_BH_parm >> U_[i];
   in_BH_parm.close();

   in_BH_parm.open(filename3.c_str(),std::ios::in);
   for (int i=0; i < N_x_cube; i++)       in_BH_parm >> epsilon_[i];
   in_BH_parm.close();

   in_BH_parm.open(filename5a.c_str(),std::ios::in);
   for (int i=0; i < N_x_cube; i++)       in_BH_parm >> t_x_plus_[i];
   in_BH_parm.close();

   in_BH_parm.open(filename5b.c_str(),std::ios::in);
   for (int i=0; i < N_x_cube; i++)       in_BH_parm >> t_y_plus_[i];
   in_BH_parm.close();

   in_BH_parm.open(filename5c.c_str(),std::ios::in);
   for (int i=0; i < N_x_cube; i++)       in_BH_parm >> t_z_plus_[i];
   in_BH_parm.close();

   in_BH_parm.open(filename5d.c_str(),std::ios::in);
   for (int i=0; i < N_x_cube; i++)       in_BH_parm >> t_x_minus_[i];
   in_BH_parm.close();

   in_BH_parm.open(filename5e.c_str(),std::ios::in);
   for (int i=0; i < N_x_cube; i++)       in_BH_parm >> t_y_minus_[i];
   in_BH_parm.close();

   in_BH_parm.open(filename5f.c_str(),std::ios::in);
   for (int i=0; i < N_x_cube; i++)       in_BH_parm >> t_z_minus_[i];
   in_BH_parm.close();

   return (true);
}


void BH_PARM::perform_simulation(void) {

#ifndef DEBUGMODE
   if (!obtain_BH_parm_from_file()) {
#endif

#ifndef DEBUGMODE
#ifdef SITE_DEPENDENT_T_U
   std::cout << "Starting to obtain site- and direction- dependent t and U values..." << std::endl << std::endl;
#else
   std::cout << "Starting to obtain direction-independent t and U values..." << std::endl << std::endl;
#endif
#endif
   if (N_x_even) {
#ifdef SITE_DEPENDENT_T_U
      for (int k_no = half_N_x; k_no < N_x; k_no++) {
         uint32_t s      = k_no * N_x_sq;
         uint32_t s_conj = (N_x_minus - k_no) * N_x_sq;

         double z_coord_sq  = x[k_no] * x[k_no];

         std::cout << "Completed : " << ((static_cast<double>(k_no - half_N_x) * 100)/(half_N_x)) << " percent..." << std::endl;

         for (int j_no = half_N_x; j_no < N_x; j_no++) {
            uint32_t q      = j_no * N_x;
            uint32_t q_conj = (N_x_minus - j_no) * N_x;

            double y_coord_sq = x[j_no] * x[j_no];

            for (int i_no= half_N_x; i_no < N_x; i_no++) {
               int i_no_conj = N_x_minus - i_no;

               double x_coord_sq = x[i_no] * x[i_no];

               //cout << x_coord << "     " << y_coord << "    " << z_coord << endl;
               //cout << waist_x << "     " << waist_y << "    " << waist_z << endl;

               uint32_t p = s + q + i_no;

               double r_coord_sq = x_coord_sq + y_coord_sq + z_coord_sq;

               double Vx_correction_factor = exp(-2*r_coord_sq/(waist_x*waist_x));
               double Vy_correction_factor = exp(-2*r_coord_sq/(waist_y*waist_y));
               double Vz_correction_factor = exp(-2*r_coord_sq/(waist_z*waist_z));

#else
               double Vx_correction_factor = 1;
               double Vy_correction_factor = 1;
               double Vz_correction_factor = 1;

#endif

               double t_x, t_y, t_z, U_x, U_y, U_z, U_on, ep;     // U_x , U_y , U_z are raw U !!!

               for (int i=0; i <= half_N_x; i++)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -PI+((2*PI/N_x)*i),V0_0_x*Vx_correction_factor);
               t_x = calculate_t() * Er_to_nK_x;
               setup_normalized_wannier();
               U_x = calculate_U_raw();

               for (int i=0; i <= half_N_x; i++)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -PI+((2*PI/N_x)*i),V0_0_y*Vy_correction_factor);
               t_y = calculate_t() * Er_to_nK_y;
               setup_normalized_wannier();
               U_y = calculate_U_raw();

               for (int i=0; i <= half_N_x; i++)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -PI+((2*PI/N_x)*i),V0_0_z*Vz_correction_factor);
               t_z = calculate_t() * Er_to_nK_z;
               setup_normalized_wannier();
               U_z = calculate_U_raw();

/*
               t_z = t_y;
               U_z = U_y;
*/
           
               U_on = calculate_U(U_x,U_y,U_z);

#ifdef SITE_DEPENDENT_T_U
               ep = VT_x * x_coord_sq + VT_y * y_coord_sq + VT_z * z_coord_sq;

               if ((i_no == (N_x-1)) || (j_no == (N_x-1)) || (k_no == (N_x-1)))  { ep = 100000000.; } 

               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s + q + i_no_conj;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s + q_conj + i_no;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s + q_conj + i_no_conj;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s_conj + q + i_no;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s_conj + q + i_no_conj;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s_conj + q_conj + i_no;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s_conj + q_conj + i_no_conj;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

#else
               std::cout << "V0 (x) = " << V0_0_x << "  ;  lambda_x = " << lambda_x << "  ;  t (x) = " << t_x << std::endl;
               std::cout << "V0 (y) = " << V0_0_y << "  ;  lambda_y = " << lambda_y << "  ;  t (y) = " << t_y << std::endl;
               std::cout << "V0 (z) = " << V0_0_z << "  ;  lambda_z = " << lambda_z << "  ;  t (z) = " << t_z << std::endl;
               std::cout << "U = " << U_on << std::endl;
               std::cout << "U/t (x) = " << (U_on / t_x) << std::endl;
               std::cout << "U/t (y) = " << (U_on / t_y) << std::endl;
               std::cout << "U/t (z) = " << (U_on / t_z) << std::endl;
               std::cout << std::endl;

               for (uint32_t p=0; p < N_x_cube; p++) {
                  t_x_raw[p] = t_x;
                  t_y_raw[p] = t_y;
                  t_z_raw[p] = t_z;
                  U_[p]      = U_on;
               }

               for (int k_no = 0; k_no < N_x; k_no++) {
                  uint32_t s      = k_no * N_x_sq;
                  double z_coord_sq  = x[k_no] * x[k_no];
                  for (int j_no = 0; j_no < N_x; j_no++) {
                     uint32_t q      = j_no * N_x;
                     double y_coord_sq = x[j_no] * x[j_no];
                     for (int i_no= 0; i_no < N_x; i_no++) {
                        double x_coord_sq = x[i_no] * x[i_no];
                        uint32_t p = s + q + i_no;
                        ep = VT_x * x_coord_sq + VT_y * y_coord_sq + VT_z * z_coord_sq;
                        if ((i_no == 0) || (j_no == 0) || (k_no == 0) || (i_no == (N_x-1)) || (j_no == (N_x-1)) || (k_no == (N_x-1)))  { ep = 100000000.; }
                        epsilon_[p]= ep;
                     }
                  }
               }

#endif
#ifdef SITE_DEPENDENT_T_U
            }
         }
      }
#endif
   }

   else {  // N_x == odd
#ifdef SITE_DEPENDENT_T_U
      for (int k_no = half_N_x; k_no < N_x; k_no++) {
         uint32_t s      = k_no * N_x_sq;
         uint32_t s_conj = (N_x_minus - k_no) * N_x_sq;

         double z_coord_sq  = x[k_no] * x[k_no];

         std::cout << "Completed : " << ((static_cast<double>(k_no - half_N_x) * 100)/(half_N_x)) << " percent..." << std::endl;

         for (int j_no = half_N_x; j_no < N_x; j_no++) {
            uint32_t q      = j_no * N_x;
            uint32_t q_conj = (N_x_minus - j_no) * N_x;

            double y_coord_sq = x[j_no] * x[j_no];

            for (int i_no= half_N_x; i_no < N_x; i_no++) {
               int i_no_conj = N_x_minus - i_no;

               double x_coord_sq = x[i_no] * x[i_no];

               //std::cout << x_coord << "     " << y_coord << "    " << z_coord << endl;
               //std::cout << waist_x << "     " << waist_y << "    " << waist_z << endl;

               double r_coord_sq = x_coord_sq + y_coord_sq + z_coord_sq;

               uint32_t p = s + q + i_no;

               double Vx_correction_factor = exp(-2*r_coord_sq/(waist_x*waist_x));
               double Vy_correction_factor = exp(-2*r_coord_sq/(waist_y*waist_y));
               double Vz_correction_factor = exp(-2*r_coord_sq/(waist_z*waist_z));

#else
               double Vx_correction_factor = 1;
               double Vy_correction_factor = 1;
               double Vz_correction_factor = 1;

#endif

               double t_x, t_y, t_z, U_x, U_y, U_z, U_on, ep;     // U_x , U_y , U_z are raw U !!!

               for (int i=0; i <= half_N_x; i++)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -PI+(PI/N_x)+((2*PI/N_x)*i),V0_0_x*Vx_correction_factor);
               t_x = calculate_t() * Er_to_nK_x;
               setup_normalized_wannier();
               U_x = calculate_U_raw();

               for (int i=0; i <= half_N_x; i++)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -PI+(PI/N_x)+((2*PI/N_x)*i),V0_0_y*Vy_correction_factor);
               t_y = calculate_t() * Er_to_nK_y;
               setup_normalized_wannier();
               U_y = calculate_U_raw();

               for (int i=0; i <= half_N_x; i++)     solve_bloch_eigenstate_and_store_into_global_bloch(i, -PI+(PI/N_x)+((2*PI/N_x)*i),V0_0_z*Vz_correction_factor);
               t_z = calculate_t() * Er_to_nK_z;
               setup_normalized_wannier();
               U_z = calculate_U_raw();

/*
               t_z = t_y;
               U_z = U_y;
*/

               U_on = calculate_U(U_x,U_y,U_z);

#ifdef SITE_DEPENDENT_T_U
               ep = VT_x * x_coord_sq + VT_y * y_coord_sq + VT_z * z_coord_sq;

               if ((i_no == (N_x-1)) || (j_no == (N_x-1)) || (k_no == (N_x-1)))  { ep = 100000000.; }

               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s + q + i_no_conj;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s + q_conj + i_no;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s + q_conj + i_no_conj;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s_conj + q + i_no;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s_conj + q + i_no_conj;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s_conj + q_conj + i_no;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

               p = s_conj + q_conj + i_no_conj;
               t_x_raw[p] = t_x;
               t_y_raw[p] = t_y;
               t_z_raw[p] = t_z;
               U_[p]      = U_on;
               epsilon_[p]= ep;

#else
               std::cout << "V0 (x) = " << V0_0_x << "  ;  lambda_x = " << lambda_x << "  ;  t (x) = " << t_x << std::endl;
               std::cout << "V0 (y) = " << V0_0_y << "  ;  lambda_y = " << lambda_y << "  ;  t (y) = " << t_y << std::endl;
               std::cout << "V0 (z) = " << V0_0_z << "  ;  lambda_z = " << lambda_z << "  ;  t (z) = " << t_z << std::endl;
               std::cout << "U = " << U_on << std::endl;
               std::cout << "U/t (x) = " << (U_on / t_x) << std::endl;
               std::cout << "U/t (y) = " << (U_on / t_y) << std::endl;
               std::cout << "U/t (z) = " << (U_on / t_z) << std::endl;
               std::cout << std::endl;

               for (uint32_t p=0; p < N_x_cube; p++) {
                  t_x_raw[p] = t_x;
                  t_y_raw[p] = t_y;
                  t_z_raw[p] = t_z;
                  U_[p]      = U_on;
               }

               for (int k_no = 0; k_no < N_x; k_no++) {
                  uint32_t s      = k_no * N_x_sq;
                  double z_coord_sq  = x[k_no] * x[k_no];
                  for (int j_no = 0; j_no < N_x; j_no++) {
                     uint32_t q      = j_no * N_x;
                     double y_coord_sq = x[j_no] * x[j_no];
                     for (int i_no= 0; i_no < N_x; i_no++) {
                        double x_coord_sq = x[i_no] * x[i_no];
                        uint32_t p = s + q + i_no;
                        ep = VT_x * x_coord_sq + VT_y * y_coord_sq + VT_z * z_coord_sq;
                        if ((i_no == 0) || (j_no == 0) || (k_no == 0) || (i_no == (N_x-1)) || (j_no == (N_x-1)) || (k_no == (N_x-1)))  { ep = 100000000.; }
                        epsilon_[p]= ep;
                     }
                  }
               }

#endif
#ifdef SITE_DEPENDENT_T_U
            }
         }
      }
#endif
   }

#ifndef DEBUGMODE
   std::cout << "Completed totally!" << std::endl;
#endif

#ifndef DEBUGMODE
   for (int k_no = 0; k_no < N_x; k_no++) {
      uint32_t s      = k_no * N_x_sq;
      for (int j_no = 0; j_no < N_x; j_no++) {
         uint32_t q      = j_no * N_x;
         for (int i_no = 0; i_no < N_x; i_no++) {
            uint32_t p = s + q + i_no;
            t_x_plus_[p] = (i_no == N_x_minus) ? 0 : 0.5*(t_x_raw[p] + t_x_raw[p+1]); 
            t_y_plus_[p] = (j_no == N_x_minus) ? 0 : 0.5*(t_y_raw[p] + t_y_raw[p+N_x]);
            t_z_plus_[p] = (k_no == N_x_minus) ? 0 : 0.5*(t_z_raw[p] + t_z_raw[p+N_x_sq]); 
            t_x_minus_[p] = (i_no == 0) ? 0 : 0.5*(t_x_raw[p] + t_x_raw[p-1]);
            t_y_minus_[p] = (j_no == 0) ? 0 : 0.5*(t_y_raw[p] + t_y_raw[p-N_x]);
            t_z_minus_[p] = (k_no == 0) ? 0 : 0.5*(t_z_raw[p] + t_z_raw[p-N_x_sq]);
         }
      }
   }
#endif

   //for (int i=0; i < N_x; i++)   print_bloch(i);
   //print_wannier();
   //print_BH_parm(); 


#ifndef DEBUGMODE
   }
   else {
   std::cout << "Obtained BH parameters from file!" << std::endl;
   //print_BH_parm();
   }
#endif

   // setting informtion regarding U_0, t_0 and centre_site
   centre_site = site_no(half_N_x,half_N_x,half_N_x);
}
