/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Lode Pollet <pollet@itp.phys.ethz.ch>
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


#include "lp_sse.h"

#include <cstdlib>
#include <stdexcept>
#include <boost/throw_exception.hpp>


void LP_Solve_SSE::init(const vector<double>& v) {
    int n = v.size();
    Nvar = n*(n-1)/2;
    Nslack = n;
    mWeight.resize(Nslack);
    C.resize(Nvar+Nslack);
    B.resize(Nslack);
    mVar.resize(Nslack);
    mSolution.resize(Nvar + Nslack);
    mBounce = 0.;
    A.resize(Nslack);
    for (int i = 0; i < Nslack; i++) A[i].resize(Nvar+Nslack);
        
    // set weight
    for (int i = 0; i < Nslack; i++) mWeight[i] = v[i];
    
    // set variables
    for (int i = 0; i < Nslack; i++) mVar[i] = Nvar+i;
        
    // set right hand site ( = 1)
    for (int i = 0; i < Nslack; i++) B[i] = 1.;
    
    // set solution. Assumed is canonical feasible form
    for (int i = 0; i < Nvar; i++) mSolution[i] = 0.;
  for (int i = 0; i < Nslack; i++) mSolution[Nvar+i] = B[i];
    
    // set Matrix A
    for (int i = 0; i < Nslack; i++) {
        for (int j = 0; j < Nvar+Nslack; j++) A[i][j] = 0.;
    }
    int inc_var = 0;
    for (int i = 0; i < Nslack; i++) {
        for (int j = i+1; j < Nslack; j++) {
            // we have two indices i and j with j > i.
            // then there is a new variable x_k at [i,j] and an entrance of x_k weight[i] / weight[j] at [j,i] in the stoch matrix
            // so for matrix A there are 2 non-zero entries in the k-th column, at row i (which is 1) and at row j (which is w[i]/w[j])
            A[i][inc_var] = 1.;
            A[j][inc_var] = mWeight[i]/mWeight[j];
            inc_var++;
        }
        // identity matrix for slack variables:
        A[i][Nvar + i] = 1.;
    }
    
  // set object function C (= -Trace of stochastic matrix), which is minus the sum of the column of A
    // for the real variables, and zero for the slack variables. Let's do brute force without caring much.
    for (int j = 0; j < Nvar; j++) {
        C[j] = 0.;
        for (int i = 0; i < Nslack; i++) {
            C[j] -= A[i][j];
        }
    }
    for (int j = 0; j < Nslack; j++) C[Nvar+j] = 0.;

    // all zero's for the variables is a feasible solution, so set Bounce to 0.
    mBounce = 0.;
    
    
    /*
    // Test Example:
    C[0] = -3; C[1] = -1; C[2] = -3;
    B[0] =  2; B[1] = +5; B[2] = +6;
    A[0][0] = 2; A[0][1] = 1; A[0][2] = 1;
    A[1][0] = 1; A[1][1] = 2; A[1][2] = 3;
    A[2][0] = 2; A[2][1] = 2; A[2][2] = 1;
    */  
#ifdef DEBUGMODE
    print_tableau();
#endif
    
}

void LP_Solve_SSE::calc_solution(vector<double>& transprob) {
    int r,s;
    int inc_iter = 0;
    const int max_iteration = 5*(Nslack+Nvar);
    for (;;) {
        // 1. check whether optimal solution has been found
        if (check_optimal()) {
            // put solution into stochastic matrix form
            copy_solution(transprob);
            return;
        }
#ifdef DEBUGMODE        
        cout << "# Solution is not optimal.\n";
#endif
        // 2. find new pivot column
        s = find_pivot_column();
        
        // 3. find new pivot row
        r = find_pivot_row(s);
#ifdef DEBUGMODE
        cout << "# Entering variable : " << s << "\n";
        cout << "# Leaving variable  : " << mVar[r] << "\n";
        cout << "# New pivot element (row, col) : " << r << "\t" << s << "\tvalue: " << A[r][s] << "\n";
#endif
        
        // 4. do pivot operation
        do_pivot(r,s);
#ifdef DEBUGMODE        
        print_tableau();
#endif
        inc_iter++;
        if (inc_iter >= max_iteration) {
            cerr << "# Simplex method failed to converge. Please check.\n";
            return;
        }
#ifdef DEBUGMODE
        char ch; cin >> ch;
#endif
    }
}


void LP_Solve_SSE::print_tableau() const {
    cout << "\n\n\n# Tableau :\n";
    cout << "# Variables : ";
    for (int i = 0; i < Nvar; i++) cout << mVar[i] << "\t";
    cout << "\n";
    cout << "# Solution : ";
    for (int i = 0; i < Nvar + Nslack; i++) cout << mSolution[i] << "\t";
    cout << "\n\n";
  for (int i = 0 ; i < Nslack; i++) {
        cout << mVar[i] << " | ";
  for (int j = 0; j < Nvar; j++) {
    cout << A[i][j] << "\t";
  }
  cout << "|\t";
  for (int j = 0; j < Nslack; j++) {
     cout << A[i][j+Nvar] << "\t";
  }
  cout << "|\t" << B[i] << "\n";
  }
  for (int j = 0; j < Nvar + Nslack+1; j++) cout << "----------";
  cout << "\n";
    cout << " " << " | ";
  for (int j = 0; j < Nvar; j++) cout << C[j] << "\t";
  cout << "|\t";
    for (int j = 0; j < Nslack; j++) cout << C[Nvar + j] << "\t";
    cout << "|\t";
  cout << mBounce << "\n";
}

bool LP_Solve_SSE::check_optimal() {
  for (int j = 0; j < Nvar+Nslack; j++) {
  if (C[j] < 0) return(false);
  }
  return(true);
}


void LP_Solve_SSE::do_pivot(const int r, const int s) {
    //value of pivot element
    double p = A[r][s];
    
    // normalize row r
    for (int k = 0; k < Nvar+Nslack; k++) A[r][k] /= p;
    B[r] /= p;
             
    // update minbounce
    mBounce -= C[s]*B[r];
    
    // update C
    for (int k =0; k < Nvar+Nslack; k++) {
        if(k != s) C[k] -= A[r][k]*C[s];
    }
    C[s] = 0.;
    
    // update B
    for (int m = 0; m < Nslack; m++) {
        if (m != r) B[m] -= B[r]*A[m][s];
    }
    
    
    // update Matrix A
    for (int i = 0; i < Nslack; i++) {
        if (i != r) {
          for (int j = 0; j < Nvar+Nslack; j++) {
                if (j != s) A[i][j] -= A[r][j]*A[i][s];
          }
        }
    }
    
    for (int i = 0; i < Nslack; i++) {
        A[i][s] = (i == r ? 1. : 0.); 
    }

    // update variables
    mVar[r] = s;
        
    // update Solution
    for (int i = 0; i < Nvar + Nslack; i++) mSolution[i] = 0.;
    for (int j = 0; j < Nslack; j++) mSolution[mVar[j]] = B[j];    
    
}


int LP_Solve_SSE::find_pivot_column() const {
    // convention : take most negative coefficient in C
    //cout << "# Which column for pivot? "; int k; cin >> k;
    //return(k);
    int imin = -1;
    double valmin = 1.0;
    for (int i = 0; i < Nvar + Nslack; i++) {
        if (C[i] < valmin) {
            valmin = C[i];
            imin = i;
        }
    }
    return(imin);
}


int LP_Solve_SSE::find_pivot_row(const int s ) const {
    // compute the ratios
    double ratio = 0.;
    double minratio = -1.;
    int min_row = -1;
    for (int i = 0; i < Nslack; i++) {
        if (A[i][s] > 0) {
            ratio = B[i] / A[i][s];
            if (min_row < 0 || ratio < minratio) {
                min_row = i;
                minratio = ratio;
            }
        }
    }
    if (min_row < 0)
        boost::throw_exception(std::runtime_error("# Error in find_pivot_row"));
    return(min_row);
}


void LP_Solve_SSE::copy_solution(vector<double>& transprob) {
    transprob.resize(Nslack*Nslack);
    
    int inc = 0;
    
    // upper triangle
    for (int i = 0; i < Nslack; i++) {
        for (int j = i+1; j < Nslack; j++) {
            transprob[i * Nslack + j] = mSolution[inc];
            inc++;
        }
    }
    
    // lower triangle
    for (int i = 0; i < Nslack; i++) {
        for (int j = 0; j < i; j++) {
            transprob[i * Nslack + j] = transprob[j*Nslack + i] * mWeight[j] / mWeight[i];
        }
    }
    
    // diagonal
    for (int i = 0; i < Nslack; i++) {
        double s = 0.;
        for (int j = 0; j < Nslack; j++) {
            if (i != j) s += transprob[i*Nslack + j];
        }
        transprob[i*Nslack + i] = 1. - s;
    }
}


double LP_Solve_SSE::calc_loc_opt_bounce() {
    // this is for test purposes only and not part of the program
    
    // order weights
    int n = Nslack;
    vector<double> w(n);
    for (int i =0; i < n; i++) w[i] = mWeight[i];
    for (int i =0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (w[j] > w[i]) {
                double s = w[j];
                w[j] = w[i];
                w[i] = s;
            }
        }
    }
    
    // normalize weights
    double s = 0.;
    for (int i = 0; i < n; i++) s += w[i];
    for (int i = 0; i < n; i++) w[i] /= s;
    
    // check
    for (int i =0; i < n-1; i++)
        if (w[i] > w[i+1]) 
            boost::throw_exception(std::runtime_error("Order weights"));
    
    
    if (n == 2) {
        return(1.-w[0]/w[1]);
    }
    else if (n==3) {
        double y1 = w[0]/(1-w[0]);
        double y2 = (1-y1)*w[1]/(1-w[0]-w[1]);
        return(1-y1-y2);
    }
    else if (n==4) {
        double y1 = w[0]/(1-w[0]);
        double y2 = (1-y1)*w[1]/(1-w[0]-w[1]);
        double y3 = (1-y1-y2)*w[2]/w[3];
        return(1-y1-y2-y3);
    }
    else
        boost::throw_exception(std::runtime_error("# Not implemented"));
}


