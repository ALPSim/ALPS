#include <bitset>
#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <boost/array.hpp>

using std::cerr;
using std::cout;
using std::endl;

// Hamiltonian of a simple one-particle system
template<class Vector>
class Hamiltonian
{
    typedef unsigned int state;
    
public:
    Hamiltonian(unsigned int L_) : L(L_) { }
    
    void mult(Vector const & vin, Vector & vout) const
    {
        ++count;
        
        vout = vin;
        fill(vout.begin(), vout.end(), 0);
        
        for (unsigned int s = 0; s < L; ++s) {
            vout[(s+1)%L] += -vin[s];
            vout[((s+L)-1)%L] += -vin[s];
        }
    }
    
    int get_count() const
    {
        return count;
    }
    
    void reset_count()
    {
        count = 0;
    }
    
    unsigned int basis_size() const { return L; }
    
private:
    unsigned int L;
    mutable int count;
};

// The ietl eigensolvers require an operator ietl::mult( A, x, b ) calculating b = A x,
// it has to be declared before the eigensolver is included.
// Look at the documentation to see a full list of requirements.
//[ example1
namespace ietl {
    template<class Vector>
    void mult(Hamiltonian<Vector> const & H, Vector& x, Vector& y)
    {
        H.mult(x, y);
    }
}

#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/jd.h>
#include <boost/random.hpp>

typedef boost::numeric::ublas::vector<double> vector_t;
typedef ietl::vectorspace<vector_t> vecspace_t;
typedef boost::lagged_fibonacci607 gen_t;
typedef boost::lagged_fibonacci44497 gen2_t;
typedef Hamiltonian<vector_t> ham_t;

int main(int argc, char **argv)
{
    int L = 100;    // dimension of the Hamiltonian
    //<-
    if( argc > 1 )
        L = atoi(argv[1]);
    //->
    
    vecspace_t vs(L);
    ham_t H(L);

    // random generator to create a starting vector x0
    // to find an eigenvector v, we require dot(x0,v) != 0
    gen_t gen;
    gen2_t gen2;
    
    int n_evals = 10;   // number of eigenpairs to be calculated
    int max_iter = L*10; // maximal number of iterations

    int m_max = 40;
    int m_min = 20;
 
    // tolerance
    double rel_tol = sqrt(std::numeric_limits<double>::epsilon());
    double abs_tol = rel_tol;

    // maximal iterations for the correction equation solver
    unsigned max_cor_iter = 10;
    // on default 5 steps are used

    ietl::jd_iteration<double> iter(max_iter, m_min, m_max, rel_tol, abs_tol);
    ietl::jd<Hamiltonian<vector_t>, vecspace_t> jd(H, vs /*, 2 for verbose mode */ );

    // the correction equation solver must be an function object
    ietl::gmres_wrapper gmres(max_cor_iter);
    
    // to find degenerated (or not well seperated) eigenvalues in the right order
    // the correction equation has to be solved exactly, this is very expensive. 
    // here we use different starting vectors instead.
    try {
            jd.eigensystem(iter, gen, n_evals/2, gmres);
            jd.eigensystem(iter, gen2, n_evals-n_evals/2, gmres);
    } catch (std::runtime_error& e) {
            cerr << "Error in eigenpair calculation: " << e.what() << "\n";
            return 1;
    }
    // if you use std=c++0x gmres is used as default

    std::vector<double> evals = jd.eigenvalues(); // copy eigenvalues
    // cout << "First eigenvector: \n" << jd.eigenvector(0) << "\n";
    std::cout.precision(10);
    std::sort(evals.begin(), evals.end());
    cout << "Sorted Eigenvalues: \n";
    std::copy(evals.begin(), evals.end(), std::ostream_iterator<double>(cout, "\n"));
    //<-
    // calculate some more eigenvalues
    ietl::bicgstab_wrapper<double,2> bicgstab;

    try {
        jd.eigensystem(iter, gen, n_evals, bicgstab);
    } catch (std::runtime_error& e) {
            cerr << "Error in eigenpair calculation: " << e.what() << "\n";
            return 1;
    }

    evals = jd.eigenvalues();
    cout << "Unsorted Eigenvalues: \n";
    std::copy(evals.begin(), evals.end(), std::ostream_iterator<double>(cout, "\n"));

    return 0;
    //->
}
    //] example1
