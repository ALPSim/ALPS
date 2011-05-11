// g++ -o test test_jd.cpp -I.. -llapack -fopenmp
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/iteration.h>
#include <ietl/vectorspace.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <exception>
#include <cassert>


namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> vector_type;
typedef ublas::matrix<double, ublas::row_major > matrix_type;


namespace ietl  {
    void mult (const matrix_type& A, const vector_type& x, vector_type& b)
    {
        assert( A.size2() == x.size() );
        if( b.size() != A.size1() )
            b = vector_type( A.size1() );

        #pragma omp parallel
        {
            #pragma omp for
            for(std::size_t i = 0; i<b.size(); ++i)
                b(i) = inner_prod( row(A, i) , x);
        }
    }
}

#include <ietl/jd.h>

class crandom {
    public:
        typedef boost::lagged_fibonacci44497 engine_type;
        typedef boost::normal_distribution<double> dist_type;
        typedef boost::uniform_real<double> udist_type;

        crandom()
            : dist_(0./*mean*/,10./*var*/), udist_(0.,1.), engine_(), rng_(engine_,dist_), 
                urng_(engine_,udist_) {}

        double operator() ()
            {
                return rng_();
            }
        double u ()
            {   return urng_();     }
    private:
        engine_type engine_;
        dist_type dist_;
        udist_type udist_;
        boost::variate_generator<engine_type&,dist_type> rng_;
        boost::variate_generator<engine_type&,udist_type> urng_;
};


int main ()
{

    int N = 1000;
    matrix_type mat(N, N);

    crandom rng;
    int n = 1;
    for(int i=0;i<N;i++)
        for(int j=0;j<=i;j++)
            //mat(i,j) = n++;
            mat(i,j) = mat(j,i) = rng();  

    //std::cout << "Matrix: " << mat << std::endl;
  
    ietl::vectorspace<vector_type> vs(N);
    boost::lagged_fibonacci607 gen;  
    int max_iter = 1000;
    double rel_tol = 1e-7;
    double abs_tol = 1e-8;
    int m_min = 20, m_max = 40;

    ietl::jd_iteration<double> iter(m_min, m_max, max_iter, rel_tol, abs_tol );

    ietl::basic_iteration<double> giter(max_iter, rel_tol, 1e-2);

    ietl::jd<matrix_type, ietl::vectorspace<vector_type> > jd_test(mat, vs, 2);

    int k = 3;

    //search lowest eigenvalue
    try{
        jd_test.eigensystem(iter, gen, 1);
    }
    catch (std::runtime_error& e) {
        std::cerr << "Something went wrong: " << e.what() << "\n";
    }
    //search some more eigenvalues
    try{
        jd_test.eigensystem(iter, gen, k);
    }
    catch (std::runtime_error& e) {
        std::cerr << "Something went wrong: " << e.what() << "\n";
    } 

    std::cout.precision(10);  
    for(int i = 0; i< jd_test.eigenvalues().size(); ++i)
        std::cout <<"eigenpair # "<<i <<"\n"<<jd_test.eigenvalue(i)<<"\n"/*<<jd_test.eigenvector(i)<<*/"\n\n";

    jd_test.reset();

    //some upper eigenvalues
    try{
        jd_test.eigensystem(iter, gen, k, true);
    }
    catch (std::runtime_error& e) {
        std::cerr << "Something went wrong: " << e.what() << "\n";
    } 

    //std::cout <<"filter ghosts...\n";
    //jd_test.filter_ghosts(giter);

    for(int i = 0; i< jd_test.eigenvalues().size(); ++i)
        std::cout <<"eigenpair # "<<i <<"\n"<<jd_test.eigenvalue(i)<<"\n"/*<<jd_test.eigenvector(i)*/<<"\n\n";

    return 0;
}
