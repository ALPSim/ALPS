#ifndef JACOBI_DAVIDSON_H
#define JACOBI_DAVIDSON_H
#include <ietl/traits.h>
#include <ietl/complex.h>
#include <ietl/gmres.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <exception>
#include <limits>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
//for lapack::heev
#include <boost/numeric/bindings/lapack/driver.hpp>
//for lapack::getrs
#include <boost/numeric/bindings/lapack/computational.hpp>
//#include <boost/numeric/bindings/detail/config/fortran.hpp>

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/hermitian.hpp>
#include <boost/numeric/bindings/std.hpp>

namespace ietl{

    namespace ublas=boost::numeric::ublas;
    namespace lapack=boost::numeric::bindings::lapack;

//jd_iterator//////////////////////////////////////////////////////
    template <class T>
    class jd_iteration : public basic_iteration<T> {
        public: 

        jd_iteration(size_t m_min, size_t m_max, size_t max_iter, T reltol = 0., T abstol = 0.) 
            : basic_iteration<T>(max_iter, reltol, abstol), m_min_(m_min), m_max_(m_max) {}

        inline size_t m_min() const
            {   return m_min_;  }
        inline size_t m_max() const
            {   return m_max_;  }
        const std::string& error_msg() 
            { return basic_iteration<T>::err_msg; }    
        private:
        size_t m_min_, m_max_;
    };
//jd_iterator//////////////////////////////////////////////////////
    namespace detail{

        template <class MATRIX, class VS> 
            class deflated_matrix;
        template <class MATRIX, class VS, class VECTOR>
            void mult(const deflated_matrix<MATRIX,VS> & Adef, const VECTOR& v, VECTOR& r);

        template <class MATRIX, class VS>
        class deflated_matrix   {
            public:
                typedef typename vectorspace_traits<VS>::vector_type vector_type;
                typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
                typedef typename real_type<scalar_type>::type real_type;

                deflated_matrix(const MATRIX& A, double theta, const std::vector<vector_type>& Q) 
                    : A_(A) , theta_(theta) , Q_(Q) {}
                ~deflated_matrix() {}

                void set_theta(double theta) 
                    { theta_ = theta; }

                friend void mult<>(const deflated_matrix& A_def, const vector_type& v, vector_type& r);

            private:
                const MATRIX& A_;
                real_type theta_;
                const std::vector<vector_type>& Q_;    
        };

        template <class MATRIX, class VS, class VECTOR>
        void mult(  const deflated_matrix<MATRIX,VS>& Adef, 
                    const VECTOR& v, 
                    VECTOR& r)
        {// (I-QQ*)(A-thetaI)(I-QQ*)v = r 
                VECTOR temp = v;

                for(size_t i = 0; i < Adef.Q_.size() ; ++i) 
                    temp -= ietl::dot(Adef.Q_[i], temp) * Adef.Q_[i];

                ietl::mult(Adef.A_, temp, r);
                r -= Adef.theta_ * temp;

                for(size_t i = 0; i < Adef.Q_.size() ; ++i) 
                    r -= ietl::dot(Adef.Q_[i], r) * Adef.Q_[i];
        }

        template <class MATRIX, class VS, class PREC> 
            class deflated_matrix_prec;
        template <class MATRIX, class VS, class PREC, class VECTOR>
            void mult(const deflated_matrix_prec<MATRIX,VS,PREC> & Adef, const VECTOR& v, VECTOR& r);

        template <class MATRIX, class VS, class PREC>
        class deflated_matrix_prec   {
            public:
                typedef typename vectorspace_traits<VS>::vector_type vector_type;
                typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
                typedef typename real_type<scalar_type>::type real_type;

                deflated_matrix_prec(const PREC& K, const MATRIX& A, double theta, const std::vector<vector_type>& Q) 
                    : A_(A) , theta_(theta) , Q_(Q), K_(K) {}
                ~deflated_matrix_prec() {}

                void set_theta(real_type theta) 
                    { theta_ = theta; }

                friend void mult<>(const deflated_matrix_prec& A_def, const vector_type& v, vector_type& r);

            private:
                const MATRIX& A_;
                const PREC& K_;
                real_type theta_;
                const std::vector<vector_type>& Q_;    
        };

        template <class MATRIX, class VS, class PREC, class VECTOR>
        void mult(  const deflated_matrix_prec<MATRIX,VS,PREC>& Adef, 
                    const VECTOR& v, 
                    VECTOR& r)
        {// (I-QQ*)K(A-thetaI)(I-QQ*)v = r 
                VECTOR temp = v;

                for(size_t i = 0; i < Adef.Q_.size() ; ++i) 
                    temp -= ietl::dot(Adef.Q_[i], temp) * Adef.Q_[i];

                ietl::mult(Adef.A_, temp, r);
                r -= Adef.theta_ * temp;

                ietl::mult(Adef.K_, r, r);

                for(size_t i = 0; i < Adef.Q_.size() ; ++i) 
                    r -= ietl::dot(Adef.Q_[i], r) * Adef.Q_[i];
        }

        template <class REAL>
        class mycmp{
            public:
                typedef std::pair<size_t, REAL> pair_type;

                mycmp(REAL t)
                    : tau_(t) {}

                inline bool operator() (pair_type i, pair_type j)
                    { return ( std::abs(i.second-tau_) < std::abs(j.second-tau_) ); }
            private:
            REAL tau_;
        };

        template <class REAL>
        inline REAL mysecond (std::pair<size_t,REAL> i)
            { return i.second; }

}//end namespace detail
//solver//////////////////////////////////////////////////////////////////////////////////////
namespace solver {

    template <class MATRIX, class VS, class PREC>
    class left_prec_solver {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename real_type<scalar_type>::type real_type;
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        left_prec_solver ( const MATRIX& A, const std::vector<vector_type>& Q, const VS& vs, PREC& K) 
            : A_(A) , Q_(Q) , vspace_(vs) , K_(K) , x0_(new_vector(vs))  {}

        template <class IT>
        void operator() ( scalar_type theta, vector_type& r, vector_type& t, IT& iter );

    protected:
        const MATRIX& A_;
        const std::vector<vector_type>& Q_;
        const VS& vspace_;
        std::vector<vector_type> Q_hat_;
        PREC& K_;
        vector_type x0_;
    };

    // K preconditioner for (A - theta *I)
    template <class MATRIX, class VS, class PREC> template <class IT>
    void left_prec_solver<MATRIX,VS,PREC>::operator() ( scalar_type theta, vector_type& r, vector_type& t, IT& iter )
    {
        for(size_t i = 0; i < vec_dimension(vspace_); ++i)
            K_(i,i) = 1./(A_(i,i) - theta); // inv(K) instead of K

        const size_t max_iter = 5;
        size_t m = Q_.size();
        size_t m_old = 0;//Q_hat_.size();
        double abs_tol = iter.absolute_tolerance();
        abs_tol *= abs_tol;
        int info;
        matrix_type M(m,m);
        matrix_type LU;
        ublas::vector<scalar_type> b(m);
        std::vector<int> pivot(m);
        ublas::vector<scalar_type> gamma(m);
        detail::deflated_matrix_prec<MATRIX,VS,PREC> Adef(K_, A_, theta, Q_);

        Q_hat_.resize(m);

        for(size_t i = m_old; i < m; ++i)
        {
            // starting with x0 = b is good if K_ approx. I
            //Q_hat_[i] = ietl_gmres(K_, Q_[i], Q_[i], max_iter, abs_tol, false);
            ietl::mult(K_, Q_[i], Q_hat_[i]);
        }

        for(size_t i = 0; i < m; ++i)
            for(size_t j = ( (i<m_old) ? m_old : 0 ); j < m; ++j)
                M(i,j) = ietl::dot(Q_[i], Q_hat_[j]);

        LU = M;
        info = lapack::getrf(LU, pivot);
            if(info != 0)   throw std::runtime_error("lapack::getrf failed.");

        //compute r_tilde
        //t = ietl_gmres(K_, r, r, max_iter, abs_tol, false);
        ietl::mult(K_, r, t)
;
        for(size_t i = 0; i < m; ++i)
            gamma(i) = ietl::dot(Q_[i],t);

        //calculate alpha from M alpha = gamma
        info = lapack::getrs(LU, pivot, gamma);
            if(info != 0)   throw std::runtime_error("lapack::getrs failed.");

        for(size_t i = 0; i < m; ++i)
            t -= Q_hat_[i]*gamma(i);

         t *= -1;

        //apply krylov solver for A_tilde v = -r_tilde
        r = ietl_gmres(Adef, t, x0_, max_iter, abs_tol, false);

        //r = y, t = y_hat
        ietl::mult(A_, r, t);
        r *= -theta;
        r += t;

        //t = ietl_gmres(K_, r, r, max_iter, abs_tol, false);
        ietl::mult(K_, r, t);

        for(size_t i = 0; i < m; ++i)
            gamma(i) = ietl::dot(Q_[i],t);

        info = lapack::getrs(LU, pivot, gamma);
            if(info != 0)   throw std::runtime_error("lapack::getrs failed.");

        //z = y_hat - Q_hat alpha
        for(size_t i = 0; i < m; ++i)
            t -= Q_hat_[i]*gamma(i);

    }//left_prec_solver::void()

///////simple left preconditioned solver///////////////////////////////////////////
//solving (I-QQ*)K(A-theta*I)(I-QQ*)t = -Kr
    template <class MATRIX, class VS, class PREC>
    class left_prec_simple {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename real_type<scalar_type>::type real_type;
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        left_prec_simple ( const MATRIX& A, const std::vector<vector_type>& Q, const VS& vs, PREC& K) 
            : A_(A) , Q_(Q) , vspace_(vs) , K_(K), Adefp_(K, A, 0., Q) , x0_(new_vector(vs))  {}

        template <class IT>
        void operator() ( scalar_type theta, vector_type& r, vector_type& t, IT& iter );

    protected:
        const MATRIX& A_;
        PREC& K_;
        detail::deflated_matrix_prec<MATRIX, VS, PREC> Adefp_;
        const std::vector<vector_type>& Q_;
        const VS& vspace_;
        std::vector<vector_type> Q_hat_;
        vector_type x0_;
    };

    template <class MATRIX, class VS, class PREC> template <class IT>
    void left_prec_simple<MATRIX,VS,PREC>::operator() ( scalar_type theta, vector_type& r, vector_type& t, IT& iter )
    {
        Adefp_.set_theta(theta);

        for(size_t i = 0; i < vec_dimension(vspace_); ++i)
        {
            K_(i,i) = 1./(A_(i,i) - theta);
        }

        r *= -1;
        //ietl::mult(K_, r, r);

        t = ietl_gmres(Adefp_, r, x0_, 5, iter.absolute_tolerance(), false);

    }


//////////simple solver with gmres/////////////////////////
    template <class MATRIX, class VS>
    class jd_solver {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename real_type<scalar_type>::type real_type;

        jd_solver (const MATRIX& A, const std::vector<vector_type>& Q, const VS& vspace)
            : Adef_(A, 0, Q), vspace_(vspace)
             {      x0_ = new_vector(vspace);    }
        
        template <class IT>
        void operator() ( scalar_type theta, vector_type& r, vector_type& t, IT& iter );

    protected:
        detail::deflated_matrix<MATRIX, VS> Adef_;
        const VS& vspace_;
        vector_type x0_;
    };
        template <class MATRIX, class VS>
        template <class IT>
        void jd_solver<MATRIX, VS>::operator()( scalar_type theta, vector_type& r, vector_type& t, IT& iter)
        {
            Adef_.set_theta(theta);
            r *= -1;
            t = ietl_gmres(Adef_, r, x0_, 5, iter.absolute_tolerance(), false); 
        }
}//end namespace solver/////////////////////////////////////////////////////////////////////////////////

    template <class MATRIX, class VS>
    class jd {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename vectorspace_traits<VS>::size_type size_type;
        typedef typename vectorspace_traits<VS>::magnitude_type magnitude_type;
        typedef typename real_type<scalar_type>::type real_type;
        typedef std::vector<vector_type> vector_set_type;
        typedef std::vector<real_type> real_set_type;

        jd (MATRIX& A, VS& vspace, size_t v = 0)
            : A_(A), vspace_(vspace), n_(vec_dimension(vspace)), verbose_(v) {}

        //default: search for lowest eigenvalues
        template <class SOLVER, class IT, class GEN>
        void eigensystem(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, bool search_highest);
        //with a target value tau
        template <class SOLVER, class IT, class GEN>
        void eigensystem(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, real_type tau);

        //without preconditioning
        template <class IT, class GEN>
        void eigensystem(IT& iter, GEN& gen, size_type k_max, bool search_highest = false)    
            {
                typedef solver::jd_solver<MATRIX,VS> solver_type;
                solver_type default_solver(A_, X_, vspace_);
                eigensystem <solver_type,IT,GEN> (default_solver, iter, gen, k_max, search_highest);
            }

        template <class IT, class GEN>
        void eigensystem(IT& iter, GEN& gen, size_type k_max, real_type tau)    
            {
                typedef solver::jd_solver<MATRIX,VS> solver_type;
                solver_type default_solver(A_, X_, vspace_);
                eigensystem <solver_type,IT,GEN> (default_solver, iter, gen, k_max, tau );
            }

        //with left preconditioner K, for (A - theta * I)
        template <class IT, class GEN, class PREC>
        void eigensystem( IT& iter, GEN& gen, size_type k_max, PREC& K, bool search_highest = false)
            {
                typedef solver::left_prec_simple<MATRIX,VS,PREC> solver_type;
                solver_type lp_solver(A_, X_, vspace_, K);
                eigensystem <solver_type,IT,GEN> (lp_solver, iter, gen, k_max, search_highest);
            }

        //jd with harmonic ritz values
        template <class SOLVER, class IT, class GEN>
        void eigensystem_harmonic(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, real_type tau);
        template <class IT, class GEN>
        void eigensystem_harmonic(IT& iter, GEN& gen, size_type k_max, real_type tau)    
            {
                typedef solver::jd_solver<MATRIX,VS> solver_type;
                solver_type default_solver(A_, X_, vspace_);
                eigensystem_harmonic <solver_type,IT,GEN> (default_solver, iter, gen, k_max, tau);
            }

        //access functions
        real_type eigenvalue(size_type k)
            {
                assert (k < Lambda_.size());
                return Lambda_[k];
            }
        real_set_type eigenvalues()
            {
                return Lambda_;
            }
        const vector_type& eigenvector (size_type k)
            {
                assert (k < X_.size());
                return X_[k];
            }
        void reset()
            {
                X_.clear();
                Lambda_.clear();
            }
        template <class IT>
        void filter_ghosts(IT& iter);

    protected:
        MATRIX& A_;
        VS& vspace_;
        size_type n_;
        vector_set_type X_;
        real_set_type Lambda_;
        size_t verbose_;    //verbosity levels: 0 = no, 1 = yes, 2 = all
    };// end of class jd

    // filter ghosts after(!) applying the jd-algorithm, because one way to prevent new
    // ghosts is using gram-schmidt repeatedly, and having the ghost-eigenvectors in X_,
    // does just that
    // might not filter all ghosts, especially if the eigenvector is close to 0
    template <class MATRIX, class VS>
    template <class IT>
    void jd<MATRIX,VS>::filter_ghosts(IT& iter)
    {
        for (size_t i = 1; i < Lambda_.size(); ++i)
            for (size_t j = 0; j < i; ++j)
                if( std::abs( Lambda_[i] - Lambda_[j] ) < iter.absolute_tolerance() ) 
                    if( ietl::two_norm( X_[i] - X_[j] ) < iter.absolute_tolerance() )   //eigenpair i is equal to j
                    {
                        if( verbose_ > 1)
                            std::cout <<"removing eigenvalue... "<< Lambda_[i] << "\n";
                        //remove i
                        Lambda_.erase( Lambda_.begin() + i);
                        for( size_t k = i; k < Lambda_.size()-1 ; ++k)
                        {
                            std::swap( X_[k], X_[k+1] );
                        }
                        X_.pop_back();
                        --i;
                    }
    }

//////jd for exterior eigenvalues
    template <class MATRIX, class VS>
    template <class SOLVER, class IT, class GEN>
    void jd<MATRIX,VS>::eigensystem(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, bool search_highest)
    {
        typedef std::vector<vector_type> vector_set_type;
        typedef std::vector<real_type> real_set_type;
        typedef ublas::hermitian_matrix< scalar_type, ublas::upper > herm_matrix_type;
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        const magnitude_type kappa = 0.25; // kappa < 1
        magnitude_type norm_t, norm_r;
        const size_type m_min = iter.m_min(), m_max = iter.m_max();
        size_type k = 0, m = 0;
        int info;

        vector_set_type V;
        vector_set_type VA;
        X_.resize( X_.size()+1 ); //X_.back() == u

        herm_matrix_type M; //might be better with preset size (k_max, k_max)

        real_set_type Theta;
        matrix_type S;

        vector_type uA, r;
        vector_type t = new_vector(vspace_);

        generate(t, gen);
        project(t, vspace_);

        while(true) 
        {
            //modified Gram-Schmidt
            norm_t = ietl::two_norm(t);

            if(norm_t == 0/*< iter.absolute_tolerance()*/) {
            //if the correction solver returns a t approx 0, numerical errors might create u non orthogonal to X_i
            //tau might be too far away? which seems not to matter for the first eigenvalue!
                    std::cerr << "Correction vector t is (almost) zero:\n";
                    if(iter.first())
                        throw std::runtime_error("generating a starting vector failed.");
                    else
                        throw std::runtime_error("try to solve the correction equation more accurately.");
                    }
            for(size_type i = 0; i < m; ++i) 
                t -= ietl::dot(V[i],t) * V[i];
            if(ietl::two_norm(t) <= kappa * norm_t)
                for(size_type i = 0; i < m; ++i)
                    t -= ietl::dot(V[i],t) * V[i];

            norm_t = ietl::two_norm(t);
            V.push_back(t/norm_t);
            ++m;
            VA.resize(m);

            ietl::mult(A_, V[m-1], VA[m-1]);
            assert(m == V.size() && m == VA.size());

            M.resize(m,m);

            for(size_type i = 0; i < m; ++i)
                M(i, m-1) = ietl::dot(V[i], VA[m-1]);

            //eigendecomposition M = S THETA S*
            Theta.resize(m);

            S = M; //copy M

            info = lapack::heev('V', S, Theta); //eigenvectors in columns of S
                if(info < 0) throw std::runtime_error("lapack::heev - illegal value.");
                if(info > 0) throw std::runtime_error("lapack::heev - failed to converge.");

            real_type theta = (search_highest) ? Theta.back() : Theta[0];

            // u = V s_1
            size_t row = (search_highest) ? Theta.size()-1 : 0;
            X_.back() = V[0] * S(0,row);
            for(size_type j = 1; j < m; ++j) 
                X_.back() += V[j] * S(j,row);

            ietl::mult(A_, X_.back(), uA);
            r = uA - theta * X_.back();

            // accept eigenpairs
            norm_r = ietl::two_norm(r);     //use 1-norm

            while(iter.finished(norm_r,theta)) 
            {
                if(iter.error_code() == 1)
                    throw std::runtime_error(iter.error_msg());

                if(verbose_ > 0)
                    std::cout << "Accepted eigenpair #"<< k+1 << "\t" << theta << "\tResidual: " << norm_r << "\n";

                Lambda_.push_back( theta ); //Lambda_[k] = Theta[0];

                if(++k == k_max) return;

                --m;
                M = ublas::zero_matrix<scalar_type>(m,m);

                //use VA temp as V
                VA.resize(m);
                for(size_type i = 0; i < m; ++i)
                {
                        size_t row = (search_highest) ? (Theta.size()-2-i) : (i+1);
                        VA[i] = V[0] * S(0,row);
                        for(size_type j = 1; j < m+1; ++j)
                            VA[i] += V[j] * S(j,row);
                }

                V.resize(m);
                for(size_type i = 0; i < m; ++i)
                {
                    std::swap(VA[i], V[i]);
                    ietl::mult(A_, V[i], VA[i]);
                    size_t row = (search_highest) ? (Theta.size()-2-i) : (i+1);
                    M(i,i) = Theta[row];
                    ublas::column(S, row + ( search_highest ? +1 : -1) ) = ublas::unit_vector<scalar_type> (S.size1(), i);
                }

                Theta.erase( (search_highest) ? --Theta.end() : Theta.begin() );

                theta = (search_highest) ? Theta.back() : Theta[0];

                X_.resize(X_.size()+1);
                X_.back() = V[0];
                r = VA[0] - theta * X_.back();
                norm_r = ietl::two_norm(r);

            } // accept eigenpairs

            // restart
            if(m >= m_max) {
                if(verbose_ > 1) 
                    std::cout<<"restart...\n";
                M = ublas::zero_matrix<scalar_type>(m_min,m_min);

                VA.resize(m_min);

                for(size_type i = 1; i < m_min; ++i) {
                    size_t row = (search_highest) ? Theta.size()-1-i : i;
                    VA[i] = V[0] * S(0,row);
                    for(size_type j = 1; j < m_max; ++j)
                        VA[i] += V[j] * S(j,row);
                }
                    
                V.resize(m_min);
                for(size_type i = 1; i < m_min; ++i)
                {
                    std::swap(V[i], VA[i]);
                    ietl::mult(A_, V[i], VA[i]);
                }

                V[0] = X_.back();
                ietl::mult(A_, X_.back(), uA);
                VA[0] = uA;

                for(size_type i = 0; i < m_min; ++i)
                    M(i,i) = Theta[(search_highest) ? Theta.size()-1-i : i];

                Theta.resize(m_min);
                m = m_min;
            }// restart

            // correction equation
            solver(theta, r, t, iter); //solver is allowed to change r

            ++iter;
            if(iter.error_code() == 1)
                throw std::runtime_error(iter.error_msg());
            if(verbose_ > 1)
                std::cout << "JD iteration " << iter.iterations() << "\tresid = " << norm_r << "\n";
        }// main loop

    }// eigensystem

//////again jd for exterior eigenvalues, but with target value tau/////////////
    template <class MATRIX, class VS>
    template <class SOLVER, class IT, class GEN >
    void jd<MATRIX,VS>::eigensystem(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, real_type tau)
    {
        typedef std::vector<vector_type> vector_set_type;
        typedef std::vector<real_type> real_set_type;
        typedef ublas::hermitian_matrix< scalar_type, ublas::upper > herm_matrix_type; //todo: ublas::lower is faster
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        const magnitude_type kappa = 0.25; // kappa < 1
        magnitude_type norm_t, norm_r;
        const size_type m_min = iter.m_min(), m_max = iter.m_max();
        size_type k = 0, m = 0;
        size_type row;
        int info;

        vector_set_type V;
        vector_set_type VA;
        X_.resize( X_.size()+1 ); //X_.back() == u

        herm_matrix_type M; //might be better with preset size (k_max, k_max)

        std::vector<std::pair<size_t, real_type> > pivot;
        detail::mycmp<real_type> cmp(tau);

        real_set_type Theta;
        matrix_type S;

        vector_type uA, r;
        vector_type t = new_vector(vspace_);

        generate(t, gen);
        project(t, vspace_);

        while(true) 
        {
            //modified Gram-Schmidt
            norm_t = ietl::two_norm(t);

            if(norm_t < iter.absolute_tolerance() /*== 0*/) {
            //if the correction solver returns a t approx 0, numerical errors might create u non orthogonal to X_i
            //tau might be too far away? which seems not to matter for the first eigenvalue!
                    std::cerr << "Correction vector t is (almost) zero:\n";
                    if(iter.first())
                        throw std::runtime_error("generating a starting vector failed.");
                    else
                        throw std::runtime_error("try to solve the correction equation more accurately.");
                    }
            for(size_type i = 0; i < m; ++i) 
                t -= ietl::dot(V[i],t) * V[i];
            if(ietl::two_norm(t) <= kappa * norm_t)
                for(size_type i = 0; i < m; ++i)
                    t -= ietl::dot(V[i],t) * V[i];

            norm_t = ietl::two_norm(t);
            V.push_back(t/norm_t);
            ++m;
            VA.resize(m);

            ietl::mult(A_, V[m-1], VA[m-1]);
            assert(m == V.size() && m == VA.size());

            M.resize(m,m);

            for(size_type i = 0; i < m; ++i)
                M(i, m-1) = ietl::dot(V[i], VA[m-1]);

            //eigendecomposition M = S THETA S*
            Theta.resize(m);

            S = M; //copy M

            info = lapack::heev('V', S, Theta); //eigenvectors in columns of S
                if(info < 0) throw std::runtime_error("lapack::heev - illegal value.");
                if(info > 0) throw std::runtime_error("lapack::heev - failed to converge.");

            //sort
            pivot.resize(m);
            for(size_type i = 0; i < m; ++i) 
                pivot[i] = std::make_pair(i,Theta[i]);
            std::sort(pivot.begin(), pivot.end(), cmp);
            std::transform(pivot.begin(), pivot.end(), Theta.begin(), detail::mysecond<real_type>);

            // u = V s_1
            row = pivot[0].first;
            X_.back() = V[0] * S(0,row);
            for(size_type j = 1; j < m; ++j) 
                X_.back() += V[j] * S(j,row);

            ietl::mult(A_, X_.back(), uA);
            r = uA - Theta[0] * X_.back();

            // accept eigenpairs
            norm_r = ietl::two_norm(r);    

            while(iter.finished(norm_r,Theta[0])) 
            {
                if(iter.error_code() == 1)
                    throw std::runtime_error(iter.error_msg());

                if(verbose_ > 0)
                    std::cout << "Accepted eigenpair #"<< k+1 << "\t" << Theta[0] << "\tResidual: " << norm_r << "\n";

                Lambda_.push_back( Theta[0] ); //Lambda_[k] = Theta[0];

                if(++k == k_max) return;

                --m;
                M = ublas::zero_matrix<scalar_type>(m,m);

                //use VA temp as V
                VA.resize(m);
                for(size_type i = 0; i < m; ++i)
                {
                        row = pivot[i+1].first;
                        VA[i] = V[0] * S(0,row);
                        for(size_type j = 1; j < m+1; ++j)
                            VA[i] += V[j] * S(j,row);
                }
                V.resize(m);
                for(size_type i = 0; i < m; ++i)
                {
                    std::swap(VA[i], V[i]);
                    ietl::mult(A_, V[i], VA[i]);
                    M(i,i) = Theta[i+1];
                    ublas::column(S, pivot[i].first) = ublas::unit_vector<scalar_type> (S.size1(), i);
                }
                Theta.erase( Theta.begin() );

                X_.resize(X_.size()+1);
                X_.back() = V[0];
                r = VA[0] - Theta[0] * X_.back();
                norm_r = ietl::two_norm(r);
            } // accept eigenpairs

            // restart
            if(m >= m_max) {
                if(verbose_ > 1) 
                    std::cout<<"restart...\n";
                M = ublas::zero_matrix<scalar_type>(m_min,m_min);
                assert(M.size1() == m_min);

                VA.resize(m_min);

                for(size_type i = 1; i < m_min; ++i) {
                    row = pivot[i].first;
                    VA[i] = V[0] * S(0,row);
                    for(size_type j = 1; j < m_max; ++j)
                        VA[i] += V[j] * S(j,row);
                }
                    
                V.resize(m_min);
                for(size_type i = 1; i < m_min; ++i)
                {
                    std::swap(V[i], VA[i]);
                    ietl::mult(A_, V[i], VA[i]);
                }

                V[0] = X_.back();
                VA[0] = uA;
                for(size_type i = 0; i < m_min; ++i)    M(i,i) = Theta[i];
                Theta.resize(m_min);
                m = m_min;
            }// restart

            // correction equation
            solver(Theta[0], r, t, iter); //solver is allowed to change r

            ++iter;
            if(iter.error_code() == 1)
                throw std::runtime_error(iter.error_msg());
            if(verbose_ > 1)
                std::cout << "JD iteration " << iter.iterations() << "\tresid = " << norm_r << "\n";
        }// main loop

    }// eigensystem

//////jd for interior eigenvalues with harmonic ritz values////////////////
    template <class MATRIX, class VS>
    template <class SOLVER, class IT, class GEN >
    void jd<MATRIX,VS>::eigensystem_harmonic(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, real_type tau)
    {
        typedef std::vector<vector_type> vector_set_type;
        typedef std::vector<real_type> real_set_type;
        typedef ublas::hermitian_matrix< scalar_type, ublas::lower > herm_matrix_type;
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        herm_matrix_type M;
        matrix_type S;

        real_set_type Theta;

        vector_set_type V, W;
        vector_type t = new_vector(vspace_);

        X_.resize(X_.size()+1);

        generate(t, gen);
        project(t, vspace_);

        size_t m = 0, k = 0;

        while(true) //k<k_max
        {
            W.resize(m+1);
            ietl::mult(A_, t, W.back() );
            W.back() -= tau * t;

            for(size_t i=0; i < m; ++i)
            {
                scalar_type gamma = ietl::dot(W[i], W.back());
                W.back() -= gamma * W[i];
                t -= gamma * V[i];
            }

            ++m;
            scalar_type norm_w = ietl::two_norm(W.back());
            W.back() /= norm_w;
            V.push_back( t/norm_w );

            M.resize(m,m);
            for(size_t i=0; i < m; ++i)
                M(i,m-1) = ietl::dot(W[i],V.back());
            
            //eigendecomposition M = S THETA S*
            Theta.resize(m);
            S = M; //copy M

            int info = lapack::heev('V', S, Theta); //eigenvectors in columns of S
                if(info < 0) throw std::runtime_error("lapack::heev - illegal value.");
                if(info > 0) throw std::runtime_error("lapack::heev - failed to converge.");
            //ascending instead of descending

            vector_type u_tilde = V[0] * S(0,m-1);
            for(size_t i=1; i < m; ++i)
                u_tilde += V[i] * S(i,m-1);

            magnitude_type mu = ietl::two_norm(u_tilde);

            X_.back() = u_tilde/mu;
            real_type theta = Theta.back()/mu/mu;

            vector_type w_tilde = W[0] * S(0,m-1);
            for(size_t i=1; i < m; ++i)
                w_tilde += W[i] * S(i,m-1);

            vector_type r = w_tilde/mu - theta * X_.back();

            magnitude_type norm_r = ietl::two_norm(r);

            while(iter.finished(norm_r,Theta[0])) 
            {
                if(iter.error_code() == 1)
                    throw std::runtime_error(iter.error_msg());

                if(verbose_ > 0)
                    std::cout << "Accepted eigenpair #"<< k+1 << "\t" << theta + tau << "\tResidual: " << norm_r << "\n";
                Lambda_.push_back( theta + tau );

                if(++k == k_max) return;

                --m;
                M = ublas::zero_matrix<scalar_type>(m,m);

                vector_set_type tempvs (m);
                for(size_t i = 0; i < m; ++i)
                {
                    tempvs[i] = V[0] * S(0,S.size2()-2-i);
                    for(size_t j = 1; j < S.size1(); ++j)
                        tempvs[i] += V[j] * S(j,S.size2()-2-i);
                }

                V.pop_back();
                for(size_t i = 0; i < m; ++i)
                    std::swap(V[i],tempvs[i]);
                for(size_t i = 0; i < m; ++i)
                {
                    tempvs[i] = W[0] * S(0,S.size2()-2-i);
                    for(size_t j = 1; j < S.size1(); ++j)
                        tempvs[i] += W[j] * S(j,S.size2()-2-i);
                }

                W.pop_back();
                for(size_t i = 0; i < m; ++i)
                    std::swap(W[i],tempvs[i]);
                Theta.pop_back();
                S.resize(m,m);
                for(size_t i = 0; i < m; ++i)
                {
                    M(i,i) = Theta[m-1-i];
                    ublas::column(S, m-1-i) = ublas::unit_vector<scalar_type> (S.size1(), i);
                }

                mu = ietl::two_norm(V[0]);
                theta = Theta.back() / mu / mu;
                X_.push_back(V[0]/mu);

                r = W[0] / mu - theta * X_.back();

            }//accept
            // restart
            if(m >= iter.m_max()) {
                if(verbose_ > 1) 
                    std::cout<<"restart...\n";
                M = ublas::zero_matrix<scalar_type>(iter.m_min(),iter.m_min());

                vector_set_type tempvs (iter.m_min());
                for(size_t i = 1; i < iter.m_min(); ++i)
                {
                    tempvs[i] = V[0] * S(0,m-1-i);
                    for(size_t j = 1; j < m; ++j)
                        tempvs[i] += V[j] * S(j,m-1-i);
                }
                V.resize(iter.m_min());

                for(size_t i = 1; i < iter.m_min(); ++i)
                    std::swap(V[i],tempvs[i]);

                for(size_t i = 1; i < iter.m_min(); ++i)
                {
                    tempvs[i] = W[0] * S(0,m-1-i);
                    for(size_t j = 1; j < m; ++j)
                        tempvs[i] += W[j] * S(j,m-1-i);
                }
                W.resize(iter.m_min());
                for(size_t i = 1; i < iter.m_min(); ++i)
                    std::swap(W[i],tempvs[i]);

                Theta.resize(iter.m_min());
                for(size_t i = 0; i < iter.m_min(); ++i)
                    M(i,i) = Theta[m-1-i];

                W[0] = w_tilde;
                V[0] = u_tilde;
                m = iter.m_min();

            }// restart

            theta += tau;

            // correction equation
            solver(theta, r, t, iter); //solver is allowed to change r

            ++iter;
            if(iter.error_code() == 1)
                throw std::runtime_error(iter.error_msg());
            if(verbose_ > 1)
                std::cout << "JD iteration " << iter.iterations() << "\tresid = " << norm_r << "\n";         
        } //main loop

    }// eigensystem exterior

}// end namespace ietl
#endif
