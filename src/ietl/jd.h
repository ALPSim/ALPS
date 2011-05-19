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

        jd_iteration(size_t m_min, size_t m_max, size_t max_iter, T reltol = 0., T abstol = 0., size_t smi = 5) 
            : basic_iteration<T>(max_iter, reltol, abstol), m_min_(m_min), m_max_(m_max), smi_(smi) {}

        inline size_t m_min() const
            {   return m_min_;  }
        inline size_t m_max() const
            {   return m_max_;  }
        //maximal number iterations for the correction equation solver
        inline size_t sub_max_iter() const
            {   return smi_;  }
        const std::string& error_msg() 
            { return basic_iteration<T>::err_msg; }    
        private:
        size_t m_min_, m_max_, smi_;
    };
//jd_iterator//////////////////////////////////////////////////////
    namespace detail{

//deflated matrix multiplication (I-QQ*)(A-\theta I)(I-QQ*)////////

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

}//end namespace detail///////////////////////////////////////////////////////
//solver//////////////////////////////////////////////////////////////////////
namespace solver {

//left preconditioned correction-equation-solver
    template <class MATRIX, class VS, class PREC>   class left_prec_solver;

    template <class MATRIX, class VS, class PREC, class VECTOR>
    void mult (const left_prec_solver<MATRIX,VS,PREC>&, const VECTOR&, VECTOR&);

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
        void operator() ( real_type theta, vector_type& r, vector_type& t, IT& iter );

        friend void mult<>(const left_prec_solver& A_def, const vector_type& v, vector_type& z);

    protected:
        const MATRIX& A_;
        const std::vector<vector_type>& Q_;
        const VS& vspace_;
        std::vector<vector_type> Q_hat_;
        PREC& K_;
        vector_type x0_;
        real_type theta_;
        mutable ublas::vector<scalar_type> gamma_;
        std::vector<int> pivot_;
        matrix_type LU_, M_;
    };

    template <class MATRIX, class VS, class PREC, class VECTOR>
    void mult (const left_prec_solver<MATRIX,VS,PREC>& LPS, const VECTOR& v, VECTOR& z)
    {
        VECTOR temp;

        z = LPS.theta_ * v;

        ietl::mult(LPS.A_, v, temp);

        temp -= z;
        // K y_tilde = y
        ietl::mult(LPS.K_, temp, z);

        for(size_t i = 0; i < LPS.Q_.size(); ++i)
            LPS.gamma_(i) = ietl::dot(LPS.Q_[i],z);

        int info = lapack::getrs(LPS.LU_, LPS.pivot_, LPS.gamma_);
            if(info != 0)   throw std::runtime_error("lapack::getrs failed.");

        for(size_t i = 0; i < LPS.Q_hat_.size(); ++i)
            z -= LPS.Q_hat_[i]* LPS.gamma_(i);
    }

    // jacobi preconditioner for (A - theta *I)
    template <class MATRIX, class VS, class PREC> template <class IT>
    void left_prec_solver<MATRIX,VS,PREC>::operator() ( real_type theta, vector_type& r, vector_type& t, IT& iter )
    {
        theta_ = theta;
        const unsigned int m = Q_.size(), m_old = Q_hat_.size();//if K changes =0;
        const size_t max_iter = iter.sub_max_iter();
        double abs_tol = iter.absolute_tolerance();
        abs_tol *= abs_tol;

        M_.resize(m,m);
        gamma_.resize(m);
        pivot_.resize(m);
        Q_hat_.resize(m);

        for(size_t i = m_old; i < m; ++i)
            ietl::mult(K_, Q_[i], Q_hat_[i]);

        for(size_t i = 0; i < m; ++i)
            for(size_t j = ( (i<m_old) ? m_old : 0 ); j < m; ++j)
                M_(i,j) = ietl::dot(Q_[i], Q_hat_[j]);

        LU_ = M_;

        int info = lapack::getrf(LU_, pivot_);
            if(info != 0)   throw std::runtime_error("lapack::getrf failed.");

        ietl::mult(K_, r, t);

        for(size_t i = 0; i < m; ++i)
            gamma_(i) = ietl::dot(Q_[i],t);

        //calculate alpha from M alpha = gamma
        info = lapack::getrs(LU_, pivot_, gamma_);
            if(info != 0)   throw std::runtime_error("lapack::getrs failed.");

        r = -t;

        for(size_t i = 0; i < m; ++i)
            r += gamma_(i)*Q_hat_[i];

        t = ietl_gmres(*this, r, x0_, max_iter, abs_tol, false);

    }//left_prec_solver::void()

//simple solver with gmres//////////////////////////////////////////////////////
    template <class MATRIX, class VS>
    class jd_solver {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename real_type<scalar_type>::type real_type;

        jd_solver (const MATRIX& A, const std::vector<vector_type>& Q, const VS& vspace)
            : Adef_(A, 0, Q), vspace_(vspace), Q_(Q), x0_(new_vector(vspace))
             {}
        
        template <class IT>
        void operator() ( scalar_type theta, vector_type& r, vector_type& t, IT& iter );

    protected:
        detail::deflated_matrix<MATRIX, VS> Adef_;
        const VS& vspace_;
        vector_type x0_;
        const std::vector<vector_type>& Q_;
    };
        template <class MATRIX, class VS>
        template <class IT>
        void jd_solver<MATRIX, VS>::operator()( scalar_type theta, vector_type& r, vector_type& t, IT& iter)
        {
            Adef_.set_theta(theta);

            //one step approximation
            /*x0_ = -r;
            for(size_t i = 0; i < Q_.size(); ++i)
                x0_ += ietl::dot(Q_[i],r)/ietl::dot(Q_[i],Q_[i]) * Q_[i];*/

            r *= -1;

            //starting with x0_ = 0
            t = ietl_gmres(Adef_, r, x0_, iter.sub_max_iter(), iter.absolute_tolerance()*iter.absolute_tolerance(), false);
        }
}//end namespace solver///////////////////////////////////////////////////////////

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

        //without preconditioning
        template <class IT, class GEN>
        void eigensystem(IT& iter, GEN& gen, size_type k_max, bool search_highest = false)    
            {
                typedef solver::jd_solver<MATRIX,VS> solver_type;
                solver_type solver(A_, X_, vspace_);
                eigensystem <solver_type,IT,GEN> (solver, iter, gen, k_max, search_highest);
            }

        //with left preconditioner K with K (A - \theta I) ~= I
        template <class IT, class GEN, class PREC>
        void eigensystem( IT& iter, GEN& gen, size_type k_max, PREC& K, bool search_highest = false)
            {
                typedef solver::left_prec_solver<MATRIX,VS,PREC> solver_type;
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
                solver_type solver(A_, X_, vspace_);
                eigensystem_harmonic <solver_type,IT,GEN> (solver, iter, gen, k_max, tau);
            }

        //access functions
        std::pair<real_type,vector_type> eigenpair(size_type k)
            {
                assert (k < Lambda_.size());
                return std::make_pair (Lambda_[k], X_[k]);
            }
        real_type eigenvalue(size_type k)
            {
                assert (k < Lambda_.size());
                return Lambda_[k];
            }
        real_set_type eigenvalues()
            {
                return Lambda_;
            }
        vector_type eigenvector (size_type k)
            {
                assert (k < X_.size());
                return X_[k];
            }
        void reset()
            {
                X_.clear();
                Lambda_.clear();
            }

    protected:
        MATRIX& A_;
        VS& vspace_;
        size_type n_;
        vector_set_type X_;
        real_set_type Lambda_;
        size_t verbose_;    //verbosity levels: 0 = no, 1 = yes, 2 = all
    };// end of class jd

//jd for exterior eigenvalues ///////////////////////////////////////////////////
    template <class MATRIX, class VS>
    template <class SOLVER, class IT, class GEN>
    void jd<MATRIX,VS>::eigensystem(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, bool search_highest)
    {
        typedef std::vector<vector_type> vector_set_type;
        typedef std::vector<real_type> real_set_type;
        typedef ublas::hermitian_matrix< scalar_type, ublas::upper > herm_matrix_type;
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        assert(k_max <= n_);

        const magnitude_type kappa = 0.25; // kappa < 1
        magnitude_type norm_t, norm_r;
        const size_type m_min = iter.m_min(), m_max = iter.m_max();
        size_type k = 0, m = 0;
        int info;

        vector_set_type V;
        vector_set_type VA;
        X_.resize( X_.size()+1 ); //X_.back() == u

        herm_matrix_type M;

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

            for(size_type i = 0; i < m; ++i) 
                t -= ietl::dot(V[i],t) * V[i];
            if(ietl::two_norm(t) <= kappa * norm_t)
                for(size_type i = 0; i < m; ++i)
                    t -= ietl::dot(V[i],t) * V[i];

            norm_t = ietl::two_norm(t);

            if(norm_t == 0/*< iter.absolute_tolerance()*/) {
                    std::cerr << "Correction vector t is zero:\n";
                    if(iter.first())
                        throw std::runtime_error("generating a starting vector failed.");
                    else
                        throw std::runtime_error("try to solve the correction equation more accurately.");
                    }

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
            norm_r = ietl::two_norm(r);

            while(iter.finished(norm_r,theta)) 
            {

                if(iter.error_code() == 1)
                    throw std::runtime_error(iter.error_msg());

                if(verbose_ > 0)
                    std::cout << "Accepted eigenpair #"<< k+1 << "\t" << theta << "\tResidual: " << norm_r << "\n";

                Lambda_.push_back( theta ); //Lambda_[k] = Theta[0];

                if(++k == k_max) return;

                if(m < 1) {//could start anew with a random vector
                    throw std::runtime_error("search space is depleted, try a different generator.");
                }
                --m;
                M = ublas::zero_matrix<scalar_type>(m,m);

                //use VA temporary as V
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

            //assure t is orthogonal to Q
            norm_t = ietl::two_norm(t);

            for(size_t i = 0; i< X_.size();++i)
                t -= ietl::dot(t,X_[i])*X_[i];

            if(ietl::two_norm(t) <= kappa * norm_t)
                for(size_t i = 0; i< X_.size();++i)
                    t -= ietl::dot(t,X_[i])*X_[i];

            ++iter;
            if(iter.error_code() == 1)
                throw std::runtime_error(iter.error_msg());
            if(verbose_ > 1)
                std::cout << "JD iteration " << iter.iterations() << "\t resid = " << norm_r << "\n";
        }// main loop

    }// eigensystem

//jd for interior eigenvalues with harmonic ritz values//////////////////////////////////////////////////

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
            if(norm_w == 0)
                throw std::runtime_error("New search vector is zero.");

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
                    std::cout << "Accepted eigenpair #"<< k+1 << "\t" << theta + tau << "\t Residual: " << norm_r << "\n";
                Lambda_.push_back( theta + tau );

                if(++k == k_max) return;

                if(m < 1)
                    throw std::runtime_error("search space is depleted, try a different generator.");
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

            for(size_t i = 0; i < X_.size(); ++i)
                t -= ietl::dot(t,X_[i])*X_[i];

            ++iter;
            if(iter.error_code() == 1)
                throw std::runtime_error(iter.error_msg());
            if(verbose_ > 1)
                std::cout << "JD iteration " << iter.iterations() << "\tresid = " << norm_r << "\n";         
        } //main loop

    }// eigensystem exterior

}// end namespace ietl
#endif
