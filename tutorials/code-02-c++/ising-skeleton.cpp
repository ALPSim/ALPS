#include <alps/scheduler/montecarlo.h>
#include <alps/alea.h>
#include <alps/hdf5.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <sstream>


class Simulation
{
public:
    Simulation(double beta,size_t L)
    :    rng_(eng_, dist_)
    ,    L_(L)
    ,   beta_(beta)
    ,    sites_(L,std::vector<int>(L))
    ,    energy_("E")
    ,    magnetization_("m")
    ,    order_("|m|")
    {
        // Init exponential map
        for(int E = -4; E <= 4; E += 2)
            exp_[E] = exp(2*beta*E);
        
        // Init random spin configuration
        for(size_t i = 0; i < L; ++i)
        {
            for(size_t j = 0; j < L; ++j)
                sites_[i][j] = 2 * randint(2) - 1;
        }
    }
    
    void run(size_t ntherm,size_t n)
    {
        n_ = n;
        ntherm_ = ntherm;
        // Thermalize for ntherm steps
        while(ntherm--)
            step();
        
        // Run n steps
        while(n--)
        {
            step();
            measure();
        }
        
        // Print observables
        print_observable(energy_);
        print_observable(order_);
        print_observable(magnetization_);
    }
    void step()
    {
        for(size_t s = 0; s < L_*L_; ++s)
        {
            // Pick random site k=(i,j)
            ...
            
            // Measure local energy e = -s_k * sum_{l nn k} s_l
            ...
            
            // Flip s_k with probability exp(2 beta e)
            ...
        }
    }
    void measure()
    {
        int E = 0; // energy
        int M = 0; // magnetization
        for(size_t i = 0; i < L_; ++i)
        {
            for(size_t j = 0; j < L_; ++j)
            {
                E -= ...
                M += ...
            }
        }
        
        // Add sample to observables
        energy_ << E/double(L_*L_);
        magnetization_ << M/double(L_*L_);
        order_ << fabs(M)/double(L_*L_);
    }
    
    void save(std::string const & filename) const {
        alps::hdf5::oarchive ar(filename);
        ar << alps::make_pvp("/parameters/L", L_);
        ar << alps::make_pvp("/parameters/BETA", beta_);
        ar << alps::make_pvp("/parameters/SWEEPS", n_);
        ar << alps::make_pvp("/parameters/THERMALIZATION", ntherm_);
        ar << alps::make_pvp("/simulation/results/"+energy_.representation(), energy_);
        ar << alps::make_pvp("/simulation/results/"+magnetization_.representation(), magnetization_);
        ar << alps::make_pvp("/simulation/results/"+order_.representation(), order_);
    }
    
    
protected:
    // Random int from the interval [0,max)
    int randint(int max) const
    {
        return static_cast<int>(max * rng_());
    }
    // Enforce periodic boundary conditions for lattice indices
    int wrap(int i) const
    {
        return (i + L_) % L_;
    }
    void print_observable(const alps::RealObservable& o) const
    {
        std::cout << o.name() << ":\t" << o.mean() << " +- " << o.error() << ";\ttau = " << o.tau() 
        << ";\tconverged: " << alps::convergence_to_text(o.converged_errors()) << std::endl;
    }
    
private:
    typedef boost::mt19937 engine_type;
    typedef boost::uniform_real<> distribution_type;
    typedef boost::variate_generator<engine_type&, distribution_type> rng_type;
    engine_type eng_;
    distribution_type dist_;
    mutable rng_type rng_;
    
    size_t L_;
    double beta_;
    std::vector< std::vector<int> > sites_;
    std::map< int, double > exp_;
    alps::RealObservable energy_;
    alps::RealObservable magnetization_;
    alps::RealObservable order_;
    size_t ntherm_;
    size_t n_;
};


int main(int,char**)
{
    size_t L = 8;    // Linear lattice size
    size_t N = 5000;    // # of simulation steps
    
    std::cout << "# L: " << L << " N: " << N << std::endl;
    
    // Scan beta range [0,1] in steps of 0.1
    for(double beta = 0.; beta <= 1.; beta += .1)
    {
        std::cout << "----------" << std::endl;
        std::cout << "beta = " << beta << std::endl;
        Simulation sim(beta,L);
        sim.run(N/2,N);
        std::stringstream fname;
        fname << "ising_L_" << L << "beta_" << beta << ".h5"  ;
        sim.save(fname.str());
    }
    
    return 0;
}
