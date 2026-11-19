/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2003 by Brigitte Surer
 *                       and Jan Gukelberger
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <alps/alea/observables.hpp>
#include <alps/alea/observables_io.hpp>
#include <alps/hdf5.hpp>
#include <boost/random.hpp>
#include <boost/multi_array.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>



class Simulation
{
public:
    Simulation(double beta,size_t L, std::string output_file)
    :   eng_(42)
    ,   rng_(eng_, dist_)
    ,   L_(L)
    ,   beta_(beta)
    ,   spins_(boost::extents[L][L])
    ,   filename_(output_file)
    {
        // Each observable keeps a binning analysis for its error and bins
        // for functions of results.
        measurements_ << alps::alea::observable("E")
                      << alps::alea::observable("m")
                      << alps::alea::observable("|m|")
                      << alps::alea::observable("m^2")
                      << alps::alea::observable("m^4");

        // Init exponential map
        for(int E = -4; E <= 4; E += 2)
            exp_table_[E] = exp(2*beta*E);

        // Init random spin configuration
        for(size_t i = 0; i < L; ++i)
        {
            for(size_t j = 0; j < L; ++j)
                spins_[i][j] = 2 * randint(2) - 1;
        }
    }

    void run(size_t ntherm,size_t n)
    {
        thermalization_ = ntherm;
        sweeps_=n;
        // Thermalize for ntherm steps
        while(ntherm--)
            step();

        // Run n steps
        while(n--)
        {
            step();
            measure();
        }
        
        //save the observables 
        save(filename_);
        
        // Print observables: mean +/- error, the autocorrelation time, and a
        // warning if the error has not converged
        alps::alea::result_set const results = measurements_.results();
        for (std::string const name : {"|m|", "E", "m"})
            std::cout << name << ":\t" << results[name] << std::endl;
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
        measurements_["E"] << E/double(L_*L_);
        double m = M/double(L_*L_);
        measurements_["m"] << m;
        measurements_["|m|"] << fabs(M)/double(L_*L_);
        measurements_["m^2"] << m*m;
        measurements_["m^4"] << m*m*m*m;
    }
    
    void save(std::string const & filename){
        alps::hdf5::archive ar(filename, "a");
        alps::alea::save_results(ar, "/simulation/results", measurements_.results());
        ar["/parameters/L"] << L_;
        ar["/parameters/BETA"] << beta_;
        ar["/parameters/SWEEPS"] << sweeps_;
        ar["/parameters/THERMALIZATION"] << thermalization_;
    }
    
    protected:
    // Random int from the interval [0,max)
    int randint(int max) const
    {
        return static_cast<int>(max * rng_());
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
    size_t sweeps_;
    size_t thermalization_;
    boost::multi_array<int,2> spins_;
    std::map< int, double > exp_table_;

    alps::alea::observable_set measurements_;
    
    std::string filename_;
};


int main(int,char**)
{
    size_t L = 16;      // Linear lattice size
    size_t N = 5000;    // # of simulation steps

    std::cout << "# L: " << L << " N: " << N << std::endl;

    // Scan beta range [0,1] in steps of 0.1
    for(double beta = 0.; beta <= 1.; beta += .1)
    {
        std::cout << "----------" << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::stringstream output_name;
        output_name << "ising.L_" << L << "beta_" << beta <<".h5";
        Simulation sim(beta,L, output_name.str());
        sim.run(N/2,N);
    }

    return 0;
}
