/* ndim_spin.hpp
 * adapted from heisenberg_lattice/heisenberg.hpp
 */

#ifndef HEISENBERG_HPP
#define HEISENBERG_HPP

#include "tinyvector/tinyvector.hpp"

#include <alps/mcbase.hpp>
#include <alps/ngs/numeric.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/random/uniform_on_sphere_n.h>
#include <alps/lattice.h>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/array.hpp>

#include <boost/function.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>

#include <vector>
#include <string>

/**
 * The simulation class derives from mcbase.
 */
template <int N>
class ALPS_DECL ndim_spin_sim : public alps::mcbase {

    public:
        typedef tinyvector<double, N, INTRIN_OPT> spintype;
        ndim_spin_sim(parameters_type const & parms, std::size_t seed_offset = 0);

        virtual void update();
        virtual void measure();
        virtual double fraction_completed() const;

        using alps::mcbase::save;
        virtual void save(alps::hdf5::archive & ar) const;

        using alps::mcbase::load;
        virtual void load(alps::hdf5::archive & ar);

        // convenience function to get a random spin (uniformly distancesributed direction)
        const spintype random_spin();

    private:
        
        alps::graph_helper<> lattice;
        int num_sites;
        std::vector<spintype> spins;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        alps::uniform_on_sphere_n<N, double, spintype > random_spin_gen;
        std::vector<double> distances;
};

/* implementation */

template<int N>
ndim_spin_sim<N>::ndim_spin_sim(parameters_type const & parms, std::size_t seed_offset)
    : alps::mcbase(parms, seed_offset)
    , lattice(alps::make_deprecated_parameters(parms))
    , num_sites(lattice.num_sites())
    , spins(num_sites)
    , sweeps(0)
    , thermalization_sweeps(int(parameters["THERMALIZATION"]))
    , total_sweeps(int(parameters["SWEEPS"]))
    , beta(1. / double(parameters["T"]))
    , random_spin_gen()
{
    // initialize spins with random values
    for(int i = 0; i < num_sites; ++i) {
        spins[i] = random_spin();
    }

    measurements
        << alps::accumulator::RealObservable("Energy")
        << alps::accumulator::RealVectorObservable("Magnetization")
        << alps::accumulator::RealObservable("Magnetization^2")
        << alps::accumulator::RealObservable("Magnetization^4")
        << alps::accumulator::RealVectorObservable("Correlations")
        << alps::accumulator::RealVectorObservable("Distances")
    ;

    using alps::ngs::numeric::operator-;
    std::vector<double> ref = lattice.coordinate(0);
    std::vector<double> a;
    for (int i = 0; i < num_sites; ++i) {
        double d = 0;
        a = lattice.coordinate(i) - ref;
        for (int j = 0; j < a.size(); ++j) {
            d += a[j] * a[j];
        }
        distances.push_back(std::sqrt(d));
    }
    measurements["Distances"] << distances;
}

template<int N>
void ndim_spin_sim<N>::update() {
    for (int j = 0; j < num_sites; ++j) {
        using std::exp;
        int i = int(double(num_sites) * random());

        // get a random site
        alps::graph_helper<>::site_descriptor site_i = lattice.site(i);

        // sum over all it's neighbors for delta_H later
        typename ndim_spin_sim<N>::spintype nn_sum;
        nn_sum.initialize(0);
        alps::graph_helper<>::neighbor_iterator nn_it, nn_end;
        for(boost::tie(nn_it, nn_end) = lattice.neighbors(site_i); nn_it != nn_end; ++nn_it) {
            nn_sum += spins[*nn_it];
        }

        // generate a new random spin and decide if we keep it
        typename ndim_spin_sim<N>::spintype new_spin = random_spin();
        double delta_H = dot(new_spin - spins[i], nn_sum);
        double p = exp( -beta * delta_H );
        if ( p >= 1. || random() < p )
            spins[i] = new_spin;
    }
}

template <int N>
void ndim_spin_sim<N>::measure() {
    sweeps++;
    if (sweeps > thermalization_sweeps) {
        typename ndim_spin_sim<N>::spintype magnetization;
        magnetization.initialize(0);
        double energy = 0;
        std::vector<double> correlations(num_sites, 0);
        // To measure magnetization, magnetic susceptibility and correlationselations we sum over all sites.
        for (int i = 0; i < num_sites; ++i) {
            magnetization += spins[i];
            correlations[i] = dot(spins[0], spins[i]);
        }
        // To measure the Energy we sum only over neighbored sites (bonds in terms of the lattice).
        alps::graph_helper<>::bond_iterator bond_it, bond_end;
        for(boost::tie(bond_it, bond_end) = lattice.bonds(); bond_it != bond_end; ++bond_it) {
            energy += - dot(spins[lattice.source(*bond_it)], spins[lattice.target(*bond_it)]);
        }
        // pull in operator/ for vectors
        using alps::ngs::numeric::operator/;
        energy /= num_sites;                // $\frac{1}{V} \sum_{\text{i,j nn}}{\sigma_i \sigma_j}$
        magnetization /= num_sites;         // $\frac{1}{V} \sum_{i}{\sigma_i}$
        double magnetization2 = dot(magnetization, magnetization);

        // store the measurements
        measurements["Energy"] << energy;
        measurements["Magnetization"] << spintype::vector(magnetization);
        measurements["Magnetization^2"] << magnetization2;
        measurements["Magnetization^4"] << magnetization2 * magnetization2;
        measurements["Correlations"] << correlations;
    }
}

template <int N>
double ndim_spin_sim<N>::fraction_completed() const {
    return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}

template <int N>
void ndim_spin_sim<N>::save(alps::hdf5::archive & ar) const {
    mcbase::save(ar);
    ar["checkpoint/sweeps"] << sweeps;
    ar["checkpoint/spins"] << spins;
}

template <int N>
void ndim_spin_sim<N>::load(alps::hdf5::archive & ar) {
    mcbase::load(ar);
    ar["checkpoint/sweeps"] >> sweeps;
    ar["checkpoint/spins"] >> spins;
}

template <int N>
const typename ndim_spin_sim<N>::spintype ndim_spin_sim<N>::random_spin() {
    return random_spin_gen(random.engine());
}

#endif 
