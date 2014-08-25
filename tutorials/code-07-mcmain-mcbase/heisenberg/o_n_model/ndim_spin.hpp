/* ndim_spin.hpp
 * adapted from heisenberg_lattice/heisenberg.hpp
 */

#ifndef HEISENBERG_HPP
#define HEISENBERG_HPP

#include "spintype.hpp"

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

template <int N>
class ALPS_DECL ndim_spin_sim : public alps::mcbase {

    public:
        typedef spintype<N> spintype;
        ndim_spin_sim(parameters_type const & parms, std::size_t seed_offset = 0);

        virtual void update();
        virtual void measure();
        virtual double fraction_completed() const;

        const spintype random_spin_accept_reject();
        const spintype random_spin();

        using alps::mcbase::save;
        virtual void save(alps::hdf5::archive & ar) const;

        using alps::mcbase::load;
        virtual void load(alps::hdf5::archive & ar);

    private:
        
        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        alps::uniform_on_sphere_n<N, double, std::vector<double> > random_spin_gen;
        alps::graph_helper<> lattice;
        std::vector<spintype> spins;
};

/* implementation */

template<int N>
ndim_spin_sim<N>::ndim_spin_sim(parameters_type const & parms, std::size_t seed_offset)
    : alps::mcbase(parms, seed_offset)
    , length(parameters["L"])
    , sweeps(0)
    , thermalization_sweeps(int(parameters["THERMALIZATION"]))
    , total_sweeps(int(parameters["SWEEPS"]))
    , beta(1. / double(parameters["T"]))
    , random_spin_gen()
    , lattice(alps::make_deprecated_parameters(parms))
    , spins(lattice.num_sites())
{
    /* initialize spins */
    alps::graph_helper<>::site_iterator site_it, site_end;
    for(boost::tie(site_it, site_end) = lattice.sites(); site_it != site_end; ++site_it)
    {
        spins[*site_it] = random_spin();
    }
    measurements
        << alps::accumulator::RealObservable("Energy")
        << alps::accumulator::RealVectorObservable("Magnetization")
        << alps::accumulator::RealObservable("Magnetization^2")
        << alps::accumulator::RealObservable("Magnetization^4")
        << alps::accumulator::RealObservable("Magnetic Susceptibility")
    ;
}

//TODO make 3d
template<int N>
void ndim_spin_sim<N>::update() {
    for (int j = 0; j < lattice.num_sites(); ++j) {
        using std::exp;
        int i = int(double(lattice.num_sites()) * random());
        /* get a random site */
        alps::graph_helper<>::site_descriptor site_i = lattice.site(i);
        /* sum all neighbors */
        typename ndim_spin_sim<N>::spintype nn_sum;
        nn_sum.initialize(0);
        alps::graph_helper<>::neighbor_iterator nn_it, nn_end;
        for(boost::tie(nn_it, nn_end) = lattice.neighbors(site_i); nn_it != nn_end; ++nn_it) {
            nn_sum += spins[*nn_it];
        }
        /* generate a new random spin for update and accept or not */
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
        typename ndim_spin_sim<N>::spintype tmag;
        tmag.initialize(0);
        double tmag_sq = 0; // for susceptibility
        double ten = 0;
        for (int i = 0; i < lattice.num_sites(); ++i) {
            tmag += spins[i];
            tmag_sq += dot(tmag, tmag);
        }
        alps::graph_helper<>::bond_iterator bond_it, bond_end;
        for(boost::tie(bond_it, bond_end) = lattice.bonds(); bond_it != bond_end; ++bond_it) {
            ten += - dot(spins[lattice.source(*bond_it)], spins[lattice.target(*bond_it)]);
        }
        // pull in operator/ for vectors
        using alps::ngs::numeric::operator/;
        ten /= lattice.num_sites();
        tmag /= lattice.num_sites();
        double tmag_sq_avg = tmag_sq / lattice.num_sites();
        double tmag_avg_sq = dot(tmag, tmag);
        measurements["Energy"] << ten;
        measurements["Magnetization"] << vector_from_spintype(tmag);
        measurements["Magnetization^2"] << dot(tmag, tmag);
        measurements["Magnetization^4"] << dot(tmag, tmag) * dot(tmag, tmag);
        measurements["Magnetic Susceptibility"] << beta * (tmag_sq_avg - tmag_avg_sq * lattice.num_sites());
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

    length = int(parameters["L"]);
    thermalization_sweeps = int(parameters["THERMALIZATION"]);
    total_sweeps = int(parameters["SWEEPS"]);
    beta = 1. / double(parameters["T"]);

    ar["checkpoint/sweeps"] >> sweeps;
    ar["checkpoint/spins"] >> spins;
}

template <int N>
const typename ndim_spin_sim<N>::spintype ndim_spin_sim<N>::random_spin() {
    return spin_from_vector<N>(random_spin_gen(random.engine()));
}

template <int N>
const typename ndim_spin_sim<N>::spintype ndim_spin_sim<N>::random_spin_accept_reject() {
    double r = 0;
    typename ndim_spin_sim<N>::spintype spin;
    /* get a random vector from uniformly distributed numbers in [-1,1) */
    do 
    {
        for(int i = 0; i < 3; ++i)
        {
            spin[i] = 2. * random() - 1.;
        }
        r = abs(spin);
    }
    while(r >= 1 || r == 0); /* reject and replace if not in unit sphere */
    /* normalize to get a point on the 3d unit sphere with uniformly distributed angles*/
    for(int i = 0; i < 3; ++i)
    {
        spin[i] /= r;
    } 
    return spin;
}

#endif 
