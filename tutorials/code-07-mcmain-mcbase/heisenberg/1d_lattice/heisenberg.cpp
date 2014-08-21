/* heisenberg.cpp
 * adapted from alps/tutorials/code-07-mcmain-mcbase/heisenberg.cpp
 * 1d grid, 3d spins
 */

#include "heisenberg.hpp"
#include "helper.hpp"

#include <boost/lambda/lambda.hpp>

heisenberg_sim::heisenberg_sim(parameters_type const & parms, std::size_t seed_offset)
    : alps::mcbase(parms, seed_offset)
    , length(parameters["L"])
    , sweeps(0)
    , thermalization_sweeps(int(parameters["THERMALIZATION"]))
    , total_sweeps(int(parameters["SWEEPS"]))
    , beta(1. / double(parameters["T"]))
    , random_spin_gen()
    , spins(length)
{
    /* initialize spins */
    for(int i = 0; i < length; ++i)
    {
        spins[i] = random_spin();
    }
    measurements
        << alps::accumulator::RealObservable("Energy")
        << alps::accumulator::RealVectorObservable("Magnetization")
        << alps::accumulator::RealObservable("Magnetization^2")
        << alps::accumulator::RealObservable("Magnetization^4")
        << alps::accumulator::RealVectorObservable("Correlations")
    ;
}

void heisenberg_sim::update() {
    for (int j = 0; j < length; ++j) {
        using std::exp;
        int i = int(double(length) * random());
        int right = ( i + 1 < length ? i + 1 : 0 );
        int left = ( i - 1 < 0 ? length - 1 : i - 1 );
        /* generate a new random spin for update and accept or not */
        spintype new_spin = random_spin();
        double delta_H = dot(new_spin - spins[i], spins[right] + spins[left]);
        double p = exp( -beta * delta_H );
        if ( p >= 1. || random() < p )
            spins[i] = new_spin;
    }
}

void heisenberg_sim::measure() {
    sweeps++;
    if (sweeps > thermalization_sweeps) {
        spintype tmag = {{0, 0, 0}};
        double ten = 0;
        //~ double sign = 1;
        std::vector<double> corr(length);
        for (int i = 0; i < length; ++i) {
            tmag += spins[i];
            //~ sign *= spins[i];
            ten += dot(-spins[i], spins[ i + 1 < length ? i + 1 : 0 ]);
            for (int d = 0; d < length; ++d)
                corr[d] += dot(spins[i], spins[( i + d ) % length ]);
        }
        // pull in operator/ for vectors
        using alps::ngs::numeric::operator/;
        corr = corr / double(length);
        ten /= length;
        tmag /= length;
        measurements["Energy"] << ten;
        measurements["Magnetization"] << vector_from_spintype(tmag);
        measurements["Magnetization^2"] << dot(tmag, tmag);
        measurements["Magnetization^4"] << dot(tmag, tmag) * dot(tmag, tmag);
        measurements["Correlations"] << corr;
    }
}

double heisenberg_sim::fraction_completed() const {
    return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}

void heisenberg_sim::save(alps::hdf5::archive & ar) const {
    mcbase::save(ar);
    ar["checkpoint/sweeps"] << sweeps;
    ar["checkpoint/spins"] << spins;
}

void heisenberg_sim::load(alps::hdf5::archive & ar) {
    mcbase::load(ar);

    length = int(parameters["L"]);
    thermalization_sweeps = int(parameters["THERMALIZATION"]);
    total_sweeps = int(parameters["SWEEPS"]);
    beta = 1. / double(parameters["T"]);

    ar["checkpoint/sweeps"] >> sweeps;
    ar["checkpoint/spins"] >> spins;
}

const spintype heisenberg_sim::random_spin() {
    return spin_from_vector(random_spin_gen(random.engine()));
}

const spintype heisenberg_sim::random_spin_accept_reject() {
    double r = 0;
    spintype spin = { {0, 0, 0} };
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
