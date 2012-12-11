// Copyright (C) 20012 TODO: Licence needet!

#include <ising.hpp>

ising_sim::ising_sim(parameters_type const & parameters)
    : params(parameters)
    , random(boost::mt19937((parameters["SEED"] | 42)), boost::uniform_real<>())
    , length(parameters["L"])
    , sweeps(0)
    , thermalization_sweeps(int(parameters["THERMALIZATION"]))
    , total_sweeps(int(parameters["SWEEPS"]))
    , beta(1. / double(parameters["T"]))
    , spins(length)
{
    for(int i = 0; i < length; ++i)
        spins[i] = (random() < 0.5 ? 1 : -1);
    measurements
        << alps::ngs::RealObservable("Energy")
        << alps::ngs::RealObservable("Magnetization")
        << alps::ngs::RealObservable("Magnetization^2")
        << alps::ngs::RealObservable("Magnetization^4")
        << alps::ngs::RealVectorObservable("Correlations")
    ;
}

void ising_sim::update() {
    for (int j = 0; j < length; ++j) {
        using std::exp;
        int i = int(double(length) * random());
        int right = ( i + 1 < length ? i + 1 : 0 );
        int left = ( i - 1 < 0 ? length - 1 : i - 1 );
        double p = exp( 2. * beta * spins[i] * ( spins[right] + spins[left] ));
        if ( p >= 1. || random() < p )
            spins[i] = -spins[i];
    }
};

void ising_sim::measure() {
    sweeps++;
    if (sweeps > thermalization_sweeps) {
        double tmag = 0;
        double ten = 0;
        double sign = 1;
        std::vector<double> corr(length);
        for (int i = 0; i < length; ++i) {
            tmag += spins[i];
            sign *= spins[i];
            ten += -spins[i] * spins[ i + 1 < length ? i + 1 : 0 ];
            for (int d = 0; d < length; ++d)
                corr[d] += spins[i] * spins[( i + d ) % length ];
        }
        std::transform(corr.begin(), corr.end(), corr.begin(), boost::lambda::_1 / double(length));
        ten /= length;
        tmag /= length;
        measurements["Energy"] << ten;
        measurements["Magnetization"] << tmag;
        measurements["Magnetization^2"] << tmag * tmag;
        measurements["Magnetization^4"] << tmag * tmag * tmag * tmag;
        measurements["Correlations"] << corr;
    }
};

double fraction_completed() const {
    return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}

void ising_sim::save(boost::filesystem::path const & filename) const {
    alps::hdf5::archive ar(filename, "w");
    ar << *this;
}

void ising_sim::load(boost::filesystem::path const & filename) {
    alps::hdf5::archive ar(filename);
    ar >> *this;
}

void ising_sim::save(alps::hdf5::archive & ar) const {
    ar["/parameters"] << params;
    
    // Set the current path of the archive to /simulation/realizations/0/clones/0
    // all relative path will be saved according to this path
    std::string context = ar.get_context();
    ar.set_context("/simulation/realizations/0/clones/0");
    
    ar["measurements"] << measurements;
    
    ar.set_context("checkpoint");
    ar["length"] << length;
    ar["sweeps"] << sweeps;
    ar["thermalization_sweeps"] << thermalization_sweeps;
    ar["beta"] << beta;
    ar["spins"] << spins;
    
    // also save the state of the random number generator to avoid overlapping sequences
    {
        std::ostringstream os;
        os << random.engine();
        ar["engine"] << os.str();
    }

    // put the archive back to the original state
    ar.set_context(context);
}

void ising_sim::load(alps::hdf5::archive & ar) {
    ar["/parameters"] >> params;

    std::string context = ar.get_context();
    ar.set_context("/simulation/realizations/0/clones/0");
    ar["measurements"] >> measurements;

    ar.set_context("checkpoint");
    ar["length"] >> length;
    ar["sweeps"] >> sweeps;
    ar["thermalization_sweeps"] >> thermalization_sweeps;
    ar["beta"] >> beta;
    ar["spins"] >> spins;

    {
        std::string state;
        ar["engine"] >> state;
        std::istringstream is(state);
        is >> random.engine();
    }

    ar.set_context(context);
}

bool ising_sim::run(boost::function<bool ()> const & stop_callback) {
    bool stopped = false;
    do {
        update();
        measure();
    } while(!(stopped = stop_callback()) && fraction_completed() < 1.);
    return !stopped;
}

ising_sim::result_names_type ising_sim::result_names() const {
    result_names_type names;
    for(observables_type::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
        names.push_back(it->first);
    return names;
}

ising_sim::result_names_type ising_sim::unsaved_result_names() const {
    return result_names_type(); 
}

ising_sim::results_type ising_sim::collect_results() const {
    return collect_results(result_names());
}

ising_sim::results_type ising_sim::collect_results(result_names_type const & names) const {
    results_type partial_results;
    for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it)
        partial_results.insert(*it, alps::mcresult(measurements[*it]));
    return partial_results;
}
