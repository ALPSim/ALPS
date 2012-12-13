// Copyright (C) 20012 TODO: Licence needet!

#include <isingmpi.hpp>

ising_mpi_sim::ising_mpi_sim(parameters_type const & parameters, boost::mpi::communicator const & comm, double t_min = 1, double t_max = 600)
    : ising_sim(parameters)
    , communicator(comm)
    , schedule(t_min, t_max)
    , clone(comm.rank())
    , binnumber(parameters["BINNUMBER"] | std::min(128, 2 * comm.size()))
{}

double fraction_completed() const {
    return fraction;
}

bool run(boost::function<bool ()> const & stop_callback) {
    bool done = false, stopped = false;
    do {
        update();
        measure();
        if ((stopped = stop_callback()) || schedule.pending()) {
            double local_fraction = stopped ? 1. : ising_sim::fraction_completed();
            schedule.update(fraction = boost::mpi::all_reduce(communicator, local_fraction, std::plus<double>()));
            done = fraction >= 1.;
        }
    } while(!done);
    return !stopped;
}

results_type collect_results(result_names_type const & names) const {
    results_type partial_results;
    for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it) {
        alps::mcresult result(measurements[*it]);
        if (result.count())
            partial_results.insert(*it, result.reduce(communicator, binnumber));
        else
            partial_results.insert(*it, result);
    }
    return partial_results;
}
