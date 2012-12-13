// Copyright (C) 20012 TODO: Licence needet!

#ifndef ISING_MPI_HPP
#define ISING_MPI_HPP

#include "ising.hpp"

#include <alps/ngs/boost_mpi.hpp>
#include <alps/ngs/scheduler/check_schedule.hpp>

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>

class ising_mpi_sim : public ising_sim {

    public:

        ising_mpi_sim(parameters_type const & parameters, boost::mpi::communicator const & comm, double t_min = 1, double t_max = 600)

        double fraction_completed() const;

        bool run(boost::function<bool ()> const & stop_callback);

        results_type collect_results(result_names_type const & names) const;

    protected:

        boost::mpi::communicator communicator;

        alps::check_schedule schedule;
        double fraction;
        int clone;
        std::size_t binnumber;
};

#endif
