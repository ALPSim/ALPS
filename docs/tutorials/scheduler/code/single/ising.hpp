// Copyright (C) 20012 TODO: Licence needet!

#ifndef ISING_HPP
#define ISING_HPP

// TODO: remove many ...
#include <alps/ngs/api.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/signal.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/mcresults.hpp>
#include <alps/ngs/mcobservables.hpp>

// #ifdef ALPS_NGS_USE_NEW_ALEA
//     #include <alps/ngs/alea/accumulator_set.hpp>
// #endif

#include <alps/ngs/observablewrappers.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <alps/random/mersenne_twister.hpp>

#include <boost/chrono.hpp>
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

class ising_sim {
    public:

        typedef alps::params parameters_type;
        typedef alps::mcresults results_type;
        typedef std::vector<std::string> result_names_type;

        ising_sim(parameters_type const & parameters);

        void update();
        void measure();
    
        double fraction_completed() const;

        void save(boost::filesystem::path const & filename) const;
        void load(boost::filesystem::path const & filename);

        // functions needet to save/load the simulation to/from hdf5
        void save(alps::hdf5::archive & ar) const;
        void load(alps::hdf5::archive & ar);

        bool run(boost::function<bool ()> const & stop_callback);

        result_names_type result_names() const;
        result_names_type unsaved_result_names() const;

        results_type collect_results() const;
        results_type collect_results(result_names_type const & names) const;

    protected:

        parameters_type params;
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > mutable random;
        alps::mcobservables measurements;
    
        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        std::vector<int> spins;
};

#endif
