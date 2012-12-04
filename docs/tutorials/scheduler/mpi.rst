MPI simulation
==============

synopsis
--------

In this tutorial we make ouer simulation mpi capable:

.. code-block:: c++

	class ising_mpi_sim : public ising_sim {

	    public:

	        ising_sim(parameters_type const & parameters, boost::mpi::communicator const & comm, double t_min = 1, double t_max = 600)

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

the ``alps::check_schedule`` class is used to determin when to communicate next. ``alps::check_schedule`` has the following interface:

* ``schedule(t_min, t_max)`` creats a ``check_schedule`` object with intervals in [t_min, t_max]
* ``schedule.pending()`` return ``true`` if the timer has gone off else ``false``
* ``schedule.check_interval()`` return the number of seconds in the current interval
* ``schedule.update(fraction)`` compute the new interval length and reset the timer

constructor
-----------



.. code-block:: c++

    ising_mpi_sim::ising_mpi_sim(parameters_type const & parameters, boost::mpi::communicator const & comm, double t_min = 1, double t_max = 600)
		: ising_sim(parameters)
        , communicator(comm)
        , schedule(t_min, t_max)
        , clone(comm.rank())
        , binnumber(parameters["BINNUMBER"] | std::min(128, 2 * comm.size()))
    {}

.. code-block:: c++

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

.. code-block:: c++

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


the main function
-----------------

.. code-block:: c++

	int main(int argc, char *argv[]) {

	    try {
	        boost::mpi::environment env(argc, argv);
	        boost::mpi::communicator comm;

	        args options(argc, argv, comm.rank());
	        alps::parameters_type<ising_sim>::type parameters;
	        if (comm.rank() == 0) {
	            std::string suffix = options.inputfile.substr(options.inputfile.find_last_of('.'));
	            if (suffix == ".xml")
	                parameters = alps::make_parameters_from_xml(options.inputfile);
	            else if (suffix == ".h5")
	                alps::hdf5::archive(options.inputfile)["/parameters"] >> parameters;
	            else
	                throw std::runtime_error("Unsupported input formant: " + suffix + "!");
	        }
	        broadcast(comm, parameters);

	        ising_sim sim(parameters, comm);

	        if (options.resume && boost::filesystem::exists(options.checkpointfile))
	            sim.load(options.checkpointfile);

	        sim.run(stop_callback(options.timelimit));

			// make checkpoint
            sim.save(options.checkpointfile);

	        using alps::collect_results;
	        alps::results_type<ising_sim>::type results = collect_results(sim);

	        if (comm.rank() == 0) {
	            std::cout << results << std::endl;
	            alps::hdf5::archive ar(options.outputfile, "w");
	            ar["/parameters"] << parameters;
	            ar["/simulation/results"] << results;
	        }

	    } catch (std::exception const & e) {
	        std::cerr << "Caught exception: " << e.what() << std::endl;
	        return EXIT_FAILURE;
	    } catch (...) {
	        std::cerr << "Caught unknown exception" << std::endl;
	        return EXIT_FAILURE;
	    }
	    return EXIT_SUCCESS;
	}
