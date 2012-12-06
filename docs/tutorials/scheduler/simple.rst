Simple ising sumulation
=======================

This tutoral describes how to write a simulation without using the alps scheduler base classes. As an example use a classical ising model.

concept
-------

A Simulation needs to follow the **simulaton conpept**:

* Let
    * ``S`` be a model of ''Simulation'' and ``s`` be an object of type ``S``. 
    * ``P`` be a model of ''Parameter'' and ``p`` be an object of type ``P``.
    * ```stop_callback`` is convertible to ``boost::function<bool()>```
    * ``in_path`` ``out_path`` be of type ``boost::filesystem::path`` pointing to a file, that does not necessarily have to exist, but the branch path exists
    * ``names`` be an object of type ``result_names_type``
* Associated types:
    * ``parameter_type`` is a type of concept ''Parameter'' storing the parameter informations
    * ``results_type`` is a type of concept ''Immutable Accumulator''
    * ``result_names_type`` is a sequence of std::strings
* Valid expressions:
    * ``alps::results_type<S>::type`` is the type of the associative container returned by ``collect_results(c)``
    * ``alps::result_names_type<S>::type`` is the type of the sequence returned by ``result_names(s)`` and ``unsaved_result_names(s)``
    * ``S(p)``, ``S s(p)``
    * ``s.save(out_path)`` writes a checkpoint, not necessarily containing all accumulators -- could also be empty if we always use a completely new simulation upon restart
    * ``s.load(in_path)`` reads a checkpoint, not necessarily containing all accumulators
    * ``s.run(stop_callback)``, runs the simulation, regularly calling ``stop_callback``. Returns the return value of ``stop_callback`` as a convertible to ``bool``
    * ``unsaved_result_names(s)`` returns a sequence of strings containing the names of observables not saved by calling ``s.save(path)``
    * ``result_names(s)`` returns a sequence of strings containing the names of all observables measured.
    * ``collect_results(s)`` returns an associative container of immutable accumulators
    * ``collect_results(s, names)``   returns an associative container of immutable accumulator whose keys are in the sequence ``names``
    * ``fraction_completed(s)`` returns a value convertible to floating point indicating the fraction of work done.

synopsis
--------

Given the concept above, our simulation has the following synopsis:

.. code-block:: c++
    :linenos:

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

constructor
-----------

.. code-block:: c++
    :linenos:

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

First we initialize the parameter class. The ``alps::params`` class proviedes a simple interface to access the parameters:

* ``value = params[key]`` if ``key`` exists, ``params[key]`` is assigned else an exception is thrown
* ``params[key].cast<T>()`` returns a value convertable to T
* ``params[key] | value`` the return type is the same as ``value``. If ``key`` exists, ``params[key]`` is returned else value.
* ``params.defined(key)`` result convertible to ``bool``, indicating the existance of ``key``

Now we initialize the random number generator and the state variables.

Inside the constructor we initialize the measurements:

* ``measurements << alps::ngs::RealObservable("Energy")`` initializes an measurement of double with mean, error and binning analysis.
* ``measurements << alps::ngs::RealVectorObservable("Correlations")`` initializes an measurement of vector<double> with mean, error and binning analysis.

update / measure / fraction_complete
------------------------------------

These Functions contains the actual simulation. ``update`` does one montecarlo step, ``measure`` updates the measurements and ``fraction_complete`` 
return the progress of the simulation.

.. code-block:: c++
    :linenos:

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

checkpointing
-------------

To save and load checkpoints we use the HDF5 data formant.

.. code-block:: c++
    :linenos:

    void ising_sim::save(boost::filesystem::path const & filename) const {
        alps::hdf5::archive ar(filename, "w");
        ar << *this;
    }
    
    void ising_sim::load(boost::filesystem::path const & filename) {
        alps::hdf5::archive ar(filename);
        ar >> *this;
    }

Now we need to tell the hdf5 archive where to store our data. Therefor we implement the following hooks:

.. code-block:: c++
    :linenos:

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

other functions requested by the concept
----------------------------------------

We need a run function which runs runs until we finished or have ran out of time.

.. code-block:: c++
    :linenos:

    bool ising_sim::run(boost::function<bool ()> const & stop_callback) {
        bool stopped = false;
        do {
            update();
            measure();
        } while(!(stopped = stop_callback()) && fraction_completed() < 1.);
        return !stopped;
    }

The ``result_names`` function needs to tell us which results are checkpointed

.. code-block:: c++
    :linenos:

    ising_sim::result_names_type ising_sim::result_names() const {
        result_names_type names;
        for(observables_type::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
            names.push_back(it->first);
        return names;
    }

Since wie save all measurements to the checkpoint, we have no unsaved results:

.. code-block:: c++
    :linenos:

    ising_sim::result_names_type ising_sim::unsaved_result_names() const {
        return result_names_type(); 
    }

If the simulation has finished we want to be able to further process the results:

.. code-block:: c++
    :linenos:

    ising_sim::results_type ising_sim::collect_results() const {
        return collect_results(result_names());
    }
    ising_sim::results_type ising_sim::collect_results(result_names_type const & names) const {
        results_type partial_results;
        for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it)
            partial_results.insert(*it, alps::mcresult(measurements[*it]));
        return partial_results;
    }

the main function
-----------------

.. code-block:: c++
    :linenos:

	int main(int argc, char *argv[]) {

	    try {
	        args options(argc, argv);

	        alps::parameters_type<ising_sim>::type parameters;

	        std::string suffix = options.inputfile.substr(options.inputfile.find_last_of('.'));
	        if (suffix == ".xml")
	            parameters = alps::make_parameters_from_xml(options.inputfile);
	        else if (suffix == ".h5")
	            alps::hdf5::archive(options.inputfile)["/parameters"] >> parameters;
	        else
	            throw std::runtime_error("Unsupported input format: " + suffix + "!");

	        ising_sim sim(parameters);

	        if (options.resume && boost::filesystem::exists(options.checkpointfile))
	            sim.load(options.checkpointfile);

	        sim.run(stop_callback(options.timelimit));
        
			// make checkpoint
	        sim.save(options.checkpointfile);
        
	        using alps::collect_results;
	        alps::results_type<ising_sim>::type results = collect_results(sim);

	        std::cout << results << std::endl;
	        alps::hdf5::archive ar(options.outputfile, "w");
	        ar["/parameters"] << parameters;
	        ar["/simulation/results"] << results;

	    } catch (std::exception const & e) {
	        std::cerr << "Caught exception: " << e.what() << std::endl;
	        return EXIT_FAILURE;
	    } catch (...) {
	        std::cerr << "Caught unknown exception" << std::endl;
	        return EXIT_FAILURE;
	    }
	    return EXIT_SUCCESS;
	}

write a build script
--------------------

TBD:
