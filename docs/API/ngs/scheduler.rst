Scheduler
=========

Supported mc simulation types
-----------------------------

  * Serial simulation - running a mc simulation in one process
  * MPI simulation - running n identical mc simulations on n mpi processes
  * Multithreaded simulation - running n identical mc simulations on n threads
  * Multithreaded MPI simulation - running n * m identical mc simulations on n mpi processes with n threads each

  * Serial clone - running 1 clone on 1 process
  * MPI clone - running 1 clone on m mpi processes
  * Multithreaded clone - running 1 clone on n threads of one process
  * Multithreaded MPI clone - running 1 clone on n mpi processes with m threads each

  * Processes can optionally have separate master/interface and worker threads

Stop Callback concept
---------------------

    * Let
        * ``CB`` be a model of ''Stop Callback'' and ``stop_callback`` be an object of type ``CB``. 
    * Associated types:
        * ``alps::ngs::external_signal`` is derived from ``boost::thread_interrupted``
    * Valid expressions:
        * ``static_cast<boost::function<void()> >(stop_callback)`` stop_callback is convertible to void()
        * ``stop_callback()`` checks if the simulations needs to stop:
            * ``throws boost::thread_interrupted``: the time limit has expired or the simulation has finished
            * ``throws alps::ngs::external_signal``: the process has received a system signal

Simulation concepts
-------------------

    * Let
        * ``S`` be a model of ''Simulation'' and ``s`` be an object of type ``S``. 
        * ``P`` be a model of ''Parameter'' and ``p`` be an object of type ``P``.
        * ``H`` be a model of ''Hirarchical Archive'' and ``ar`` be an object of type ``H``
        * ``CB`` be a model of ''Stop Callback'' and ``stop_callback`` be an object of type ``CB``
        * ``path`` be a ``boost::filesystem::path`` pointing to a file, that does not necessarily have to exist, but the branch path exists;
        * ``h5path`` be a std::string
        * ``names`` is a sequence of std::strings
    * Associated types:
        * ``parameter_type``  is a type of concept "Parameter" storing the parameter informations
        * ``results_type``  is a type of concept "Immutable Accumulator"
    * Valid expressions:
        * ``alps::results_type<S>::type`` is the type of the associative container returned by ``collect_results(c)``
        * ``alps::result_names_type<S>::type`` is the type of the sequence returned by ``result_names(s)`` and ``unsaved_result_names(s)``
        * ``S(p)``, ``S s(p)``
        * ``ar[h5path] << s``, ``s.save(path)`` writes a checkpoint, not necessarily containing all accumulators -- could also be empty if we always use a completely new simulation upon restart
        * ``ar[h5path] >> s``, ``s.load(path)`` reads a checkpoint, not necessarily containing all accumulators
        * ``s.run(stop_callback)``, runs the simulation, regularly calling ``stop_callback``. Returns a value convertible to ``bool``:
            * ``true``: the simulation has finished 
            * ``false``: the simulation is still running '''TBD: in which case is this possible?'''
        * ``unsaved_result_names(s)`` returns a sequence of strings containing the names of observables not saved by calling ``s.save(path)``
        * ``result_names(s)`` returns a sequence of strings containing the names of all observables measured.
        * ``collect_results(s)`` returns an associative container of immutable accumulators
        * ``collect_results(s, names)``   returns an associative container of immutable accumulator whose keys are in the sequence ``names``
        * ``fraction_completed(s)`` returns a value convertible to floating point indicating the fraction of work done.

MPI Simulation concepts
-----------------------

    * The ''MPI Simulation'' concept is a refinement of the ''Simulation'' concept
    * Let
        * ``comm`` be a ``boost::mpi::communicator``
    * Valid expressions:
        * ``S s(comm, parms)``, ``S(comm, parms)``
        * ``desired_num_processes<S>(p)`` returns a pair of values convertible to an unsigned integral type indicating the minimum and maximum number of processes desired.
