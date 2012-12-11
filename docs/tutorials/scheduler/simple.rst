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


.. literalinclude:: code/single/ising.hpp
   :language: c++
   :lines: 33-74

constructor
-----------

.. literalinclude:: code/single/ising.cpp
   :language: c++
   :lines: 5-24

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

.. literalinclude:: code/single/ising.cpp
   :language: c++
   :lines: 26-65

checkpointing
-------------

To save and load checkpoints we use the HDF5 data formant.

.. literalinclude:: code/single/ising.cpp
   :language: c++
   :lines: 67-75

Now we need to tell the hdf5 archive where to store our data. Therefor we implement the following hooks:

.. literalinclude:: code/single/ising.cpp
   :language: c++
   :lines: 77-127

other functions requested by the concept
----------------------------------------

We need a run function which runs runs until we finished or have ran out of time.

.. literalinclude:: code/single/ising.cpp
   :language: c++
   :lines: 129-136

The ``result_names`` function needs to tell us which results are checkpointed

.. literalinclude:: code/single/ising.cpp
   :language: c++
   :lines: 138-143

Since wie save all measurements to the checkpoint, we have no unsaved results:

.. literalinclude:: code/single/ising.cpp
   :language: c++
   :lines: 145-147

If the simulation has finished we want to be able to further process the results:

.. literalinclude:: code/single/ising.cpp
   :language: c++
   :lines: 149-158

the main function
-----------------

.. literalinclude:: code/single/main.cpp
   :language: c++
   :lines: 5-

write the build script
----------------------

In the ``Makefile`` we can use options from the alps library

.. literalinclude:: code/single/Makefile
   :language: makefile
   :lines: 3-
