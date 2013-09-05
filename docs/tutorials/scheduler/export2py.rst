Export C++ simulation to python
===============================

We export the interface using ``boost::python``

.. code-block:: c++

	BOOST_PYTHON_MODULE(exported_ising_c) {
	    boost::python::class_<ising_sim>(
	          "sim",
	          boost::python::init<ising_sim::parameters_type const &>()
	    )
	        .def("fraction_completed", &ising_sim::fraction_completed)
	        .def("run", &py_run)
	        .def("resultNames", &ising_sim::result_names)
	        .def("unsavedResultNames", &ising_sim::unsaved_result_names)
	        .def("collectResults", &py_collect_results, py_collect_results_overloads(boost::python::args("names")))
	        .def("save", &py_save)
	        .def("load", &py_load)
	    ;
	}

Not all member funtions can be exported directly, so we need some wrapper.

The ``run`` function take a callback, so we need to create a wrapper around the actual run member function 
to be able to call a python callback function from c++

.. code-block:: c++

	bool py_run_helper(boost::python::object stop_callback) {
	    return boost::python::call<bool>(stop_callback.ptr());
	}
	bool py_run(ising_sim & self, boost::python::object stop_callback) {
	    return self.run(boost::bind(py_run_helper, stop_callback));
	}

There exists two ``collect_results`` functions in the ``ising_sim`` class. ``boost::python`` can handle this:

.. code-block:: c++

	ising_sim::results_type py_collect_results(ising_sim & self, ising_sim::result_names_type const & names = ising_sim::result_names_type()) {
	    return names.size() ? self.collect_results(names) : self.collect_results();
	}
	BOOST_PYTHON_FUNCTION_OVERLOADS(py_collect_results_overloads, py_collect_results, 1, 2)

We only want to export the save/load function from the concept:

.. code-block:: c++

	void py_save(ising_sim const & self, std::string & filename) {
	    self.save(boost::filesystem::path(filename));
	}
	void py_load(ising_sim & self, std::string & filename) {
	    self.load(boost::filesystem::path(filename));
	}


write a build script
--------------------

TBD:


The Python script
-----------------

To use the exported class we need to import the ``pyalps.ngs`` framework ``numpy`` and the exported ising simulation ``exported_ising_c``

.. code-block:: python

	import pyalps.hdf5 as hdf5
	import pyalps.ngs as ngs
	import numpy as np
	import sys, time, getopt

	import exported_ising_c as ising


And a main function could look like

.. code-block:: python

	if __name__ == '__main__':

	    try:
	        optlist, positional = getopt.getopt(sys.argv[1:], 'T:c')
	        args = dict(optlist)
	        try:
	            limit = float(args['-T'])
	        except KeyError:
	            limit = 0
	        resume = True if 'c' in args else False
	        outfile = positional[0]
	    except (IndexError, getopt.GetoptError):
	        print 'usage: [-T timelimit] [-c] outputfile'
	        exit()

	    sim = ising.sim(ngs.params({
	        'L': 100,
	        'THERMALIZATION': 1000,
	        'SWEEPS': 10000,
	        'T': 2
	    }))

	    if resume:
			sim.load(outfile[0:outfile.rfind('.h5')] + '.clone0.h5')

	    if limit == 0:
	        sim.run()
	    else:
	        start = time.time()
	        sim.run(lambda: time.time() > start + float(limit))

		sim.save(outfile[0:outfile.rfind('.h5')] + '.clone0.h5')

	    results = sim.collectResults()
	    print results

	    with hdf5.archive(outfile, 'w') as ar:
	        ar['/parameters'] = sim.parameters
	        ar['/simulation/results'] = results
