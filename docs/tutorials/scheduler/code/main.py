
	import pyalps.hdf5 as hdf5
	import pyalps.ngs as ngs
	import numpy as np
	import sys, time, getopt

	import ising_sim as ising

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