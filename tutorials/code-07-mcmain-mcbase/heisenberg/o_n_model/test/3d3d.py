import sys, time, getopt

sys.path.extend(['../build/', '/opt/alps/lib/'])
import pyalps
import ndsim_c as ndsim

if __name__ == '__main__':

    outfile = '3d3d_result.h5'

    sim = ndsim.heisenberg_sim(pyalps.ngs.params({
        'LATTICE_LIBRARY'   :   '/opt/alps/lib/xml/lattices.xml',
        'LATTICE'           :   'simple cubic lattice',
        'L'                 :   '20',
        'THERMALIZATION'    :   '1000',
        'SWEEPS'            :   '100000',
        'T'                 :   '2'
    }))

    sim.run(lambda: False)

    results = sim.collectResults()
    print results

    with pyalps.hdf5.archive(outfile, 'w') as ar:
        ar['parameters'] = sim.parameters
        ar['simulation/results'] = results
