 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   #
 #                                                                                 #
 # This software is part of the ALPS libraries, published under the ALPS           #
 # Library License; you can use, redistribute it and/or modify it under            #
 # the terms of the license, either version 1 or (at your option) any later        #
 # version.                                                                        #
 #                                                                                 #
 # You should have received a copy of the ALPS Library License along with          #
 # the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       #
 # available from http://alps.comp-phys.org/.                                      #
 #                                                                                 #
 #  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        #
 # FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       #
 # SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       #
 # FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     #
 # ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     #
 # DEALINGS IN THE SOFTWARE.                                                       #
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pyalps.ngs as ngs
import mpi
import numpy as np
import sys, time

import ising

def main(limit, resume, output):
    #implement nice argv parsing ...

    sim = ising.sim({
        'L': 100,
        'SWEEPS': 1000,
        'T': 2,
        'THERMALIZATION': 100
    }, mpi.world if mpi.size > 1 else None)

    if resume == 't':
        ar = ngs.h5ar(sim.params.valueOrDefault('DUMP', 'dump.h5'), 'r')
        sim.load(ar)
        del ar

    if limit == 0:
        sim.run()
    else:
        start = time.time()
        sim.run(lambda: time.time() > start + int(limit))

    ar = ngs.h5ar(sim.params.valueOrDefault('DUMP', 'dump.h5') + ('.' + str(mpi.rank) if mpi.size > 1 else ''), 'a')
    sim.save(ar)
    del ar

    results = ngs.collectResults(sim)

    if mpi.rank == 0:
        print "Correlations:           ", results["Correlations"]
        print "Energy:                 ", results["Energy"]

        print "Mean of Energy:         ", results["Energy"].mean
        print "Error of Energy:        ", results["Energy"].error
        print "Mean of Correlations:   ", results["Correlations"].mean

        print "-2 * Energy / 13:        ", -2. * results["Energy"] / 13.
        print "1 / Correlations        ", 1. / results["Correlations"]
        print "Energy - Magnetization: ", results["Energy"] - results["Magnetization"]

        print "Sin(Energy):            ", results["Energy"].sin()
        print "Tanh(Correlations):     ", results["Correlations"].tanh()

        ngs.saveResults(results, sim.params, ngs.h5ar(output, 'a'), "/simulation/results")

if __name__ == "__main__":
    apply(main, sys.argv[1:])
