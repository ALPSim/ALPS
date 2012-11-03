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
import numpy as np
import sys, time

#TODO: enable mpi

import ngsising_c as ising

def main(limit, resume, output):
    #implement nice argv parsing ...

    sim = ising.sim(ngs.params({
        'L': 100,
        'THERMALIZATION': 1000,
        'SWEEPS': 10000,
        'T': 2
    }))

    if resume == 't':
        ar = ngs.h5ar(sim.params.valueOrDefault('DUMP', 'dump'), 'r')
        sim.load(ar)
        del ar

    if limit == 0:
        sim.run()
    else:
        start = time.time()
        sim.run(lambda: time.time() > start + int(limit))

    p = sim.params

    ar = ngs.h5ar(sim.params.valueOrDefault('DUMP', 'dump'), 'w')
    sim.save(ar)
    del ar

    results = ngs.collectResults(sim)
    print results
    ngs.saveResults(results, sim.params, ngs.h5ar(output, 'w'), "/simulation/results")

if __name__ == "__main__":
    apply(main, sys.argv[1:])