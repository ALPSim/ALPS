 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   #
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
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs
import sys, time, traceback, getopt
import sys, time

import ising

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

    sim = ising.sim({
        'L': 100,
        'THERMALIZATION': 100,
        'SWEEPS': 1000,
        'T': 2
    })

    if resume:
        try:
            with hdf5.archive(outfile[0:outfile.rfind('.h5')] + '.clone0.h5', 'r') as ar:
                sim.load(ar)
        except ArchiveNotFound: pass

    if limit == 0:
        sim.run(lambda: False)
    else:
        start = time.time()
        sim.run(lambda: time.time() > start + float(limit))

    with hdf5.archive(outfile[0:outfile.rfind('.h5')] + '.clone0.h5', 'w') as ar:
        ar['/'] = sim

    results = ngs.collectResults(sim)
    print results
    with hdf5.archive(outfile, 'w') as ar:
        ngs.saveResults(results, sim.paramters, ar, "/simulation/results")
