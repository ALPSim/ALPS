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

class ising(ngs.base):
    def __init__(self, par, seed = 0):
        ngs.base.__init__(self, par)
        self.completed = 0
        self.length = int(par['L']);
        self.sweeps = 0
        self.thermalization_sweeps = long(par["THERMALIZATION"])
        self.total_sweeps = long(par["SWEEPS"])
        self.beta = 1. / float(par["T"])
        self.tmag = 0.
        self.ten = 0.
        self.sign = 0.
        self.spins = np.ones(self.length)
        for i in range(self.length):
            if self.random() < 0.5:
                self.spins[i] = -1
        # TODO: use shift observables: self.measurements << RealObservable("Energy")
        self.measurements.createRealObservable("Energy")
        self.measurements.createRealObservable("Magnetization")
        self.measurements.createRealObservable("Magnetization^2")
        self.measurements.createRealObservable("Magnetization^4")
        self.measurements.createRealVectorObservable("Correlations")
    def do_update(self):
        for j in range(self.length):
            i = int(float(self.length) * self.random());
            right = 0
            if i + 1 < self.length:
                right = i + 1 
            left = i - 1
            if i - 1 < 0:
                left = self.length - 1
            p = np.exp(2. * self.beta * self.spins[i] * (self.spins[right] + self.spins[left]))
            if p >= 1. or self.random() < p:
                self.spins[i] =- self.spins[i]
    def do_measurements(self):
        self.sweeps += 1
        if self.sweeps > self.thermalization_sweeps:
            tmag = 0;
            ten = 0;
            sign = 1;
            corr = np.zeros(self.length)
            for i in range(self.length):
                tmag += self.spins[i]
                sign *= self.spins[i]
                if i + 1 < self.length:
                    ten += -self.spins[i] * self.spins[i + 1]
                else:
                    ten += -self.spins[i] * self.spins[0]
                for d in range(self.length):
                    corr[d] += self.spins[i] * self.spins[(i + d) % self.length]
            corr / float(self.length)
            ten /= self.length
            tmag /= self.length
            # TODO: use shift observables: self.measurements["energy"] << ten
            self.measurements["Energy"].append(ten)
            self.measurements["Magnetization"].append(tmag)
            self.measurements["Magnetization^2"].append(tmag * tmag)
            self.measurements["Magnetization^4"].append(tmag * tmag * tmag * tmag)
            self.measurements["Correlations"].append(corr)
    def fraction_completed(self):
        if self.sweeps < self.thermalization_sweeps:
            return 0
        else:
            return (self.sweeps - self.thermalization_sweeps) / float(self.total_sweeps)

#TODO: implement, mpi: nur eine base, die alles kann und je nach dem soll er mpi nehmen
#ising.__bases__ = (ngs.mpisim, ) + ising.__bases__

def main(limit, mode, resume, input, output):
    #implement nice argv parsing ...

    sim = ising({
        'L': 100,
        'SWEEPS': 1000,
        'T': 2,
        'THERMALIZATION': 100
    })

    if resume == 't':
        ar = ngs.h5ar(sim.params.valueOrDefault('DUMP', 'dump'), ngs.h5ar.READ)
        sim.load(ar)
        del ar

    if limit == 0:
        sim.run()
    else:
        start = time.time()
        sim.run(lambda: time.time() > start + int(limit))

    ar = ngs.h5ar(sim.params.valueOrDefault('DUMP', 'dump'), ngs.h5ar.REPLACE)
    sim.save(ar)
    del ar

    results = ngs.collectResults(sim)

# TODO: implement
    print "Correlations:           ", results["Correlations"]
    print "Energy:                 ", results["Energy"]

#    print "Mean of Energy:         ", results["Energy"].mean<double>()
#    print "Error of Energy:        ", results["Energy"].mean<double>()
#    print "Mean of Correlations:   ", results["Correlations"].mean<std::vector<double> >()

    print "-2 * Energy / 13:        ", -2. * results["Energy"] / 13.
    print "1 / Correlations        ", 1. / results["Correlations"]
#    print "Energy - Magnetization: ", results["Energy"] - results["Magnetization"]

#    print "Sin(Energy):            ", sin(results["Energy"])
#    print "Tanh(Correlations):     ", tanh(results["Correlations"])

    ngs.saveResults(results, sim.params, ngs.h5ar(output, ngs.h5ar.REPLACE), "/simulation/results")

if __name__ == "__main__":
    apply(main, sys.argv[1:])
