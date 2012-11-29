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
import sys, time, traceback, getopt

class isingSim:

    # TODO: how do we deal with typedefs?

    def __init__(self, params):
        self.random = ngs.random01(params.valueOrDefault('SEED', 42))
        self.parameters = params
        self.measurements = {
            #jan: constructor syntax would be nicer: ngs.RealObservable('X')
            'Energy':  ngs.createRealObservable('Energy'),
            'Magnetization': ngs.createRealObservable('Magnetization'),
            'Magnetization^2': ngs.createRealObservable('Magnetization^2'),
            'Magnetization^4': ngs.createRealObservable('Magnetization^4'),
            'Correlations': ngs.createRealVectorObservable('Correlations')
        }

        self.length = int(params['L'])
        self.sweeps = 0
        self.thermalization_sweeps = long(params['THERMALIZATION'])
        self.total_sweeps = long(params['SWEEPS'])
        self.beta = 1. / float(params['T'])
        self.spins = np.array([(-x if self.random() < 0.5 else x) for x in np.ones(self.length)])
        
        self.realization = '0'
        self.clone = '0'

    def update(self):
        for j in range(self.length):
            i = int(float(self.length) * self.random())
            right = i + 1 if i + 1 < self.length else 0
            left = self.length - 1 if i - 1 < 0 else i - 1
            p = np.exp(2. * self.beta * self.spins[i] * (self.spins[right] + self.spins[left]))
            if p >= 1. or self.random() < p:
                self.spins[i] =- self.spins[i]

    def measure(self):
        self.sweeps += 1
        if self.sweeps > self.thermalization_sweeps:
            tmag = 0
            ten = 0
            sign = 1
            corr = np.zeros(self.length)
            for i in range(self.length):
                tmag += self.spins[i]
                sign *= self.spins[i]
                ten += -self.spins[i] * self.spins[i + 1 if i + 1 < self.length else 0]
            for d in range(self.length):
                corr[d] = np.inner(self.spins, np.roll(self.spins, d)) / float(self.length)
            ten /= self.length
            tmag /= self.length
            self.measurements['Energy'] << ten
            self.measurements['Magnetization'] << tmag
            self.measurements['Magnetization^2'] << tmag**2
            self.measurements['Magnetization^4'] << tmag**4
            self.measurements['Correlations'] << corr

    def fraction_completed(self):
        return 0 if self.sweeps < self.thermalization_sweeps else (self.sweeps - self.thermalization_sweeps) / float(self.total_sweeps)

    def save(self, filename):
        with ngs.archive(filename, 'w') as ar:
            ar['/'] = self

    def load(self, filename):
        with ngs.archive(filename, 'r') as ar:
            self = ar['/']

    def run(self, stopCallback):
        stopped = False
        while True:
            self.update()
            self.measure()
            stopped = stopCallback()
            if (stopped or self.fraction_completed() >= 1.):
                return not stopped

    def result_names(self):
        return self.measurements.keys()

    def unsaved_result_names(self):
        return self.result_names_type(self)

    def collectResults(self, names = None):
        if names == None:
            names = self.result_names()
        partial_results = {}
        for name in names:
            partial_results[name] = ngs.observable2result(self.measurements[name])
        return partial_results

    def save(self, ar):
    
        try:

            ar["/parameters"] = self.parameters
            context = ar.context
            ar.set_context("/simulation/realizations/" + self.realization + "/clones/" + self.clone)

            ar["length"] = self.length # TODO: where to put the checkpoint informations?
            ar["sweeps"] = self.sweeps
            ar["thermalization_sweeps"] = self.thermalization_sweeps
            ar["beta"] = self.beta
            ar["spins"] = self.spins
            ar["measurements"] = self.measurements
            ar["engine"] = self.random

            ar.set_context(context)

        except:
            traceback.print_exc(file=sys.stderr)
            raise

    def load(self,  ar):
    
        try:
        
            params.load(ar["/parameters"]) # TODO: do we want to load the parameters?

            context = ar.context
            ar.set_context("/simulation/realizations/" + self.realization + "/clones/" + self.clone)

            self.length = ar["length"]
            self.sweeps = ar["sweeps"]
            self.thermalization_sweeps = ar["thermalization_sweeps"]
            self.beta = ar["beta"]
            self.spins = ar["spins"]
            self.measurements = ar["measurements"]
            self.random.load(ar["engine"])

            ar.set_context(context)

        except:
            traceback.print_exc(file=sys.stderr)
            raise

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

    sim = isingSim(ngs.params({
        'L': 100,
        'THERMALIZATION': 100,
        'SWEEPS': 1000,
        'T': 2
    }))

    if resume:
        try:
            with ngs.archive(outfile[0:outfile.rfind('.h5')] + '.clone0.h5', 'r') as ar:
                sim.load(ar)
        except ArchiveNotFound: pass

    if limit == 0:
        sim.run()
    else:
        start = time.time()
        sim.run(lambda: time.time() > start + float(limit))

    if resume:
        with ngs.archive(outfile[0:outfile.rfind('.h5')] + '.clone0.h5', 'w') as ar:
            ar['/'] = sim

    results = sim.collectResults() # TODO: how should we do that?
    for key, value in results.iteritems():
        print "{}: {}".format(key, value)

    with ngs.archive(outfile, 'w') as ar: # TODO: how sould we name archive? ngs.hdf5.archive?
        ar['/parameters'] = sim.parameters
        ar['/simulation/results'] = results
