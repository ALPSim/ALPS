 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   #
 #                                                                                 #
 # SPDX-License-Identifier: MIT
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pyalps.ngs as ngs
import numpy as np
import sys

class sim(ngs.mcbase):
    def __init__(self, par, seed = 42):
        ngs.mcbase.__init__(*(self, par, seed))

        self.measurements << ngs.RealObservable("Energy")
        self.measurements << ngs.RealObservable("Magnetization")
        self.measurements << ngs.RealObservable("Magnetization^2")
        self.measurements << ngs.RealObservable("Magnetization^4")
        self.measurements << ngs.RealVectorObservable("Correlations")

        self.length = int(self.parameters['L'])
        self.sweeps = 0
        self.thermalization_sweeps = long(self.parameters['THERMALIZATION'])
        self.total_sweeps = long(self.parameters['SWEEPS'])
        self.beta = 1. / float(self.parameters['T'])
        self.spins = np.array([(-x if self.random() < 0.5 else x) for x in np.ones(self.length)])

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

    def save(self, ar):
    
        try:
            ngs.mcbase.save(self, ar)
            ar["checkpoint/sweeps"] = self.sweeps
            ar["checkpoint/spins"] = self.spins

        except:
            traceback.print_exc(file=sys.stderr)
            raise

    def load(self,  ar):

        try:
            ngs.mcbase.load(self, ar)

            self.length = int(self.parameters["L"]);
            self.thermalization_sweeps = int(self.parameters["THERMALIZATION"]);
            self.total_sweeps = int(self.parameters["SWEEPS"]);
            self.beta = 1. / double(self.parameters["T"]);

            self.sweeps = ar["checkpoint/sweeps"]
            self.spins = ar["checkpoint/spins"]

        except:
            traceback.print_exc(file=sys.stderr)
            raise
