# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Copyright (C) 2009 - 2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
#                              Brigitte Surer <surerb@phys.ethz.ch>
#
# Distributed under the Boost Software License, Version 1.0. (See accompany-
# ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
#
##############################################################################

import math
import pyalps
import pyalps.alea as alpsalea
import pyalps.pytools as alpstools

class Simulation:
    # Seed random number generator: self.rng() will give a random float from the interval [0,1)
    rng = alpstools.rng(42)
    
    def __init__(self,beta,L):
        self.L = L
        self.beta = beta
        
        # Init exponential map
        self.exp_table = dict()
        for E in range(-4,5,2): 
          self.exp_table[E] = math.exp(2*beta*E)
        
        # Init random spin configuration
        self.spins = [ [2*self.randint(2)-1 for j in range(L)] for i in range(L) ]
        
        # Init observables
        self.energy = alpsalea.RealObservable('E')
        self.magnetization = alpsalea.RealObservable('m')
        self.abs_magnetization = alpsalea.RealObservable('|m|')
    
    def run(self,ntherm,n):
        # Thermalize for ntherm steps
        self.n = n
        self.ntherm = ntherm
        while ntherm > 0:
            self.step()
            ntherm = ntherm-1
        
        # Run n steps
        while n > 0:
            self.step()
            self.measure()
            n = n-1
        
        # Print observables
        print '|m|:\t', self.abs_magnetization.mean, '+-', self.abs_magnetization.error, ',\t tau =', self.abs_magnetization.tau
        print 'E:\t', self.energy.mean, '+-', self.energy.error, ',\t tau =', self.energy.tau
        print 'm:\t', self.magnetization.mean, '+-', self.magnetization.error, ',\t tau =', self.magnetization.tau
    
    def step(self):
        for s in range(self.L*self.L):
            # Pick random site k=(i,j)
            i = self.randint(self.L)
            j = self.randint(self.L)
            
            # Measure local energy e = -s_k * sum_{l nn k} s_l
            e = self.spins[(i-1+self.L)%self.L][j] + self.spins[(i+1)%self.L][j] + self.spins[i][(j-1+self.L)%self.L] + self.spins[i][(j+1)%self.L]
            e *= -self.spins[i][j]
            
            # Flip s_k with probability exp(2 beta e)
            if e > 0 or self.rng() < self.exp_table[e]:
                self.spins[i][j] = -self.spins[i][j]
                
    def save(self, filename):
        pyalps.save_parameters(filename, {'L':self.L, 'BETA':self.beta, 'SWEEPS':self.n, 'THERMALIZATION':self.ntherm})
        self.abs_magnetization.save(filename)
        self.energy.save(filename)
        self.magnetization.save(filename)
       
        
    def measure(self):
        E = 0.    # energy
        M = 0.    # magnetization
        for i in range(self.L):
            for j in range(self.L):
                E -= self.spins[i][j] * (self.spins[(i+1)%self.L][j] + self.spins[i][(j+1)%self.L])
                M += self.spins[i][j]
            
        # Add sample to observables
        self.energy << E/(self.L*self.L)
        self.magnetization << M/(self.L*self.L)
        self.abs_magnetization << abs(M)/(self.L*self.L)
        
    # Random int from the interval [0,max)
    def randint(self,max):
        return int(max*self.rng())
        



