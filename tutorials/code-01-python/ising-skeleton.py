import math
import pyalps.pyalea as alpsalea
import pyalps.pytools as alpstools

class Simulation:
    # Seed random number generator: self.rng() will give a random float from the interval [0,1)
    rng = alpstools.rng(42)
    
    def __init__(self,beta,L):
        self.L = L
        
        # Init exponential map
        self.exp_table = dict()
        for E in range(-4,5,2): self.exp_table[E] = math.exp(2*beta*E)
        
        # Init random spin configuration
        self.spins = [ [2*self.randint(2)-1 for j in range(L)] for i in range(L) ]
        
        # Init observables
        self.energy = alpsalea.RealObservable('E')
        self.magnetization = alpsalea.RealObservable('m')
        self.abs_magnetization = alpsalea.RealObservable('|m|')

    def run(self,ntherm,n):
        # Thermalize for ntherm steps
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
            ...
            
            # Measure local energy e = -s_k * sum_{l nn k} s_l
            ...
            
            # Flip s_k with probability exp(2 beta e)
            ...
        
    def measure(self):
        E = 0.	# energy
        M = 0.	# magnetization
        for i in range(self.L):
            for j in range(self.L):
                E -= ...
                M += ...

        # Add sample to observables
        self.energy << E/(self.L*self.L)
        self.magnetization << M/(self.L*self.L)
        self.abs_magnetization << abs(M)/(self.L*self.L)
        
    # Random int from the interval [0,max)
    def randint(self,max):
        return int(max*self.rng())

if __name__ == '__main__':
    L = 4       # Linear lattice size
    N = 5000    # of simulation steps

    print '# L:', L, 'N:', N

    # Scan beta range [0,1] in steps of 0.1
    for beta in [0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]:
        print '-----------'
        print 'beta =', beta
        sim = Simulation(beta,L)
        sim.run(N/2,N)
