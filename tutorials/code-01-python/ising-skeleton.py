import math
import numpy as np
from pyalps.pyalea import *

class Simulation:
	# Init random number generator: 
	# self._rng.random() will give a random float from the interval [0,1)
	_eng = engine()
	_dist = uniform()
	_rng = random(_eng,_dist)
	
	def __init__(self,beta,L):
		self._L = L
		
		# Init exponential map
		self._exp = dict()
		for E in range(-4,5,2): self._exp[E] = math.exp(2*beta*E)
		
		# Init random spin configuration
		self._sites = np.empty([L,L],dtype=np.int)
		for i in range(self._L):
			for j in range(self._L):
				self._sites[i,j] = 2*self.randint(2)-1
		
		# Init observables
		self._energy = RealObservable('E')
		self._magnetization = RealObservable('m')
		self._order = RealObservable('|m|')

    def save(self, filename):
        pyalps.save_parameters(filename, {'L':self.L, 'BETA':self.beta, 'SWEEPS':self.n, 'THERMALIZATION':self.ntherm})
        self.abs_magnetization.save(filename)
        self.energy.save(filename)
        self.magnetization.save(filename)
   
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
		self.print_observable(self._energy)
		self.print_observable(self._order)
		self.print_observable(self._magnetization)
		
	def step(self):
		for s in range(self._L*self._L):
			# Pick random site k=(i,j)
			...
			
			# Measure local energy e = -s_k * sum_{l nn k} s_l
			...
			
			# Flip s_k with probability exp(2 beta e)
			...
		
	def measure(self):
		E = 0.	# energy
		M = 0.	# magnetization
		for i in range(self._L):
			for j in range(self._L):
				E -= ...
				M += ...

		# Add sample to observables
		self._energy << E/(self._L*self._L)
		self._magnetization << M/(self._L*self._L)
		self._order << abs(M)/(self._L*self._L)
		
	# Random int from the interval [0,max)
	def randint(self,max):
		return int(max*self._rng.random())

	# Enforce periodic boundary conditions for lattice indices
	def wrap(self,i):
		return (i+self._L)%self._L
		
	def print_observable(self,o):
		print o, ':\t', o.mean, '+-', o.error, ',\t tau =', o.tau



L = 4	# Linear lattice size
N = 5000	# of simulation steps

print '# L:', L, 'N:', N

# Scan beta range [0,1] in steps of 0.1
for beta in [0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]:
    print '-----------'
    print 'beta =', beta
    sim = Simulation(beta,L)
    sim.run(N/2,N)
    sim.save('ising_L_'+str(L)+'beta_'+str(beta)+'.h5')
