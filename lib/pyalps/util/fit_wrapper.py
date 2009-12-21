#!/opt/local/bin/python2.6

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

class Parameter:	
	def __init__(self, value):
		self.value = value

	def set(self, value):
		self.value = value

	def get(self):
		return self.value

	def __call__(self):
		return self.value

def fit(self,function, parameters, y, x = None):
	def f(params):
		i = 0
		for p in parameters:
			p.set(params[i])
			i += 1
		return y - function(self,x)

	if x == None: x = np.arange(y.shape[0])
	p = [param() for param in parameters]
	optimize.leastsq(f, p)
