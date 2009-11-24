import urllib, copy, h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from dataset_core import *

class Hdf5Loader:
	# Pre: file is a h5py file descriptor
	# Post: returns a dictionary of all parameters saved in the file
	def ReadParameters(self,file):
		LOP = []
		self.h5param = h5py.File(file)
		pgrp = self.h5param.require_group("/parameters")
		pgrp.visit(LOP.append)
		dict = {}
		for m in LOP:
			try:
				dict[m] = pgrp[m].value
			except AttributeError:
				pass
		return dict 
		
	# Pre: file is a h5py file descriptor
	# Post: returns a list of all measurements saved in the file
	def FindAllMeasurements(self,file):
		list_of_paths = []
		self.h5f = h5py.File(file)
		sgrp = self.h5f.require_group("/simulation")
		sgrp.visit(list_of_paths.append)
		return list_of_paths
	
	# Pre: file is a h5py file descriptor
	# Post: returns DataSet with all parameters set
	def ReadMeasurementFromFile(self,file):
		LOM = self.FindAllMeasurements(file)
		grp = self.h5f.require_group("/simulation")
		params = self.ReadParameters(file)
		
		sets = []
		for m in LOM:
			splitpath = m.rpartition("/")
			last_entry = splitpath[2]
			obs = splitpath[0].rpartition("/")[2]
			if last_entry == "mean" :
				try:
					d = DataSet()
					try:
						d.y = grp[m].value
						d.x = np.arange(0,len(d.y))
					except TypeError : 
						d.y = [grp[m].value]
						d.x = [0]
					d.x = np.array(d.x)
					d.y = np.array(d.y)	
					d.props['observable'] = obs
					d.props['hdf5_file'] = file
					d.props.update(params)
					sets.append(d)
				except AttributeError:
					pass
		return sets		