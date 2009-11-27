import urllib, copy, h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from dataset_core import *
import alea.floatwitherror as fwe

class Hdf5Loader:
	def GetFileNames(self, flist):
		self.files = [f.replace('xml','run1.h5') for f in flist] #will be updated once we have aggregated hdf5-files
		return self.files
		
	# Pre: file is a h5py file descriptor
	# Post: returns a dictionary of all parameters saved in the file
	def ReadParameters(self,file):
		LOP = []
		self.h5param = h5py.File(file)
		pgrp = self.h5param.require_group("/parameters")
		pgrp.visit(LOP.append)
		dict = {'hdf5_file' : file}
		for m in LOP:
			try:
				dict[m] = pgrp[m].value
			except AttributeError:
				pass
		return dict 
		
	def GetResultsPath(self, file):
		path = "/simulation/realizations/0/clones/"
		self.h5f = h5py.File(file)
		sgrp = self.h5f.require_group(path)
		newpath =  path + sgrp.keys()[0] + "/results"
		return newpath
		
	def GetObservableList(self,file):
		p = self.GetResultsPath(file)
		obsgrp = self.h5f.require_group(p)
		return obsgrp.keys()
		
	# Pre: file is a h5py file descriptor
	# Post: returns DataSet with all parameters set
	def ReadMeasurementFromFile(self,flist):
		fs = self.GetFileNames(flist)
		sets = []
		for f in fs:
			path = self.GetResultsPath(f)
			grp = self.h5f.require_group(path)
			params = self.ReadParameters(f)
			obs_list = self.GetObservableList(f)
			kwd = "mean"
			for m in obs_list:
				if kwd in grp[m].keys():
					p_mean = m + "/mean"
					p_error = m + "/error"
					try:
						d = DataSet()
						subset = []
						all_m = grp[p_mean].value
						all_e = grp[p_error].value
						try:
							size = len(all_m)
						except:
							size=0
							subset.append(fwe.FloatWithError(all_m,all_e))
						for i in range(0,size):
							subset.append(fwe.FloatWithError(all_m[i],all_e[i]))
						d.y = np.array(subset)
						d.x =	 np.arange(0,len(d.y))
						d.props['hdf5_path'] = m
						d.props.update(params)
						sets.append(d)
					except AttributeError:
						pass
		return sets