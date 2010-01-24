import urllib, copy, h5py
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import optimize

from dataset import ResultFile
from dataset import DataSet
from floatwitherror import FloatWithError as fwe

import pyalps.pytools as pt # the C++ conversion functions

# or the C++ class as alternative
#from pyalps.pyalea import value_with_error as fwe

class Hdf5Loader:
    def GetFileNames(self, flist):
        self.files = [f.replace('xml','h5') for f in flist] 
        return self.files
        
    # Pre: file is a h5py file descriptor
    # Post: returns a dictionary of all parameters saved in the file
    def ReadParameters(self,file):
        LOP = []
        self.h5param = h5py.File(file)
        pgrp = self.h5param.require_group("/parameters")
        pgrp.visit(LOP.append)
        dict = {'filename' : file}
        for m in LOP:
            try:
                dict[m] = pgrp[m].value
            except AttributeError:
                pass
        return dict 
        
    def GetProperties(self, flist):
        fs = self.GetFileNames(flist)
        resultfiles = []
        for f in fs:
            rfile = ResultFile()
            rfile.props = self.ReadParameters(f)
            rfile.props["ObservableList"] = self.GetObservableList(f)
            resultfiles.append(rfile)
        return resultfiles
        
    def GetResultsPath(self, file):
        self.h5f = h5py.File(file)
        path = "/simulation/results/"
        return path
        
    def GetObservableList(self,file):
        p = self.GetResultsPath(file)
        obsgrp = self.h5f.require_group(p)
        return obsgrp.keys()
        
    # Pre: file is a h5py file descriptor
    # Post: returns DataSet with all parameters set
    def ReadMeasurementFromFile(self,flist,statvar,measurements=None):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            path = self.GetResultsPath(f)
            grp = self.h5f.require_group(path)
            params = self.ReadParameters(f)
            list = self.GetObservableList(f)
            obslist = []
            if measurements == None:
                obslist = list
            else:
                obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list]
            try:
                d = DataSet()
                subset=[]
                for m in obslist:
                    if "mean" in grp[m].keys() and "error" in grp[m+"/mean"].keys():
                        mean = grp[m+"/mean/value"].value
                        error = grp[m+"/mean/error"].value
                        print "reading count"
                        d.props['count'] = grp[m+"/count"].value
                        print "value of count", grp[m+"/count"].value 
                        try:
                            size = len(mean)
                            subset = [fwe(mean[i],error[i]) for i in range(0,size)]
                        except:
                            size=0
                            subset = [fwe(mean,error)]
                    elif "mean" in grp[m].keys():
                        value = grp[m+"/mean/value"].value
                        try:
                            size=len(value)
                            subset = [float(value[i]) for i in range(0,size)]
                        except:
                            size=0
                            subset = [float(value)]
                    d.y = np.array(subset)
                    if "labels" in grp[m].keys():
                        d.x = np.array(grp[m+"/labels"].value)
                    else:
                        d.x = np.arange(0,len(d.y))
                    d.props['hdf5_path'] = path + m
                    d.props['observable'] = pt.hdf5_name_decode(m)
                    d.props.update(params)
                    if "timeseries" in statvar:
                        tslist = grp[m+"/timeseries"].keys()
                        for l in tslist:
                            d.props[l] = grp[m+"/timeseries/"+l].value
                    for s in statvar:
                        if s in grp[m].keys():
                            try:
                                d.props[s] = grp[m+"/"+s+"/value"].value
                            except:
                                pass
                    sets.append(d)
            except AttributeError:
                print "could not create DataSet"
                pass
        return sets
