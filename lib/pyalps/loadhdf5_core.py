import urllib, copy, h5py, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from dataset import ResultFile
from dataset import DataSet
from floatwitherror import FloatWithError as fwe

import pyalps.pytools as pt # the C++ conversion functions

# or the C++ class as alternative
#from pyalps.pyalea import value_with_error as fwe

class Hdf5Missing(Exception):
    def __init__(self,what):
        self.what = what
    
    def __str__(self):
        return 'Failed to find ' + self.what

class Hdf5Loader:
    def GetFileNames(self, flist):
        self.files = [f.replace('xml','h5') for f in flist] 
        return self.files
        
    # Pre: file is a h5py file descriptor
    # Post: returns a dictionary of all parameters saved in the file
    def ReadParameters(self,proppath):
        LOP = []
        dict = {'filename' : file}
        try:
            pgrp = self.h5f.require_group(proppath)
            pgrp.visit(LOP.append)
            for m in LOP:
                try:
                    dict[m] = pgrp[m].value
                except AttributeError:
                    pass
                except KeyError:
                    pass
        except KeyError:
            raise Hdf5Missing(proppath + ' in ReadParameters()')
        return dict 
        
    def GetProperties(self,flist,proppath):
        fs = self.GetFileNames(flist)
        resultfiles = []
        for f in fs:
            rfile = ResultFile()
            rfile.props = self.ReadParameters(f,proppath)
            rfile.props["ObservableList"] = self.GetObservableList(f)
            resultfiles.append(rfile)
        return resultfiles
        
    def GetObservableList(self,respath):
        try:
            obsgrp = self.h5f.require_group(respath)
        except KeyError:
            raise Hdf5Missing(respath + ' in GetObservableList()')
        olist = [pt.hdf5_name_encode(obs) for obs in obsgrp.keys()]
        return olist
        
    # Pre: file is a h5py file descriptor
    # Post: returns DataSet with all parameters set
    def ReadMeasurementFromFile(self,flist,statvar,proppath,respath,measurements=None):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            self.h5f = h5py.File(f)
            self.h5fname = f
            
            list_ = self.GetObservableList(respath)
            # this is exception-safe in the sense that it's also required in the line above
            grp = self.h5f.require_group(respath)
            params = self.ReadParameters(proppath)
            obslist = []
            if measurements == None:
                obslist = list_
            else:
                obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
            subset=[]
            for m in obslist:
                try:
                    d = DataSet()
                    if "mean" in grp[m].keys() and "error" in grp[m+"/mean"].keys():
                        mean = grp[m+"/mean/value"].value
                        error = grp[m+"/mean/error"].value
                        d.props['count'] = grp[m+"/count"].value
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
                    d.props['hdf5_path'] = respath + m
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
                    print "Could not create DataSet"
                    pass
        return sets
