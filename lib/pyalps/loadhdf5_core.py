import urllib, copy, h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from dataset import ResultFile
from dataset import DataSet
from floatwitherror import FloatWithError as fwe
# or the C++ class as alternative
#from pyalps.pyalea import value_with_error as fwe

class Hdf5Loader:
    def GetFileNames(self, flist):
        self.files = [f.replace('.out.xml','.clone1.h5') for f in flist] #will be updated once we have aggregated hdf5-files
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
        path = "/simulation/realizations/0/clones/0/results/0/worker/0/"
        self.h5f = h5py.File(file)
        sgrp = self.h5f.require_group(path)
        # catch some exceptions and write more reasonable error
        # right now, the error is quite useless
        newpath =  path + sgrp.keys()[0] + "/results"
        return path #newpath
        
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
            list_ = self.GetObservableList(f)
            obslist = []
            if measurements == None:
                obslist = list_
            else:
                obslist = [obs for obs in measurements if obs in list_]
            kwd = "mean"
            for m in obslist:
                if kwd in grp[m].keys():
                    p_mean = m + "/mean/value"
                    p_error = m + "/mean/error"
                    try:
                        d = DataSet()
                        subset = []
                        all_m = grp[p_mean].value
                        all_e = grp[p_error].value
                        try:
                            size = len(all_m)
                        except:
                            size=0
                            subset.append(fwe(all_m,all_e))
                        for i in range(0,size):
                            subset.append(fwe(all_m[i],all_e[i]))
                        d.y = np.array(subset)
                        d.x =     np.arange(0,len(d.y))
                        d.props['hdf5_path'] = path + m
                        d.props['observable'] = m
                        d.props.update(params)
                        for s in statvar:
                            if s in grp[m].keys():
                                d.props[s] = grp[m+"/"+s+"/value"].value
                        sets.append(d)
                    except AttributeError:
                        pass
        return sets