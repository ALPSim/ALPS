# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2010 by Bela Bauer <bauerb@phys.ethz.ch>
#                            Brigitte Surer <surerb@phys.ethz.ch> 
# 
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

import urllib, copy, h5py, os
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import datetime

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
        self.files = [f[:-3]+"h5" for f in flist if f[-3:]=='xml']
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
            pass
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
    def ReadSpectrumFromFile(self,flist,proppath,respath):
        print "I'm here spectrum", flist
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            print "Reading ", f
            self.h5f = h5py.File(f)
            self.h5fname = f
            params = self.ReadParameters(proppath)
            print params
            grp = self.h5f.require_group(respath)
            print 'path ', respath, ' has keys ',grp.keys()
            if 'sectors' in grp.keys():
                print "Has sectors"
                sectors_grp = self.h5f.require_group(respath+'/sectors')
                for secnum in sectors_grp.keys():
                    try:
                        d = DataSet()
                        secpath = respath+'/sectors/'+secnum
                        d.props['hdf5_path'] = secpath 
                        d.y = np.array(sectors_grp[secnum+'/energies'].value )
                        d.props.update(params)
                        d.props.update(self.ReadParameters(secpath+'/quantumnumbers' ))
                        sets.append(d)
                    except AttributeError:
                        print "Could not create DataSet"
                        pass
        return sets
        
    # Pre: file is a h5py file descriptor
    # Post: returns DataSet with all parameters set
    def ReadMeasurementFromFile(self,flist,statvar,proppath,respath,measurements=None):
        print "I'm here", len(flist)
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
                            d.props["timeseries_"+l] = np.array(grp[m+"/timeseries/"+l].value)
                    if "jacknife" in statvar:
                        jklist = grp[m+"/jacknife"].keys()
                        for l in jklist:
                            d.props["jacknife_"+l] = np.array(grp[m+"/jacknife/"+l].value)
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

 