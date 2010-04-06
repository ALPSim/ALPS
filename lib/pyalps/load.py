# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Bela Bauer <bauerb@phys.ethz.ch>
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
import numpy as np

from dataset import ResultFile
from dataset import DataSet
#from floatwitherror import FloatWithError as fwe

import pyalps.pytools as pt # the C++ conversion functions

# or the C++ class as alternative
from pyalps.pyalea import value_with_error as fwe
from pyalps.pyalea import vector_with_error as vwe
from pyalps.pyalea import vector_of_value_with_error as vfwe
from pyalps.pyalea import convert2vector_of_value_with_error as convert2vfwe

def parse_label(label):
    if '--' in label:
      vals = label.rsplit('--')
      return (eval(vals[0]),eval(vals[1]))
    else:
      return eval(str(label))
 
 
def parse_labels(labels):
    if type(labels)==int: 
      return np.array([labels])
    larr=[]
    allsame = True
    first = None
    for x in labels:
      v = parse_label(x)
      larr.append(v)
      if '--' in x:      
        if first==None:
          first = v[0]
        else:
          if first != v[0]:
            allsame = False
      else:
        allsame = False
    if allsame:
      larr = [x[1] for x in larr]
    return np.array(larr)
       
class Hdf5Missing(Exception):
    def __init__(self,what):
        self.what = what
    
    def __str__(self):
        return 'Failed to find ' + self.what

class Hdf5Loader:
    def GetFileNames(self, flist):
        files = []
        for f in flist:
          if f[-4:]=='.xml':
            f = f[:-3]+'h5'
          else:
            if f[-3:]!='.h5':
              f += '.h5'
          files.append(f)
        return files
        
    # Pre: file is a h5py file descriptor
    # Post: returns a dictionary of all parameters saved in the file
    def ReadParameters(self,proppath):
        LOP = []
        dict = {'filename' : self.h5fname}
        try:
            pgrp = self.h5f.require_group(proppath)
            pgrp.visit(LOP.append)
            for m in LOP:
                try:
                    dict[m] = pgrp[m].value
                    try:
                        dict[m] = float(dict[m])
                    except ValueError:
                        pass
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
        olist = obsgrp.keys()
        return olist

  # Pre: file is a h5py file descriptor
    # Post: returns DataSet with all parameters set
    
    def read_one_spectrum(self,path):
        pass
        
    def ReadSpectrumFromFile(self,flist,proppath='/parameters',respath='/spectrum'):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            fileset=[]
            self.h5f = h5py.File(f)
            self.h5fname = f
            params = self.ReadParameters(proppath)
            grp = self.h5f.require_group(respath)
            if 'energies' in grp.keys():
                    try:
                        d = DataSet()
                        d.props['hdf5_path'] = respath 
                        d.props['observable'] = 'spectrum'
                        d.y = np.array(grp['energies'].value )
                        d.x = range(len(d.y))
                        d.props.update(params)
                        d.props.update(self.ReadParameters('quantumnumbers' ))
                        fileset.append(d)
                    except AttributeError:
                        print "Could not create DataSet"
                        pass
            if 'sectors' in grp.keys():
                sectors_grp = self.h5f.require_group(respath+'/sectors')
                for secnum in sectors_grp.keys():
                    try:
                        d = DataSet()
                        secpath = respath+'/sectors/'+secnum
                        d.props['hdf5_path'] = secpath 
                        d.props['observable'] = 'spectrum'
                        d.y = np.array(sectors_grp[secnum+'/energies'].value )
                        d.x = range(len(d.y))
                        d.props.update(params)
                        d.props.update(self.ReadParameters(secpath+'/quantumnumbers' ))
                        fileset.append(d)
                    except AttributeError:
                        print "Could not create DataSet"
                        pass
            sets.append(fileset)
        return sets

    def ReadDiagDataFromFile(self,flist,proppath='/parameters',respath='/spectrum', measurements=None, index=None):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            fileset=[]
            self.h5f = h5py.File(f)
            self.h5fname = f
            params = self.ReadParameters(proppath)
            grp = self.h5f.require_group(respath)
            if 'results' in grp.keys():
                list_ = self.GetObservableList(respath+'/results')
                if measurements == None:
                    obslist = list_
                else:
                    obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                for m in obslist:
                    if "mean" in grp['results/'+m].keys():
                        try:
                            d = DataSet()
                            secresultspath = respath+'/results/'+m
                            d.props['hdf5_path'] = secresultspath 
                            d.props['observable'] = pt.hdf5_name_decode(m)
                            if index == None:
                                d.y = np.array(grp['results/'+m+'/mean/value'].value )
                                d.x = range(len(d.y))
                            else:
                                try:
                                    d.y = np.array(grp['results/'+m+'/mean/value'].value[index] )
                                except:
                                    pass
                            if "labels" in grp['results/'+m].keys():
                                d.x = parse_labels(grp['results/'+m+'/labels'].value)
                            else:
                                d.x = np.arange(0,len(d.y))
                            d.props.update(params)

                        except AttributeError:
                            print "Could not create DataSet"
                            pass
                    fileset.append(d)
            if 'sectors' in grp.keys():
                sectors_grp = self.h5f.require_group(respath+'/sectors')
                list_ = self.GetObservableList(respath+'/sectors/0/results')
                if measurements == None:
                    obslist = list_
                else:
                    obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                for secnum in sectors_grp.keys():
                    sector_sets=[]
                    for m in obslist:
                        if "mean" in sectors_grp[secnum+'/results/'+m].keys():
                            try:
                                d = DataSet()
                                secpath = respath+'/sectors/'+secnum
                                secresultspath = respath+'/sectors/'+secnum+'/results/'+m
                                d.props['hdf5_path'] = secresultspath 
                                d.props['observable'] = pt.hdf5_name_decode(m)
                                if index == None:
                                    d.y = np.array(sectors_grp[secnum+'/results/'+m+'/mean/value'].value )
                                    d.x = range(len(d.y))
                                else:
                                    try:
                                        d.y = np.array(sectors_grp[secnum+'/results/'+m+'/mean/value'].value[index] )
                                    except:
                                        pass
                                if "labels" in sectors_grp[secnum+'/results/'+m].keys():
                                    d.x = parse_labels(sectors_grp[secnum+'/results/'+m+'/labels'].value)
                                else:
                                    d.x = np.arange(0,len(d.y))
                                d.props.update(params)
                                d.props.update(self.ReadParameters(secpath+'/quantumnumbers'))
                                sector_sets.append(d)

                            except AttributeError:
                                print "Could not create DataSet"
                                pass
                    fileset.append(sector_sets)
            sets.append(fileset)
        return sets
        
    # Pre: file is a h5py file descriptor
    # Post: returns DataSet with the evaluated binning analysis set
    def ReadBinningAnalysis(self,flist,measurements=None,proppath='/parameters',respath=None):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            fileset = []
            self.h5f = h5py.File(f)
            self.h5fname = f
            if respath == None:
              respath="/simulation/results"
            list_ = self.GetObservableList(respath)
            # this is exception-safe in the sense that it's also required in the line above
            grp = self.h5f.require_group(respath)
            params = self.ReadParameters(proppath)
            obslist = []
            if measurements == None:
                obslist = list_
            else:
                obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
            for m in obslist:
                try:
                    d = DataSet()
                    if "timeseries" in grp[m].keys():
                        k = grp[m+'/timeseries'].keys()
                        if "logbinning" in k and "logbinning2" in k and "logbinning_counts" in k:
                            bins = np.array(grp[m+"/timeseries/logbinning"].value[0:-4])
                            bins2 = np.array(grp[m+"/timeseries/logbinning2"].value[0:-4])
                            counts = np.array(grp[m+"/timeseries/logbinning_counts"].value[0:-4])
                            scale = 1
                            for i in range(len(counts)):
                                mean = bins[i]/(counts[i]*scale)
                                mean2 = bins2[i]/counts[i]
                                bins2[i] = np.sqrt((mean2-mean*mean)/counts[i])
                                scale *=2
                            d.y = bins2
                            d.x = np.arange(0,len(d.y))
                            d.props['hdf5_path'] = respath + m
                            d.props['observable'] = 'binning analysis of ' + pt.hdf5_name_decode(m)
                            d.props.update(params)
                            fileset.append(d)
                except AttributeError:
                    print "Could not create DataSet"
                    pass
            sets.append(fileset)
        return sets
                                                                            
    # Pre: file is a h5py file descriptor
    # Post: returns DataSet with all parameters set
    def ReadMeasurementFromFile(self,flist,statvar=None,proppath='/parameters',respath='/simulation/results',measurements=None):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            fileset = []
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
            for m in obslist:
                if not statvar: #if not a specific statistical variable is specified then use the default mean & error or mean
                    try:
                        d = DataSet()
                        if "mean" in grp[m].keys() and "error" in grp[m+"/mean"].keys():
                            mean = grp[m+"/mean/value"].value
                            error = grp[m+"/mean/error"].value
                            d.props['count'] = grp[m+"/count"].value
                            try:
                                size = len(mean)
                                if size == 1:
                                    d.y = np.array([fwe(mean[0],error[0])])
                                else:
                                    d.y = convert2vfwe(vwe(mean,error))
                            except:
                                size=0
                                d.y = np.array([fwe(mean,error)])
                        elif "mean" in grp[m].keys():
                            value = grp[m+"/mean/value"].value
                            try:
                                size=len(value)
                                d.y = np.array([float(x) for x in value])
                            except:
                                size=1
                                d.y = np.array([value])
                        if "labels" in grp[m].keys():
                            d.x = parse_labels(grp[m+"/labels"].value)
                        else:
                            d.x = np.arange(0,len(d.y))
                        d.props['hdf5_path'] = respath + m
                        d.props['observable'] = pt.hdf5_name_decode(m)
                        d.props.update(params)
                        fileset.append(d)
                    except AttributeError:
                        print "Could not create DataSet"
                        pass
                else:
                    for s in statvar:
                        for k in grp[m+"/"+s]:
                            try:
                                d=DataSet()
                                value = np.array(grp[m+"/"+s+"/"+k].value)
                                try:
                                    size = len(value)
                                    if size == 1:
                                        d.y = np.array([value])
                                    else:
                                        d.y = np.array(value)
                                except:
                                    size=1
                                    d.y = np.array([value])
                                
                                if "labels" in grp[m].keys():
                                    d.x = parse_labels(grp[m+"/labels"].value)
                                else:
                                    d.x = np.arange(0,len(d.y))
                                d.props['hdf5_path'] = respath + m
                                d.props['observable'] = pt.hdf5_name_decode(m+"/"+s+"/"+k)
                                d.props.update(params)
                                fileset.append(d)
                            except AttributeError:
                                print "Could not create DataSet"
                                pass
            sets.append(fileset)
        return sets

def loadBinningAnalysis(files,what=None):
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadBinningAnalysis(files,measurements=what)

def loadMeasurements(files,what=None):
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadMeasurementFromFile(files,measurements=what)
    
    
def loadEigenstateMeasurements(files,what=None):
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadDiagDataFromFile(files,measurements=what)

def loadSpectra(files):
    ll = Hdf5Loader()
    return ll.ReadSpectrumFromFile(files)
