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

import urllib, copy, os, traceback
import numpy as np

import pyalps.hdf5 as h5
import pyalps.alea as pa

from dataset import ResultFile
from dataset import DataSet
#from floatwitherror import FloatWithError as fwe

import pyalps.pytools as pt # the C++ conversion functions

# or the C++ class as alternative
from pyalps.alea import MCScalarData as fwe
from pyalps.alea import MCVectorData as vwe

in_vistrails=True
try:
  import vistrails.core.modules.basic_modules
  from core import debug
except:
  in_vistrails=False

def log(m):
    """ print a log message either to the console or through the VisTrails logger """
    if in_vistrails:
      debug.log(str(m))
    else:
      print m
    
def parse_label(label):
    if '--' in label:
      vals = label.rsplit('--')
      ret = ()
      for val in vals:
          ret = ret + (eval(val),)
      return ret
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
          if first != v[0] or len(v) != 2:
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
    """The Hdf5Loader class loads simulation parameters and observables from hdf5-files and returns them as hierarchical datasets"""
    def GetFileNames(self, flist):
        files = []
        for f in flist:
          if f[-4:]=='.xml':
            f = f[:-3]+'h5'
          else:
            if f[-3:]!='.h5':
              f += '.h5'
          if os.path.exists(f):
            files.append(f)
          else:
            log( "FILE "+ f+ "DOES NOT EXIST!")
        return files
        
    def ReadParameters(self,proppath):
        dict = {'filename' : self.h5fname}
        LOP=self.h5f.list_children(proppath)
        for m in LOP:
                try:
                    dict[m] = self.h5f[proppath+'/'+m]
                    try:
                        dict[m] = float(dict[m])
                    except:
                        dict[m] = map(float,dict[m])
                except ValueError:
                    pass
        return dict 
        
    def GetProperties(self,flist,proppath='/parameters',respath='/simulation/results',verbose=False):
        fs = self.GetFileNames(flist)
        resultfiles = []
        for f in fs:
            try:
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log( "Loading from file" + f)
                rfile = ResultFile(f)
                rfile.props = self.ReadParameters(proppath)
                try:
                    obs = self.GetObservableList(respath)
                    rfile.props["ObservableList"] = [pt.hdf5_name_decode(x) for x in obs]
                except: pass
                resultfiles.append(rfile)
            except Exception, e:
                log(e)
                log(traceback.format_exc())
        return resultfiles
        
    def GetObservableList(self,respath):
        if self.h5f.is_group(respath):
            olist = self.h5f.list_children(respath)
        else:
            olist = []
        return olist

# Pre: file is a hdf5 file descriptor
# Post: returns DataSet with all parameters set
    
    def read_one_spectrum(self,path):
        pass
        
    def ReadSpectrumFromFile(self,flist,proppath='/parameters',respath='/spectrum',verbose=False):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            try:
                fileset=[]
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log("Loading from file " + f)
                params = self.ReadParameters(proppath)
                if 'energies' in self.h5f.list_children(respath):
                        try:
                            d = DataSet()
                            d.props['hdf5_path'] = respath 
                            d.props['observable'] = 'spectrum'
                            d.y = self.h5f[respath+'/energies']
                            d.x = range(len(d.y))
                            d.props.update(params)
                            try:
                                d.props.update(self.ReadParameters('quantumnumbers'))
                            except:
                                if verbose: log("no quantumnumbers stored ")
                                pass
                            fileset.append(d)
                        except AttributeError:
                            pass
                if 'sectors' in self.h5f.list_children(respath):
                    for secnum in self.h5f.list_children(respath+'/sectors'):
                        try:
                            d = DataSet()
                            secpath = respath+'/sectors/'+secnum
                            d.props['hdf5_path'] = secpath 
                            d.props['observable'] = 'spectrum'
                            d.y = self.h5f[secpath+'/energies']
                            d.x = range(len(d.y))
                            d.props.update(params)
                            try:
                                d.props.update(self.ReadParameters(secpath+'/quantumnumbers'))
                            except:
                                if verbose: log("no quantumnumbers stored ")
                                pass
                            fileset.append(d)
                        except AttributeError:
                            log( "Could not create DataSet")
                            pass
                sets.append(fileset)
            except Exception, e:
                log(e)
                log(traceback.format_exc())
        return sets
        
    def GetIterations(self, current_path, params={}, measurements=None, index=None, verbose=False):
        iterationset=[]
        #iteration_grp = self.h5f.require_group(respath+'/iteration')
        for it in self.h5f.list_children(current_path+'/iteration'):
            obsset=[]
            iteration_props = {}
            if 'parameters' in self.h5f.list_children(current_path+'/iteration/'+it):
                iteration_props = self.ReadParameters(current_path+'/iteration/'+it+'/parameters')
            iteration_props['iteration'] = it
            
            respath = current_path+'/iteration/'+it+'/results'
            list_ = self.GetObservableList(respath)
            if measurements == None:
                obslist = list_
            else:
                obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
            for m in obslist:
                if m in self.h5f.list_children(respath):
                    if "mean" in self.h5f.list_children(respath+'/'+m):
                        try:
                            d = DataSet()
                            itresultspath = respath+'/'+m
                            if verbose: log("Loading "+ m)
                            measurements_props = {}
                            measurements_props['hdf5_path'] = itresultspath 
                            measurements_props['observable'] = pt.hdf5_name_decode(m)
                            if index == None:
                                d.y = self.h5f[itresultspath+'/mean/value']
                                d.x = np.arange(0,len(d.y))
                            else:
                                try:
                                    d.y = self.h5f[itresultspath+'/mean/value'][index]
                                except:
                                    pass
                            if "labels" in self.h5f.list_children(itresultspath):
                                d.x = parse_labels(self.h5f[itresultspath+'/labels'])
                            else:
                                d.x = np.arange(0,len(d.y))
                            d.props.update(params)
                            d.props.update(iteration_props)
                            d.props.update(measurements_props)
                        except AttributeError:
                            log( "Could not create DataSet")
                    obsset.append(d)
            iterationset.append(obsset)
        return iterationset
        
    def ReadDiagDataFromFile(self,flist,proppath='/parameters',respath='/spectrum', measurements=None, index=None, loadIterations=False,verbose=False):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            try:
                fileset=[]
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log("Loading from file"+ f)
                params = self.ReadParameters(proppath)
                if 'results' in self.h5f.list_children(respath):
                    list_ = self.GetObservableList(respath+'/results')
                    if measurements == None:
                        obslist = list_
                    else:
                        obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                    if loadIterations==True:
                        if "iteration" in self.h5f.list_children(respath):
                            fileset.append(self.GetIterations(respath+'/results', params, measurements, index, verbose))
                    else:        
                        for m in obslist:
                            if "mean" in self.h5f.list_children(respath+'/results/'+m):
                                try:
                                    if verbose: log("Loading" + m)
                                    d = DataSet()
                                    secresultspath = respath+'/results/'+m
                                    d.props['hdf5_path'] = secresultspath 
                                    d.props['observable'] = pt.hdf5_name_decode(m)
                                    if index == None:
                                        d.y = self.h5f[secresultspath+'/mean/value']
                                        d.x = np.arange(0,len(d.y))
                                    else:
                                        try:
                                            d.y = self.h5f[secresultspath+'/mean/value'][index]
                                        except:
                                            pass
                                    if "labels" in self.h5f.list_children(secresultspath):
                                        d.x = parse_labels(self.h5f[secresultspath+'/labels'])
                                    else:
                                        d.x = np.arange(0,len(d.y))
                                    d.props.update(params)

                                except AttributeError:
                                    log("Could not create DataSet")
                            fileset.append(d)
                if loadIterations==True:
                    if "iteration" in self.h5f.list_children(respath):
                        fileset.append(self.GetIterations(respath, params, measurements, index, verbose))
                if 'sectors' in self.h5f.list_children(respath):
                    list_ = self.GetObservableList(respath+'/sectors/0/results')
                    if measurements == None:
                        obslist = list_
                    else:
                        obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                    for secnum in self.h5f.list_children(respath+'/sectors'):
                        sector_sets=[]
                        for m in obslist:
                            if "mean" in self.h5f.list_children(respath+'/sectors/'+secnum+'/results/'+m):
                                try:
                                    if verbose: log("Loading" + m)
                                    d = DataSet()
                                    secpath = respath+'/sectors/'+secnum
                                    secresultspath = respath+'/sectors/'+secnum+'/results/'+m
                                    d.props['hdf5_path'] = secresultspath 
                                    d.props['observable'] = pt.hdf5_name_decode(m)
                                    if index == None:
                                        d.y = self.h5f[secresultspath+'/mean/value']
                                        d.x = np.arange(0,len(d.y))
                                    else:
                                        try:
                                            d.y = self.h5f[secresultspath+'/mean/value'][index]
                                        except:
                                            pass
                                    if "labels" in self.h5f.list_children(secresultspath):
                                        d.x = parse_labels(self.h5f[secresultspath+'/labels'])
                                    else:
                                        d.x = np.arange(0,len(d.y))
                                    d.props.update(params)
                                    try:
                                        d.props.update(self.ReadParameters(secpath+'/quantumnumbers'))
                                    except:
                                        if verbose: log("no quantumnumbers stored ")
                                        pass
                                    sector_sets.append(d)

                                except AttributeError:
                                    log( "Could not create DataSet")
                                    pass
                        fileset.append(sector_sets)
                sets.append(fileset)
            except RuntimeError:
                raise
            except Exception, e:
                log(e)
                log(traceback.format_exc())
        return sets
        
    # Pre: file is a hdf5 file descriptor
    # Post: returns DataSet with the evaluated binning analysis set
    def ReadBinningAnalysis(self,flist,measurements=None,proppath='/parameters',respath=None,verbose=False):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            try:
                fileset = []
                if verbose: log( 'loading from file ' +f)
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if respath == None:
                  respath="/simulation/results"
                list_ = self.GetObservableList(respath)
                # this is exception-safe in the sense that it's also required in the line above
                #grp = self.h5f.require_group(respath)
                params = self.ReadParameters(proppath)
                obslist = []
                if measurements == None:
                    obslist = list_
                else:
                    obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                for m in obslist:
                    try:
                        d = DataSet()
                        if "timeseries" in  self.h5f.list_children(respath+'/'+m):
                            k = self.h5f.list_children(respath+'/'+m+'/timeseries')
                            if "logbinning" in k and "logbinning2" in k and "logbinning_counts" in k:
                                if verbose: log("Loading"+ m)
                                bins = self.h5f[respath+'/'+m+'/timeseries/logbinning'][0:-4]
                                bins2 = self.h5f[respath+'/'+m+'/timeseries/logbinning2'][0:-4]
                                counts = self.h5f[respath+'/'+m+'/timeseries/logbinning_counts'][0:-4]
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
                                if verbose: log( '  loaded binnig analysis for '+m)
                                fileset.append(d)
                    except AttributeError:
                        log( "Could not create DataSet")
                sets.append(fileset)
            except Exception, e:
                log( e)
                log( traceback.format_exc())
        return sets
    
    # Pre: file is a hdf5 file descriptor
    # Post: returns DataSet with all parameters set
    def ReadMeasurementFromFile(self,flist,proppath='/parameters',respath='/simulation/results',measurements=None,verbose=False):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            try:
                fileset = []
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log("Loading from file " + f)
                list_ = self.GetObservableList(respath)
                params = self.ReadParameters(proppath)
                obslist = []
                if measurements == None:
                    obslist = list_
                else:
                    obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                for m in obslist:
                    if verbose: log( "Loading " + m)
                    size=0
                    xmin=0
                    xstep=1
                    if "histogram" in self.h5f.list_children(respath+'/'+m):
                        obs = self.h5f[respath+'/'+m+'/histogram']
                        xmin = self.h5f[respath+'/'+m+'/@min']
                        xstep = self.h5f[respath+'/'+m+'/@stepsize']
                        size = len(obs)
                    elif "error" in self.h5f.list_children(respath+'/'+m+'/mean'): 
                        if self.h5f.is_scalar(respath+'/'+m+'/mean/value'):
                            obs = pa.MCScalarData()
                            obs.load(self.h5fname, respath+'/'+m)
                            obs=np.array([obs])
                            size=1
                            if obs[0].count==0:
                              obs=None
                        else:
                            obs=None
                            if not self.h5f.is_group(respath+'/'+m+'/timeseries'): # check for simple binning
                                obs = np.array(self.h5f[respath+'/'+m+'/mean/value']);
                                if params.has_key('L'): # ugly fix... really ugly
                                  L = int(params['L'])
                                  if L == obs.size:
                                    params['origin'] = [(L-1.)/2.];
                                  if L**2 == obs.size: # dimension 2
                                    obs = obs.reshape([L,L]);
                                    params['origin'] = [(L-1.)/2., (L-1.)/2.];
                                  elif L**3 == obs.size: # dimension 3
                                    obs = obs.reshape([L,L,L]);
                                    params['origin'] = [(L-1.)/2., (L-1.)/2., (L-1.)/2.];
                                size = obs.size
                            else:
                                obs = pa.MCVectorData()
                                obs.load(self.h5fname, respath+'/'+m)
                                size=len(obs.mean)
                                if obs.count==0:
                                    obs=None
                    else:
                        if self.h5f.is_scalar(respath+'/'+m+'/mean/value'):
                            obs = self.h5f[respath+'/'+m+'/mean/value']
                            obs=np.array([obs])
                            size=1
                        else:
                            obs = self.h5f[respath+'/'+m+'/mean/value']
                            size=len(obs)
                    try:
                      if obs is not None:
                        d = DataSet()
                        d.y = obs
                        d.x = np.arange(xmin,xmin+xstep*size,xstep)
                        d.props['hdf5_path'] = respath +"/"+ m
                        d.props['observable'] = pt.hdf5_name_decode(m)
                        d.props.update(params)
                        fileset.append(d)
                    except AttributeError:
                        log( "Could not create DataSet")
                sets.append(fileset)
            except Exception, e:
                log(e)
                log(traceback.format_exc())
        return sets

    # Pre: file is a hdf5 file descriptor
    # Post: returns DataSet with all parameters set
    def ReadDMFTIterations(self,flist,observable='G_tau',measurements='0',proppath='/parameters',respath='/simulation/iteration',verbose=False):
        fs = self.GetFileNames(flist)
        fileset = []
        for f in fs:
            try:
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log("Loading from file "+ f)
                list_ = self.GetObservableList(respath+'/1/results/'+observable+'/')
                #grp = self.h5f.require_group(respath)
                params = self.ReadParameters(proppath)
                obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                iterationset=[]
                for it in  self.h5f.list_children(respath):
                    obsset=[]
                    for m in obslist:
                        try:
                            if verbose: log( "Loading "+ m)
                            d = DataSet()
                            size=0
                            path=it+'/results/'+observable+'/'+m
                            if "mean" in self.h5f.list_children(respath+'/'+path):
                                if self.h5f.is_scalar(respath+'/'+path+'/mean/value'):
                                    size=1
                                    obs=self.h5f[respath+'/'+path+'/mean/value']
                                    d.y = np.array([obs]) 
                                else:
                                    obs=self.h5f[respath+'/'+path+'/mean/value']
                                    size=len(obs)
                                    d.y = obs
                                d.x = np.arange(0,size)
                            d.props['hdf5_path'] = respath +"/"+ path
                            d.props['observable'] = pt.hdf5_name_decode(m)
                            d.props['iteration'] = it
                            d.props.update(params)
                        except AttributeError:
                            log( "Could not create DataSet")
                            pass
                        obsset.append(d)
                    iterationset.append(obsset)
                fileset.append(iterationset)
            except Exception, e:
                log( e)
                log( traceback.format_exc())
        return fileset
        
def loadBinningAnalysis(files,what=None,verbose=False):
    """ loads MC binning analysis from ALPS HDF5 result files
    
        this function loads results of a MC binning analysis from ALPS HDF5 result files
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            what (list): optional argument that is either a string or list of strings, specifying the names of the observables for which the binning analysis should be loaded
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            a list of list of DataSet objects: loaded binning analysis.
            The elements of the outer list each correspond to the file names specified as input.
            The elements of the inner list are each for a different observable.
            The x-values of the DataSet objects are the logarithmic binning level and the y-values the error estimates at that binning level.
    """
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadBinningAnalysis(files,measurements=what,verbose=verbose)

def loadMeasurements(files,what=None,verbose=False,respath='/simulation/results'):
    """ loads ALPS measurements from ALPS HDF5 result files
    
        this function loads results of ALPS simulations ALPS HDF5 result files
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            what (list): optional argument that is either a string or list of strings, specifying the names of the observables which should be loaded
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            a list of list of DataSet objects: loaded measurements.
            The elements of the outer list each correspond to the file names specified as input.
            The elements of the inner list are each for a different observable.
            The y-values of the DataSet objects are the measurements and the x-values optionally the labels (indices) of array-valued measurements
    """
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadMeasurementFromFile(files,measurements=what,verbose=verbose,respath=respath)
    
    
def loadEigenstateMeasurements(files,what=None, verbose=False):
    """ loads ALPS eigenstate measurements from ALPS HDF5 result files
    
        this function loads results of ALPS diagonalization or DMRG simulations from an HDF5 file
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            what (list): an optional argument that is either a string or list of strings, specifying the names of the observables which should be loaded
            verbose (bool): an optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            list of list of (lists of) DataSet objects: loaded measurements.
            The elements of the outer list each correspond to the file names specified as input
            The elements of the next level are different quantum number sectors, if any exists
            The elements of the inner-most list are each for a different observable
            The y-values of the DataSet objects is an array of the measurements in all eigenstates calculated in this sector, and the x-values optionally the labels (indices) of array-valued measurements
    """
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadDiagDataFromFile(files,measurements=what,verbose=verbose)
    
def loadIterationMeasurements(files,what=None,verbose=False):
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadDiagDataFromFile(files,measurements=what,loadIterations=True,verbose=verbose)

def loadMPSIterations(files,what=None,verbose=False):
    ## MPS iterations are stored in: /simulation/iteration/X/results/M with X the sweep index and M the name of the measurement
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadDiagDataFromFile(files,measurements=what,respath='/simulation',loadIterations=True,verbose=verbose)

def loadSpectra(files,verbose=False):
    """ loads ALPS spectra from ALPS HDF5 result files
    
        This function loads the spectra calculated in ALPS diagonalization or DMRG simulations from an HDF5 file.
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded.
        
        Returns:
            list of (lists of) DataSet objects: Loaded spectra.
            The elements of the outer list each correspond to the file names specified as input.
            The elements of the next level are different quantum number sectors, if any exists.
            The y-values of the DataSet objects are the energies in that quantum number sector.
    """
    ll = Hdf5Loader()
    return ll.ReadSpectrumFromFile(files,verbose=verbose)

def loadDMFTIterations(files,observable='G_tau',measurements='0',verbose=False):
    """ loads ALPS measurements from ALPS HDF5 result files
    
        this function loads results of ALPS simulations ALPS HDF5 result files
        
        Parameters:
            files (list): ALPS HDF5 result files.
            observable (str): optional argument specifying the name of the observables which should be loaded
            measurements (list): optional argument that is either a string or list of strings, specifying the names of the measurements which should be loaded
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            list of list of list of DataSet objects: loaded iteration measurements.
            The elements of the outer list each correspond to the file names specified as input.
            The elements of the next level are different iterations.
            The elements of the inner list contains a DataSet for each measurement.
            The y-values of the DataSet objects are the measurements and the x-values optionally the labels (indices) of array-valued measurements
    """
    ll = Hdf5Loader()
    if isinstance(measurements,str):
      measurements = [measurements]
    return ll.ReadDMFTIterations(files,observable=observable,measurements=measurements,verbose=verbose)

def loadProperties(files,proppath='/parameters',respath='/simulation/results',verbose=False):
    """ loads properties (parameters) of simulations from ALPS HDF5 result files
    
        this function loads the properties (parameters) of ALPS simulations ALPS HDF5 result files
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            list of dicts: properties contained in each file.
    """
    ll = Hdf5Loader()
    res = ll.GetProperties(files,proppath,respath,verbose=verbose)
    results = []
    for x in res:
      results.append(x.props)
    return results

def loadObservableList(files,proppath='/parameters',respath='/simulation/results',verbose=False):   
    """ loads lists of existing measurements from ALPS HDF5 result files
    
        The function returns a list of lists, containing the names of measurements that are stored in the result files 
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
    """
    ll = Hdf5Loader()
    res = ll.GetProperties(files,proppath,respath,verbose=verbose)
    results = []
    for x in res:
      results.append(x.props['ObservableList'])
    return results

def loadTimeEvolution( flist,globalproppath='/parameters',resroot='/timesteps/',localpropsuffix='/parameters', measurements=None):
    ll=Hdf5Loader()
    data=[]
    #loop over files
    for f in flist:
        try:
            #open the file and open the results root group
            h5file = h5.archive(f, 'r')
            #enumerate the subgroups
            L=h5file.list_children(resroot)
            #Create an iterator of length the number of subgroups
            stepper=[i+1 for i in range(len(L))]
            #Read in global props
            globalprops=ll.GetProperties([f],globalproppath)
            for d in stepper:
                #Get the measurements from the numbered subgroups
                locdata=ll.ReadMeasurementFromFile([f],proppath=resroot+str(d)+localpropsuffix, \
                respath=resroot+str(d)+'/results', measurements=measurements)
                #Append the global props to the local props
                for i in range(len(locdata[0])):
                    locdata[0][i].props.update(globalprops[0].props)
                #Extend the total dataset with this data
                data.extend(locdata)
        except Exception, e:
            log(e)
            log( traceback.format_exc())
    return data
    

   
