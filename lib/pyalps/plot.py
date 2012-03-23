# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Bela Bauer <bauerb@phys.ethz.ch>
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

import matplotlib.pyplot as plt
import numpy as np
from hlist import flatten
from dataset import DataSet
from matplotlib.font_manager import FontProperties
import platform
from plot_core import convertToText, makeGracePlot, makeGnuplotPlot

colors = ['k','b','g','m','c','y']
markers = ['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '8', '+', 'x']


def plot(data):
    """ plots a list of datasets
    
        This function takes a DataSet or a list of DataSets and creates a matplotlib plot.
        
        It creates a new plot set for each dataset, using the x and y members of the DataSet.
        It also inspects the props dictionary and uses the following key-value pairs in that dict:
        
         title
         xlabel
         ylabel
         label
         filename (used as alternative label of the set if no label is specified)
         color (can be 'k','b','g','m','c', or 'y'
         line (can be 'line' or 'scatter')
      
    """
    lines = []
    icolor = 0
    imarker = 0
    if isinstance(data,DataSet):
      s = [data]
    else:
      s = data
    for q in flatten(s):
        try:
            xmeans = np.array([xx.mean for xx in q.x])
            xerrors = np.array([xx.error for xx in q.x])
        except AttributeError:
            xmeans = [float(vvv) for vvv in q.x]
            xerrors = None
        except TypeError:
            xmeans = [q.x]
            xerrors = None
        
        try:
            ymeans = np.array([xx.mean for xx in q.y])
            yerrors = np.array([xx.error for xx in q.y])
        except AttributeError:
            ymeans = [float(vvv) for vvv in q.y]
            yerrors = None
        except TypeError: # this usually means that it's scalar
            ymeans = [q.y]
            yerrors = None

        if 'label' in q.props and q.props['label'] != 'none':
            lab = q.props['label']
        elif 'filename' in q.props:
            lab = q.props['filename']
        else:
            lab = None

        if 'xlabel' in q.props:
            plt.xlabel(q.props['xlabel'])

        if 'ylabel' in q.props:
            plt.ylabel(q.props['ylabel'])

        if 'title' in q.props:
            plt.title(q.props['title'])
            
        thiscolor = colors[icolor]
        icolor = (icolor+1)%len(colors)
        if 'color' in q.props:
            thiscolor = q.props['color']
        
        if 'line' in q.props and q.props['line'] == 'scatter':
            plt.scatter(xmeans, ymeans, c=thiscolor, marker=markers[imarker], label=lab)
            imarker = (imarker+1)%len(markers)
        else:
            line_props = None
            if 'line' in q.props:
                line_props = q.props['line']
            
            plt.errorbar(xmeans,ymeans,yerr=yerrors,xerr=xerrors,fmt=line_props,color=thiscolor,label=lab)




class MplXYPlot_core:
    colors = ['k','b','g','m','c','y']
    markers = ['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '8', '+', 'x']
    
    def __init__(self):
        self.icolor = 0
        self.output = ''
    
    def draw_lines(self):
        self.lines = []
        self.icolor = 0
        self.imarker = 0
        
        xlog = False
        ylog = False
        if 'xaxis' in self.plt and 'logarithmic' in self.plt['xaxis']:
            xlog = self.plt['xaxis']['logarithmic']
        if 'yaxis' in self.plt and 'logarithmic' in self.plt['yaxis']:
            ylog = self.plt['yaxis']['logarithmic']
        
        for q in flatten(self.plt['data']):
            try:
                xmeans = np.array([xx.mean for xx in q.x])
                xerrors = np.array([xx.error for xx in q.x])
            except AttributeError:
                xmeans = [float(vvv) for vvv in q.x]
                xerrors = None
            
            try:
                ymeans = np.array([xx.mean for xx in q.y])
                yerrors = np.array([xx.error for xx in q.y])
            except AttributeError:
                ymeans = [float(vvv) for vvv in q.y]
                yerrors = None
            
            thiscolor = self.colors[self.icolor]
            self.icolor = (self.icolor+1)%len(self.colors)
            if 'color' in q.props:
                thiscolor = q.props['color']
            
            if 'line' in q.props and q.props['line'] == 'scatter':
                thismarker = self.markers[self.imarker]
                if 'marker' in q.props:
                    thismarker = q.props['marker']
                self.imarker = (self.imarker+1)%len(self.markers)
                self.lines.append([plt.scatter(xmeans, ymeans, c=thiscolor, marker=thismarker)])
                if xerrors != None or yerrors != None:
                    plt.errorbar(xmeans, ymeans, yerr=yerrors, xerr=xerrors, fmt=None)
            else:
                line_props = thiscolor
                if 'line' in q.props:
                    line_props += q.props['line']
                
                self.lines.append(plt.errorbar(xmeans,ymeans,yerr=yerrors,xerr=xerrors,fmt=line_props))
                
                if 'linewidth' in q.props:
                    self.lines[-1][0].set_linewidth(q.props['linewidth'])
                if 'markersize' in q.props:
                    self.lines[-1][0].set_markersize(q.props['markersize'])
                
            if xlog:
                plt.xscale('log')
            if ylog:
                plt.yscale('log')
            
            if 'label' in q.props and q.props['label'] != 'none':
                self.lines[-1][0].set_label(q.props['label'])
            elif 'filename' in q.props:
                self.lines[-1][0].set_label(q.props['filename'])
            
            
            if 'legend' in self.plt:
                if 'scatter_labels' in self.plt['legend']:
                    if self.plt['legend']['scatter_labels'] == True:
                        plt.annotate(q.props['label'], (xmeans[0],ymeans[0]))
            
    def __call__(self, desc):
        self.plt = desc
        
        self.draw_lines()
        
        for ds in flatten(self.plt['data']):
            if 'xlabel' in ds.props:
                plt.xlabel(ds.props['xlabel'])
            if 'ylabel' in ds.props:
                plt.ylabel(ds.props['ylabel'])
        
        if 'xaxis' in self.plt:
            if 'label' in self.plt['xaxis']:
                if 'fontsize' in self.plt['xaxis']:
                    plt.xlabel(self.plt['xaxis']['label'],fontsize=self.plt['xaxis']['fontsize'])
                else:
                    plt.xlabel(self.plt['xaxis']['label'])
            if 'min' in self.plt['xaxis'] and 'max' in self.plt['xaxis']:
                plt.xlim(self.plt['xaxis']['min'],self.plt['xaxis']['max'])
                
        if 'yaxis' in self.plt:
            if 'label' in self.plt['yaxis']:
                if 'fontsize' in self.plt['yaxis']:
                    plt.ylabel(self.plt['yaxis']['label'],fontsize=self.plt['yaxis']['fontsize'])
                else:
                    plt.ylabel(self.plt['yaxis']['label'])
            if 'min' in self.plt['yaxis'] and 'max' in self.plt['yaxis']:
                plt.ylim(self.plt['yaxis']['min'],self.plt['yaxis']['max'])
        
        if 'legend' in self.plt:
            showlegend = True
            if 'scatter_labels' in self.plt['legend']:
                if self.plt['legend']['scatter_labels'] == True:
                    showlegend = False
            if showlegend:
                legend_prop = {}
                if platform.system() != 'Linux' and 'fontsize' in self.plt['legend']:
                    legend_prop['size'] = self.plt['legend']['fontsize']
                if 'location' in self.plt['legend']:
                    legend_loc = self.plt['legend']['location']
                else:
                    legend_loc = 0

                if len(legend_prop) > 0:
                    plt.legend(loc=legend_loc,prop=FontProperties(legend_prop))
                else:
                    plt.legend(loc=legend_loc)
        
        if 'title' in self.plt:
            plt.title(self.plt['title'])
        
        if 'grid' in self.plt:
            plt.grid(self.plt['grid'])
        
            
 
 
 
