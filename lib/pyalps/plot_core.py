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

import matplotlib.pyplot as plt
import numpy as np
from xml.etree import ElementTree
from dataset import DataSet
from floatwitherror import FloatWithError as fwe
from floatwitherror import get_mean
from hlist import flatten

def read_xml(filename):
    root = ElementTree.parse(filename).getroot()
    data = DataSet()

    data.props['xlabel'] = root.find('xaxis').attrib['label']
    data.props['ylabel'] = root.find('yaxis').attrib['label']

    x = []
    y = []
    dx = []
    dy = []
    for point in root.find('set').getchildren():
        x.append(float(point.find('x').text))
        y.append(float(point.find('y').text))
        if point.find('dx') != None:
            x[-1] = fwe(x[-1],float(point.find('dx').text))
        if point.find('dy') != None:
            y[-1] = fwe(y[-1],float(point.find('dy').text))

    data.x = np.array(x)
    data.y = np.array(y)
    
    parameters = root.find('PARAMETERS')
    for par in parameters.getchildren():
        data.props[par.attrib['name']] = par.text
    
    return data

def Axis(label=None,mmin=None,mmax=None,log=False):
    d = {}
    if label != None:
        d['label'] = label
    if mmin != None:
        d['min'] = min
    if mmax != None:
        d['max'] = max
    if log != None:
        d['log'] = log
    return d

def Legend(location=None):
    d = {}
    if location != None:
        d['location'] = location
    return d

def Plot(data,xaxis=None,yaxis=None,legend=None):
    d = {'data':data}
    if xaxis != None:
        d['xaxis'] = xaxis
    if yaxis != None:
        d['yaxis'] = yaxis
    if legend != None:
        d['legend'] = legend

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
                
            if 'line' in q.props and q.props['line'] == 'scatter':
                self.lines.append([plt.scatter(xmeans, ymeans, c=self.colors[self.icolor], marker=self.markers[self.imarker])])
                self.imarker = (self.imarker+1)%len(self.markers)
            else:
                line_props = self.colors[self.icolor]
                if 'line' in q.props:
                    line_props += q.props['line']
                
                self.lines.append(plt.errorbar(xmeans,ymeans,yerr=yerrors,xerr=xerrors,fmt=line_props))
                
            if xlog:
                plt.xscale('log')
            if ylog:
                plt.yscale('log')
            
            if 'label' in q.props and q.props['label'] != 'none':
                self.lines[-1][0].set_label(q.props['label'])
            elif 'filename' in q.props:
                self.lines[-1][0].set_label(q.props['filename'])
            
            self.icolor = (self.icolor+1)%len(self.colors)
            
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
                prop = {}
                if 'fontsize' in self.plt['legend']:
                    prop['size'] = self.plt['legend']['fontsize']
                if 'location' in self.plt['legend']:
                    plt.legend(loc=self.plt['legend']['location'],prop=prop)
                else:
                    plt.legend(prop=prop)
        
        if 'title' in self.plt:
            plt.title(self.plt['title'])
            
 
def convert_to_text(desc):
        output = ''

        if 'title' in desc:
            output += desc['title'] + '\n'           

        if 'xaxis' in desc:
            output += 'X'
            if 'label' in desc['xaxis']:
                output += ': ' + desc['xaxis']['label']
            if 'min' in desc['xaxis'] and 'max' in desc['xaxis']:
                output += ': ' + str(desc['xaxis']['min']) + ' to ' + str(data['xaxis']['max'])
            output+='\n'

        if 'yaxis' in desc:
            output += 'Y'
            if 'label' in desc['yaxis']:
                output += ': ' + desc['yaxis']['label']
            if 'min' in desc['yaxis'] and 'max' in desc['yaxis']:
                output += ': ' + str(desc['yaxis']['min']) + ' to ' + str(data['yaxis']['max'])
            output+='\n\n'
                
        
        for q in desc['data']:
            if 'label' in q.props and q.props['label'] != 'none':
                output += q.props['label']
            elif 'filename' in q.props:
                output += q.props['filename']
            output += '\n'

            for i in range(len(q.x)):
                output += str(q.x[i]) + '\t' + str(q.y[i]) + '\n'
            output+='\n\n'                
        return output
            
def convert_to_grace(desc):
        output =  '# Grace project file\n'
        output += '#\n@    g0 on\n@    with g0\n'
        output += '@     frame linewidth 2.0\n'
        output += '@    page background fill off\n'

        xrange = [0,1]
        yrange = [0,1]
        if 'xaxis' in desc and 'min' in desc['xaxis'] and 'max' in desc['xaxis']: 
            xrange = [ desc['xaxis']['min'],desc['xaxis']['max']]
        if 'yaxis' in desc and 'min' in desc['yaxis'] and 'max' in desc['yaxis']:
            yrange = [ desc['yaxis']['min'],desc['yaxis']['max']]

        output += '@    world ' + str(xrange[0])+', ' + str (yrange[0]) + ','
        output +=                 str(xrange[1])+', ' + str (yrange[1]) + '\n'

        if 'title' in desc:
            output += '@    title "'+ desc['title'] + '"\n'           
            output += '@    title size 1.500000\n'

        xlog = False
        ylog = False
        if 'xaxis' in desc and 'logarithmic' in desc['xaxis']:
            xlog = desc['xaxis']['logarithmic']
        if 'yaxis' in desc and 'logarithmic' in desc['yaxis']:
            ylog = desc['yaxis']['logarithmic']
            
        if xlog:
            output += '@    xaxes scale Logarithmic\n'
        else:
            output += '@    xaxes scale Normal\n'

        if ylog:
            output += '@    yaxes scale Logarithmic\n'
        else:
            output += '@    yaxes scale Normal\n'

        if 'xaxis' in desc:
            if 'label' in desc['xaxis']:
                output += '@    xaxis  label "' + desc['xaxis']['label'] +'"\n'
                output += '@    xaxis  label char size 1.500000\n'
        output += '@    xaxis  ticklabel char size 1.250000\n'
        output += '@    xaxis  tick minor ticks 4\n'

        if 'yaxis' in desc:
            if 'label' in desc['yaxis']:
                output += '@    yaxis  label "' + desc['yaxis']['label'] +'"\n'
                output += '@    yaxis  label char size 1.500000\n'
        output += '@    yaxis  ticklabel char size 1.250000\n'
        output += '@    yaxis  tick minor ticks 4\n'
        
        if 'legend' in desc:
            output += '@    legend on\n'
            output += '@    legend loctype view\n'
            output += '@    legend 0.85, 0.8\n'
        
        num = 0
        symnum = 0
        for q in desc['data']:
            output += '@target G0.S'+str(num)+'\n'
            output += '@    s'+str(num)+' symbol ' + str(num+1) +'\n'
            output += '@    s'+str(num)+' symbol size 0.500000\n'
            if 'line' in q.props and q.props['line'] == 'scatter':
              symnum += 1
              output += '@    s'+str(num)+' line type 0\n'
              output += '@    s'+str(num)+' symbol ' + str(symnum) + '\n'
              output += '@    s'+str(num)+' symbol size 1.000000\n'
            else:
              output += '@    s'+str(num)+' line type 1\n'
            if 'label' in q.props and q.props['label'] != 'none':
                output += '@    s'+str(num)+' legend "' + q.props['label'] + '"\n'
            elif 'filename' in q.props:
                output += '@    s'+str(num)+' legend "' + q.props['filename'] + '"\n'
            output += '\n'

            if len(q.y):
                try:
                    xerrors = np.array([xx.error for xx in q.x])
                except AttributeError:
                    xerrors = None
                
                try:
                    yerrors = np.array([xx.error for xx in q.y])
                except AttributeError:
                    yerrors = None
                    
                if xerrors == None and yerrors == None:
                    output += '@type xy\n'
                    for i in range(len(q.x)):
                        output += str(q.x[i]) + '\t' + str(q.y[i]) + '\n'
                if xerrors == None and yerrors != None:
                    output += '@type xydy\n'
                    for i in range(len(q.x)):
                        output += str(q.x[i]) + '\t' + str(q.y[i].mean) + '\t' + str(q.y[i].error) + '\n'
                if xerrors != None and yerrors == None:
                    output += '@type xydx\n'
                    for i in range(len(q.x)):
                        output += str(q.x[i]) + '\t' + str(q.y[i].mean) + '\t' + str(q.x[i].error) + '\n'
                if xerrors != None and yerrors != None:
                    output += '@type xydxdy\n'
                    for i in range(len(q.x)):
                        output += str(q.x[i]) + '\t' + str(q.y[i].mean) + '\t' + str(q.x[i].error) + '\t' + str(q.x[i].error) + '\n'
                output += '&\n'
                num+=1
                     
        return output
        
def convert_to_gnuplot(desc, outfile="output.eps", fontsize=24):
    output =  '# Gnuplot project file\n'
    output += 'set terminal postscript color eps enhanced '+str(fontsize)+'\n'
    output += 'set output "' + outfile + '"\n'

    if 'xaxis' in desc and 'min' in desc['xaxis'] and 'max' in desc['xaxis']: 
        xrange = [ desc['xaxis']['min'],desc['xaxis']['max']]
        output += 'set xrange [' + str(xrange[0])+': ' + str (xrange[1]) + ']\n'
    if 'yaxis' in desc and 'min' in desc['yaxis'] and 'max' in desc['yaxis']:
        yrange = [ desc['yaxis']['min'],desc['yaxis']['max']]
        output += 'set yrange [' + str(yrange[0])+': ' + str (yrange[1]) + ']\n'
    
    if 'title' in desc:
        output += 'set title "'+ desc['title'] + '"\n'           

    xlog = False
    ylog = False
    if 'xaxis' in desc and 'logarithmic' in desc['xaxis']:
        xlog = desc['xaxis']['logarithmic']
    if 'yaxis' in desc and 'logarithmic' in desc['yaxis']:
        ylog = desc['yaxis']['logarithmic']
        
    if xlog:
        output += 'set xlogscale \n'
    else:
        output += '# no xlogscale \n'

    if ylog:
        output += 'set ylogscale\n'
    else:
        output += '# no ylogscale\n'

    if 'xaxis' in desc:
        if 'label' in desc['xaxis']:
            output += 'set xlabel "' + desc['xaxis']['label'] +'"\n'
    
    if 'yaxis' in desc:
        if 'label' in desc['yaxis']:
            output += 'set ylabel "' + desc['yaxis']['label'] +'"\n'
                    
    if 'legend' in desc:
        output += 'set key top right\n'
        
    num = 0
    output += 'plot '
    for q in desc['data']:
        if len(q.y):
            try:
                xerrors = np.array([xx.error for xx in q.x])
            except AttributeError:
                xerrors = None
                
            try:
                yerrors = np.array([xx.error for xx in q.y])
            except AttributeError:
                yerrors = None
        if 'label' in q.props:
            if xerrors == None and yerrors == None:
                output += ' "-" using 1:2 title "' + q.props['label'] + '",'
            if xerrors == None and yerrors != None:
                output += ' "-" using 1:2:3 w yerrorbars  title "' + q.props['label'] + '",'
            if xerrors != None and yerrors == None:
                output += ' "-" using 1:2:3 w xerrorbars  title "' + q.props['label'] + '",'
            if xerrors != None and yerrors != None:
                output += ' "-" using 1:2:3:4 w xyerrorbars  title "' + q.props['label'] + '",'
        else:
            if xerrors == None and yerrors == None:
                output += ' "-" using 1:2 notitle ,"' 
            if xerrors == None and yerrors != None:
                output += ' "-" using 1:2:3 w yerrorbars  notitle ,' 
            if xerrors != None and yerrors == None:
                output += ' "-" using 1:2:3 w xerrorbars  notitle ,' 
            if xerrors != None and yerrors != None:
                output += ' "-" using 1:2:3:4 w xyerrorbars  notitle ,'
                
        output=output[:-1]
        output+='\n'
    for q in desc['data']:    
            if xerrors == None and yerrors == None:
                output += '# X Y \n'
                for i in range(len(q.x)):
                    output += str(q.x[i]) + '\t' + str(q.y[i]) + '\n'
                output += 'end \n'
            if xerrors == None and yerrors != None:
                output += '# X Y DY \n'
                for i in range(len(q.x)):
                    output += str(q.x[i]) + '\t' + str(q.y[i].mean) + '\t' + str(q.y[i].error) + '\n'
                output += 'end \n'
            if xerrors != None and yerrors == None:
                output += '# X Y DX \n'
                for i in range(len(q.x)):
                    output += str(q.x[i].mean) + '\t' + str(q.y[i]) + '\t' + str(q.x[i].error) + '\n'
                output += 'end \n'
            if xerrors != None and yerrors != None:
                output += '# X Y DXY \n'
                for i in range(len(q.x)):
                    output += str(q.x[i].mean) + '\t' + str(q.y[i].mean) + '\t' + str(q.x[i].error) + '\t' + str(q.y[i].error) + '\n'
                output += 'end \n'
            output += '\n'
            num+=1
                     
    return output


