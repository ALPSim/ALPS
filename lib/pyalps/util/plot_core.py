import matplotlib.pyplot as plt
import numpy as np

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
    
    def __init__(self):
        self.icolor = 0
        self.output = ''
    
    def draw_lines(self):
        self.lines = []
        self.icolor = 0
        
        xlog = False
        ylog = False
        if 'xaxis' in self.plt and 'logarithmic' in self.plt['xaxis']:
            xlog = self.plt['xaxis']['logarithmic']
        if 'yaxis' in self.plt and 'logarithmic' in self.plt['yaxis']:
            ylog = self.plt['yaxis']['logarithmic']
        
        print self.plt    
        for q in self.plt['data']:
            line_props = self.colors[self.icolor]
            if 'line' in q.props:
                line_props += q.props['line']
            
            try:
                xmeans = np.array([xx.mean for xx in q.x])
                xerrors = np.array([xx.error for xx in q.x])
            except AttributeError:
                xmeans = q.x
                xerrors = None
            
            try:
                ymeans = np.array([xx.mean for xx in q.y])
                yerrors = np.array([xx.error for xx in q.y])
            except AttributeError:
                ymeans = q.y
                yerrors = None
            
            print xmeans, ymeans, xerrors, yerrors
            self.lines.append(plt.errorbar(xmeans,ymeans,yerr=xerrors,xerr=yerrors,fmt=line_props))
            
            if xlog:
                plt.xscale('log')
            if ylog:
                plt.yscale('log')
            
            if 'label' in q.props and q.props['label'] != 'none':
                self.lines[-1][0].set_label(q.props['label'])
            elif 'filename' in q.props:
                self.lines[-1][0].set_label(q.props['filename'])
            
            self.icolor = (self.icolor+1)%len(self.colors)
            
    def __call__(self, desc):
        self.plt = desc
        
        self.draw_lines()
        
        if 'xaxis' in self.plt:
            if 'label' in self.plt['xaxis']:
                plt.xlabel(self.plt['xaxis']['label'])
            if 'min' in self.plt['xaxis'] and 'max' in self.plt['xaxis']:
                plt.xlim(self.plt['xaxis']['min'],self.plt['xaxis']['max'])
                
        if 'yaxis' in self.plt:
            if 'label' in self.plt['yaxis']:
                plt.ylabel(self.plt['yaxis']['label'])
            if 'min' in self.plt['yaxis'] and 'max' in self.plt['yaxis']:
                plt.ylim(self.plt['yaxis']['min'],self.plt['yaxis']['max'])
        
        if 'legend' in self.plt:
            if 'location' in self.plt['legend']:
                plt.legend(loc=self.plt['legend']['location'])
            else:
                plt.legend()
        
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
        output += '#\n@g0 on\n@with g0\n'
        output += '@ frame linewidth 2.0'
        output += '@page background fill off'

        xrange = [0,1]
        yrange = [0,1]
        if 'xaxis' in desc and 'min' in desc['xaxis'] and 'max' in desc['xaxis']: 
            xrange = [ desc['xaxis']['min'],data['xaxis']['max']]
        if 'yaxis' in desc and 'min' in desc['yaxis'] and 'max' in desc['yaxis']:
            yrange = [ desc['yaxis']['min'],data['yaxis']['max']]

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
            output += '@    xaxes scale Logarithmic\n'
        else:
            output += '@    xaxes scale Normal\n'

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
        for q in desc['data']:
            output += '@target G0.S'+str(num)+'\n'
            output += '@    s'+str(num)+' symbol ' + str(num+1) +'\n'
            output += '@    s'+str(num)+' symbol size 0.500000\n'
            output += '@    s'+str(num)+' line type 1\n'
            if 'label' in q.props and q.props['label'] != 'none':
                output += '@    s'+str(num)+' legend "' + q.props['label'] + '"\n'
            elif 'filename' in q.props:
                output += '@    s'+str(num)+' legend "' + q.props['label'] + '"\n'
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
