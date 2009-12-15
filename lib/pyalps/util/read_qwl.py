from xml.etree import ElementTree
import numpy as np

from dataset import DataSet

def read_qwl(filename):
    root = ElementTree.parse(filename).getroot()

    show_legend = (root.find('legend').attrib['show'] == 'true')

    xaxis = {'label': root.find('xaxis').attrib['label']}
    yaxis = {'label': root.find('yaxis').attrib['label']}

    x = []
    y = []
    for point in root.find('set').getchildren():
        x.append(point.find('x').text)
        y.append(point.find('y').text)

    data = DataSet()
    data.x = np.array(x)
    data.y = np.array(y)
    
    plotd = {}
    plotd['xaxis'] = xaxis
    plotd['yaxis'] = yaxis
    plotd['data'] = data
    if show_legend:
        plotd['legend'] = {}
    
    return plotd
