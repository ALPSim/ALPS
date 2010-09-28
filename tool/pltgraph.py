#  Copyright Bela Bauer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import matplotlib
matplotlib.use('macosx')

import matplotlib.pyplot as plt
from xml.etree import ElementTree

import sys

def parse(fn):
    root = ElementTree.parse(fn).getroot()
    vertices = {}
    
    for vertex in root.findall('VERTEX'):
        vid = vertex.get('id')
        vpos = tuple([float(x) for x in vertex.find('COORDINATE').text.split()])
        vertices[vid] = vpos
    
    edges = {}
    for edge in root.findall('EDGE'):
        edges[edge.get('id')] = dict([(k,edge.get(k)) for k in ['source', 'target', 'type', 'vector']])
    
    return (vertices,edges)

def showgraph(graph):
    vertices = graph[0]
    edges = graph[1]
    
    x = [v[0] for v in vertices.values()]
    y = [v[1] for v in vertices.values()]
    plt.scatter(x, y)
    
    for edge in edges.values():
        print edge
        s = edge['source']
        t = edge['target']
        
        p0 = vertices[s]
        p1 = vertices[t]
        
        plt.plot([p0[0], p1[0]], [p0[1], p1[1]])

if __name__ == '__main__':
    graph = parse(sys.argv[1])
    showgraph(graph)
    plt.show()
