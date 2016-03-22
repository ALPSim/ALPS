##############################################################################
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
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
##############################################################################

import config, license
import os, random, subprocess, sys
from math import cos, sin, pi
from xml.dom import minidom

import wx
import vtk
from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

def prog():
    return "ALPS Lattice Preview"

def systemLibraryPath():
    if os.path.exists(config.prefix() + '/lib/xml/lattices.xml'):
        return config.prefix() + '/lib/xml/lattices.xml'
    else:
        return config.builddir() + '/../lib/xml/lattices.xml'

def lattice2xml():
    if os.path.exists(config.builddir() + '/lattice2xml'):
        return config.builddir() + '/lattice2xml'
    else:
        return config.prefix() + '/bin/lattice2xml'

class LatticeLibrary:
    def __init__(self, libraryPath = ''):
        self.graphs = []
        self.params = {}
        if libraryPath:
            self.parseXML(libraryPath)
    def parseXML(self, libraryPath):
        self.graphs = []
        self.params = {}
        lattices = []
        latticeparams = {}
        top = minidom.parse(libraryPath)
        for n0 in top.childNodes:
            if n0.localName == 'LATTICES':
                for n1 in n0.childNodes:
                    if n1.localName == 'LATTICE':
                        name = n1.getAttribute("name")
                        p = {}
                        for n2 in n1.childNodes:
                            if n2.localName == 'PARAMETER':
                                p[n2.getAttribute("name")] = n2.getAttribute("default")
                        lattices.append([name, p])
                        latticeparams[name] = p
        for n0 in top.childNodes:
            if n0.localName == 'LATTICES':
                for n1 in n0.childNodes:
                    if n1.localName == 'LATTICEGRAPH':
                        name = n1.getAttribute("name")
                        p = {}
                        for n2 in n1.childNodes:
                            if n2.localName == 'FINITELATTICE':
                                for n3 in n2.childNodes:
                                    if n3.localName == 'LATTICE':
                                        ref = n3.getAttribute("ref")
                                        if ref and latticeparams.has_key(ref):
                                            for k in latticeparams[ref].keys():
                                                p[k] = latticeparams[ref][k]
                                    elif n3.localName == 'PARAMETER':
                                        p[n3.getAttribute("name")] = n3.getAttribute("default")
                                    elif n3.localName == 'EXTENT':
                                        size = n3.getAttribute("size")
                                        if not p.has_key(size) and not size.isdigit():
                                            p[size] = ""
                        self.graphs.append(name)
                        self.params[name] = p
                    elif n1.localName == 'GRAPH':
                        name = n1.getAttribute("name")
                        self.graphs.append(name)
                        self.params[name] = {}
        top.unlink()

class LatticeData:
    # def __init__(self, paramfile = ''):
    def __init__(self, param = {}):
        self.lattice2xml = lattice2xml()
        self.clear()
        if param:
            self.parseXML(param)

    def clear(self):
        self.vertices = []
        self.edges = []
        self.hasCoordinate = False
        self.vertexTypes = []
        self.maxVertexType = 0
        self.edgeTypes = []
        self.maxEdgeType = 0

    def parseXML(self, param):
        self.clear()
        p = subprocess.Popen(self.lattice2xml, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, close_fds=True)
        (cout, cerr) = p.communicate(param)
        if p.returncode == 0:
            xmltree = minidom.parseString(cout)
            self.parseTree(xmltree)
            xmltree.unlink()
            if not self.hasCoordinate:
                self.assignCoordinates()
            self.updateVertexTypes()
            self.updateEdgeTypes()
        return (p.returncode, cerr)

    def parseTree(self, node):
        if node.localName == 'EDGE':
            source = node.getAttribute("source")
            target = node.getAttribute("target")
            tp = node.getAttribute("type")
            out = node.getAttribute("outside")
            if tp == "" or tp == None:
                tp = "0"
            if out == "" or out == None:
                out = "0"
            self.edges.append([int(source), int(target), int(tp), int(out)])
        elif node.localName == 'VERTEX':
            tp = node.getAttribute("type")
            out = node.getAttribute("outside")
            if tp == ""    or tp == None:
                tp = "0"
            if out == "" or out == None:
                out = "0"
            coords = [0, 0, 0]
            for child in node.childNodes:
                if child.localName == 'COORDINATE':
                    for grandchild in child.childNodes:
                        if grandchild.nodeType == grandchild.TEXT_NODE:
                            vec = grandchild.wholeText.split(" ")
                            if len(vec) >= 1:
                                self.hasCoordinate = True
                                coords[0] = float(vec[0])
                            if len(vec) >= 2:
                                coords[1] = float(vec[1])
                            if len(vec) >= 3:
                                coords[2] = float(vec[2])
            self.vertices.append([coords, int(tp), int(out)])
        else:
            for child in node.childNodes:
                self.parseTree(child)

    def assignCoordinates(self):
        n = len(self.vertices)
        for s in range(0, n):
            self.vertices[s][0][0] = cos(2*pi*s/n)
            self.vertices[s][0][1] = sin(2*pi*s/n)

    def updateVertexTypes(self):
        types = {}
        for (coord, type, out) in self.vertices:
            types[type] = 1
        for k in types.keys():
            self.vertexTypes.append(int(k))
        self.vertexTypes.sort()
        self.maxVertexType = self.vertexTypes[len(self.vertexTypes)-1]

    def updateEdgeTypes(self):
        types = {}
        for (source, target, type, out) in self.edges:
            types[type] = 1
        for k in types.keys():
            self.edgeTypes.append(int(k))
        self.edgeTypes.sort()
        self.maxEdgeType = self.edgeTypes[len(self.edgeTypes)-1]

class LatticeParameterWindow(wx.Frame):
    def __init__(self, parent, libraryPath = ''):
        wx.Frame.__init__(self, parent, -1, "Lattice Preview", size=wx.Size(600,600))
        self.libraryPath = libraryPath
        self.useSystem = True
        if libraryPath:
            self.useSystem = False
        self.graphName = ''
        self.parameters = {}
        self.library = LatticeLibrary()

        self.sizer = wx.BoxSizer(wx.VERTICAL)

        # menu bar
        menuBar = wx.MenuBar()
        menuPreview = wx.Menu()
        menuPreviewAbout = menuPreview.Append(-1, "About...", "About")
        menuPreview.AppendSeparator()
        menuPreviewNew = menuPreview.Append(-1, "&New...\tCtrl+N", "New Preview")
        menuPreviewClose = menuPreview.Append(-1, "&Close\tCtrl+W", "Close Window")
        self.Bind(wx.EVT_MENU, self.ShowAbout, menuPreviewAbout)
        self.Bind(wx.EVT_MENU, self.OnNew, menuPreviewNew)
        self.Bind(wx.EVT_MENU, self.OnClose, menuPreviewClose)
        menuBar.Append(menuPreview, "Preview")
        self.SetMenuBar(menuBar)

        #
        # Lattice Library
        #

        self.sizer.Add(wx.StaticText(self, -1, "Lattice Library:"), 0, wx.EXPAND|wx.ALL, 5)

        self.rb_system = wx.RadioButton(self, -1, "System Lattice Library", style=wx.RB_GROUP)
        self.rb_user = wx.RadioButton(self, -1, "User Lattice XML")
        self.tx_user = wx.TextCtrl(self, -1, self.libraryPath, size=wx.Size(200,10),
                                   style=wx.TE_READONLY)
        btn_choose = wx.Button(self, -1, "Choose")
        self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio, self.rb_system)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnRadio, self.rb_user)
        self.Bind(wx.EVT_BUTTON, self.OnChoose, btn_choose)
        self.sizer.Add(self.rb_system, 0, wx.EXPAND|wx.ALL, 5)
        self.sizer.Add(self.rb_user, 0, wx.EXPAND|wx.ALL, 5)
        choose_sizer = wx.BoxSizer(wx.HORIZONTAL)
        choose_sizer.Add((10,10), 0, wx.EXPAND|wx.ALL, 5)
        choose_sizer.Add(self.tx_user, 1, wx.EXPAND|wx.ALL, 5)
        choose_sizer.Add(btn_choose, 0, wx.EXPAND|wx.ALL, 5)
        self.sizer.Add(choose_sizer, 0, wx.EXPAND|wx.ALL, 5)
        self.sizer.Add(wx.StaticLine(self), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)

        #
        # Lattice/Graph
        #

        self.sizer.Add(wx.StaticText(self, -1, "Lattice/Graph:"), 0, wx.EXPAND|wx.ALL, 5)

        graphs = []
        self.graphBox = wx.Choice(self, -1, choices=graphs)
        self.Bind(wx.EVT_CHOICE, self.OnSelectGraph, self.graphBox)
        self.sizer.Add(self.graphBox, 0, wx.EXPAND|wx.ALL, 15)
        self.sizer.Add(wx.StaticLine(self), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 5)

        self.sizer.Add(wx.StaticText(self, -1, "Parameters:"), 0, wx.EXPAND|wx.ALL, 5)
        self.fgs = wx.FlexGridSizer(0, 2, 5, 5)
        self.sizer.Add(self.fgs, 1, wx.EXPAND|wx.ALL, 15)

        btn_cancel = wx.Button(self, -1, "Cancel")
        self.btn_preview = wx.Button(self, -1, "Preview")
        self.btn_preview.SetDefault()
        self.Bind(wx.EVT_BUTTON, self.OnCancel, btn_cancel)
        self.Bind(wx.EVT_BUTTON, self.OnPreview, self.btn_preview)
        btns = wx.BoxSizer(wx.HORIZONTAL)
        btns.Add((1, 1), 1, wx.RIGHT, 5)
        btns.Add(btn_cancel, 0, wx.RIGHT, 5)
        btns.Add(self.btn_preview, 0, wx.RIGHT, 5)
        self.sizer.Add(btns, 0, wx.EXPAND|wx.ALL, 5)

        self.SetSizer(self.sizer)
        self.UpdateLibrary()

    def ShowAbout(self, evt):
        frame = license.AboutThisSoftware(self, prog())
        frame.CentreOnParent(wx.BOTH)
        frame.Show()

    def OnNew(self, evt):
        frame = LatticeParameterWindow(None)
        frame.Show()

    def OnClose(self, evt):
        self.Destroy()

    def UpdateLibrary(self):
        if self.useSystem:
            self.rb_system.SetValue(True)
            self.tx_user.Enable(False)
            self.library.parseXML(systemLibraryPath())
        else:
            self.rb_user.SetValue(True)
            self.tx_user.Enable(True)
            if self.libraryPath:
                self.library.parseXML(self.libraryPath)
        self.UpdateGraph()

    def UpdateGraph(self):
        graphs = []
        if len(self.library.graphs) > 0:
            graphs = ["Please choose a lattice/graph"]
            graphs += self.library.graphs
        self.graphBox.SetItems(graphs)
        self.graphBox.SetSelection(0)
        self.graphName = ''
        self.UpdateParameters()
        self.btn_preview.Enable(False)

    def OnRadio(self, evt):
        if evt.GetEventObject().GetLabel() == 'System Lattice Library':
            if self.useSystem == False:
                self.useSystem = True
                self.UpdateLibrary()
        else:
            if self.useSystem == True:
                if self.tx_user.GetValue() == '':
                    self.OnChoose(True)
                else:
                    self.useSystem = False
                    self.UpdateLibrary()

    def OnChoose(self, event):
        wildCard = "XML file (*.xml)|*.xml|All files (*.*)|*.*"
        path = os.getcwd()
        if self.libraryPath:
            path = os.path.dirname(self.libraryPath)
        dialog = wx.FileDialog(self.Parent, "Choose an XML file", path, "", wildCard, wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            if self.useSystem or not self.libraryPath == dialog.GetPath():
                self.useSystem = False
                self.libraryPath = dialog.GetPath()
                self.UpdateLibrary()
            self.tx_user.SetValue(self.libraryPath)
        else:
            self.UpdateLibrary()
        dialog.Destroy()

    def OnSelectGraph(self, event):
        if self.graphBox.GetSelection() > 0 :
            graph = self.library.graphs[self.graphBox.GetSelection() - 1]
            self.btn_preview.Enable(True)
        else:
            graph = ''
            self.btn_preview.Enable(False)
        if not graph == self.graphName:
            self.graphName = graph
            self.UpdateParameters()

    def UpdateParameters(self):
        self.fgs.Clear(True)
        if self.graphName:
            self.parameters = self.library.params[self.graphName]
            if len(self.parameters):
                for k in self.parameters.keys():
                    name = wx.StaticText(self, -1, k + " :  ")
                    value = wx.TextCtrl(self, -1, self.parameters[k], size=wx.Size(200,-1), name=k)
                    self.Bind(wx.EVT_TEXT, self.OnParameterInput, value)
                    self.fgs.Add(name, 0, wx.ALIGN_LEFT)
                    self.fgs.Add(value, 1, wx.ALIGN_LEFT)
            else:
                name = wx.StaticText(self, -1, "None")
                self.fgs.Add(name, 0, wx.ALIGN_LEFT)
                name.Enable(False)
        self.sizer.Layout()

    def OnParameterInput(self, evt):
        key = evt.GetEventObject().GetName()
        if self.parameters.has_key(key):
            self.parameters[key] = evt.GetEventObject().GetValue()

    def OnPreview(self, evt):
        params = ""
        if self.useSystem:
            params += "LATTICE_LIBRARY = \"" + systemLibraryPath() + "\"; "
        else:
            params += "LATTICE_LIBRARY = \"" + self.libraryPath + "\"; "
        params += "LATTICE = \"" + self.graphName + "\"; "
        params += "UNROLL_BOUNDARY = 1; "
        for k in self.parameters.keys():
            params += k + " = \"" + self.parameters[k] + "\"; "
        lattice = LatticeData()
        (ret, cerr) = lattice.parseXML(params)
        if ret == 0:
            frame = PreviewLatticeWindow(lattice, size=self.GetSize(), pos=self.GetPosition())
            frame.Show()
            self.Destroy()
        else:
            dialog = wx.MessageDialog(self, cerr, "Error (code = " + str(ret) + ")",
                                      wx.OK|wx.ICON_ERROR)
            dialog.ShowModal()

    def OnCancel(self, evt):
        self.Destroy()

class PreviewLatticeWindow(wx.Frame):
    def __init__(self, lattice, size, pos):
        wx.Frame.__init__(self, None, -1, "Lattice Preview", size=size, pos=pos)
        self.lattice = lattice

        # menu bar
        menuBar = wx.MenuBar()
        menuFile = wx.Menu()
        menuFileAbout = menuFile.Append(-1, "About...", "About")
        menuFile.AppendSeparator()
        menuFileNew = menuFile.Append(-1, "&New...\tCtrl+N", "New Preview")
        menuFileClose = menuFile.Append(-1, "&Close\tCtrl+W", "Close Window")
        self.Bind(wx.EVT_MENU, self.ShowAbout, menuFileAbout)
        self.Bind(wx.EVT_MENU, self.OnNew, menuFileNew)
        self.Bind(wx.EVT_MENU, self.OnClose, menuFileClose)
        menuBar.Append(menuFile, "Preview")
        self.menuView = wx.Menu()
        if len(lattice.vertexTypes) > 0:
            self.menuView.Append(-1, "Vertex Type").Enable(False)
            for t in lattice.vertexTypes:
                m = self.menuView.AppendCheckItem(-1, "  " + str(t))
                m.Check(True)
                self.Bind(wx.EVT_MENU, self.OnViewVertex, m)
        if len(lattice.vertexTypes) > 0 and len(lattice.edgeTypes) > 0:
            self.menuView.AppendSeparator()
        if len(lattice.edgeTypes) > 0:
            self.menuView.Append(-1, "Edge Type").Enable(False)
            for t in lattice.edgeTypes:
                m = self.menuView.AppendCheckItem(-1, "  " + str(t))
                m.Check(True)
                self.Bind(wx.EVT_MENU, self.OnViewEdge, m)
        if len(lattice.vertexTypes) > 0 or len(lattice.edgeTypes) > 0:
            self.menuView.AppendSeparator()
            m = self.menuView.Append(-1, "Show All")
            self.Bind(wx.EVT_MENU, self.OnShowAll, m)
            m = self.menuView.Append(-1, "Hide All")
            self.Bind(wx.EVT_MENU, self.OnHideAll, m)
        menuBar.Append(self.menuView, "View")
        self.SetMenuBar(menuBar)

        main = wx.BoxSizer(wx.VERTICAL)
        self.SetBackgroundColour('#eeffff')

        self.vtkwidget = wxVTKRenderWindowInteractor(self, -1)
        main.Add(self.vtkwidget, 1, wx.EXPAND)
        self.SetSizer(main)
        self.Layout()

        self.vtkwidget.Enable(1)
        self.vtkwidget.AddObserver("ExitEvent", lambda o, e, f=self: f.Close())

        self.vertexActors = []
        self.vertexView = []
        self.edgeActors   = []
        self.edgeView = []
        self.InitColors()

        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.1, 0.2, 0.4)
        self.vtkwidget.GetRenderWindow().AddRenderer(self.renderer)
        vtk.vtkOutputWindow().PromptUserOff()

        self.ShowPreview()

    def InitColors(self):
        random.seed()
        self.VertexColors = [
            [0.8,0.3,0.3], #dark red
            [0.4,0.6,0.5], #dark green 2
            [0.7,0.5,0.2]  #orange 2
            ]
        for i in range(3, self.lattice.maxVertexType+1):
            self.VertexColors.append([random.random(),random.random(),random.random()])
        self.EdgeColors = [
            [0.4,0.6,0.5], #dark green 2
            [0.8,0.3,0.5], #dark red
            [1.0,0.4,0.2]  #orange 2
            ]
        for i in range(3, self.lattice.maxEdgeType+1):
            self.EdgeColors.append([random.random(),random.random(),random.random()])

    def ShowPreview(self):
        self.vertexActors = []
        self.vertexView = []
        self.edgeActors = []
        self.edgeView = []

        for v in self.lattice.vertices:
            self.vertexActors.append(vtk.vtkActor())
            self.vertexView.append(True)
        for e in self.lattice.edges:
            self.edgeActors.append(vtk.vtkActor())
            self.edgeView.append(True)

        for i in range(len(self.lattice.vertices)):
            (coord, tp, out) = self.lattice.vertices[i]
            color = self.VertexColors[tp]
            self.drawVertex(self.vertexActors[i], coord, color)
            self.renderer.AddActor(self.vertexActors[i])
        for i in range(len(self.lattice.edges)):
            (source, target, tp, out) = self.lattice.edges[i]
            color = self.EdgeColors[tp]
            self.drawEdge(self.edgeActors[i], self.lattice.vertices[source-1][0],
                          self.lattice.vertices[target-1][0], color)
            self.renderer.AddActor(self.edgeActors[i])
        self.setOpacity()

    def OnViewVertex(self, evt):
        type = int(evt.GetEventObject().GetLabel(evt.GetId()))
        self.vertexView[type] = evt.IsChecked()
        self.setOpacity()

    def OnViewEdge(self, evt):
        type = int(evt.GetEventObject().GetLabel(evt.GetId()))
        self.edgeView[type] = evt.IsChecked()
        self.setOpacity()

    def OnShowAll(self, evt):
        for m in evt.GetEventObject().GetMenuItems():
            if m.IsCheckable():
                m.Check(True)
        for i in range(len(self.vertexView)):
            self.vertexView[i] = True
        for i in range(len(self.edgeView)):
            self.edgeView[i] = True
        self.setOpacity()

    def OnHideAll(self, evt):
        for m in evt.GetEventObject().GetMenuItems():
            if m.IsCheckable():
                m.Check(False)
        for i in range(len(self.vertexView)):
            self.vertexView[i] = False
        for i in range(len(self.edgeView)):
            self.edgeView[i] = False
        self.setOpacity()

    def ShowAbout(self, evt):
        frame = license.AboutThisSoftware(self, prog(), copyright())
        frame.CentreOnParent(wx.BOTH)
        frame.Show()

    def OnNew(self, evt):
        frame = LatticeParameterWindow(None)
        frame.Show()

    def OnClose(self, evt):
        self.Destroy()

    def drawVertex(self, actor, coord, color):
        v = vtk.vtkSphereSource()
        v.SetPhiResolution(20)
        v.SetThetaResolution(20)
        v.SetCenter(coord[0], coord[1], coord[2])
        v.SetRadius(0.05)
        actor.GetProperty().SetColor(color)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(v.GetOutputPort())
        actor.SetMapper(mapper)

    def drawEdge(self, actor, source, target, color):
        e = vtk.vtkPolyData()
        p = vtk.vtkPoints()
        p.InsertPoint(0, source[0], source[1], source[2])
        p.InsertPoint(1, target[0], target[1], target[2])
        e.SetPoints(p)
        c = vtk.vtkCellArray()
        c.InsertNextCell(2)
        c.InsertCellPoint(0)
        c.InsertCellPoint(1)
        e.SetLines(c)
        edge = vtk.vtkTubeFilter()
        edge.SetRadius(0.01)
        edge.SetNumberOfSides(10)
        edge.SetInput(e)
        actor.GetProperty().SetColor(color)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(edge.GetOutputPort())
        actor.SetMapper(mapper)

    def setOpacity(self):
        for i in range(len(self.lattice.vertices)):
            (coord, tp, out) = self.lattice.vertices[i]
            if self.vertexView[tp]:
                if out == 1:
                    self.vertexActors[i].GetProperty().SetOpacity(0.3)
                else:
                    self.vertexActors[i].GetProperty().SetOpacity(1.0)
            else:
                self.vertexActors[i].GetProperty().SetOpacity(0.0)
        for i in range(len(self.lattice.edges)):
            (source, target, tp, out) = self.lattice.edges[i]
            if self.edgeView[tp]:
                if out == 1:
                    self.edgeActors[i].GetProperty().SetOpacity(0.3)
                else:
                    self.edgeActors[i].GetProperty().SetOpacity(1.0)
            else:
                self.edgeActors[i].GetProperty().SetOpacity(0.0)
        self.vtkwidget.Render()

## main routine
if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    path = ''
    if len(sys.argv) > 1:
        path = os.path.abspath(sys.argv[1])
    frame = LatticeParameterWindow(None, path)
    frame.Show()
    app.MainLoop()
