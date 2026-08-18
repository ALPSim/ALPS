##############################################################################
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
##############################################################################

import config
import wx
from wx.lib.hyperlink import HyperLinkCtrl

alpsDescription = """The ALPS project (Algorithms and Libraries for Physics Simulations) is an open source effort aiming at providing high-end simulation codes for strongly correlated quantum mechanical systems as well as C++ libraries for simplifying the development of such code. ALPS strives to increase software reuse in the physics community."""

alpsLicense = """MIT License

Copyright 2003-2025 ALPS Collaboration

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."""

class AboutThisSoftware(wx.Frame):
    def __init__(self, parent, name, version = config.version(), copyright = config.copyright()):
        wx.Frame.__init__(self, parent, -1, size=(480,400))
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        title = wx.StaticText(self, label=name + " version " + version)
        title.SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD, False, 'Verdana'))
        sizer.Add(title, 0, wx.EXPAND|wx.ALL, 5)

        copy = wx.StaticText(self, -1, copyright)
        sizer.Add(copy, 0, wx.EXPAND|wx.ALL, 5)

        wiki = HyperLinkCtrl(self, -1, "ALPS Wiki", URL="http://alps.comp-phys.org")
        sizer.Add(wiki, 0, wx.ALL, 5)

        desc = wx.StaticText(self, -1, alpsDescription, size=(400,90), style=wx.TE_MULTILINE)
        sizer.Add(desc, 0, wx.EXPAND|wx.ALL, 5)

        lic = wx.TextCtrl(self, -1, alpsLicense, style=wx.TE_MULTILINE|wx.TE_READONLY,
                          size=(400,100))
        sizer.Add(lic, 1, wx.EXPAND|wx.ALL, 5)

        btn = wx.Button(self, -1, "Close")
        btn.SetDefault()
        sizer.Add(btn, 0, wx.ALIGN_CENTER|wx.ALL, 5)
        self.Bind(wx.EVT_BUTTON, self.OnButton)

        sizer.Layout()

    def OnButton(self, event):
        self.Destroy()

if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    frame = AboutThisSoftware(None, 'My Program')
    frame.Show()
    app.MainLoop()
