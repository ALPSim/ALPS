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

import config
import wx
from wx.lib.hyperlink import HyperLinkCtrl

alpsDescription = """The ALPS project (Algorithms and Libraries for Physics Simulations) is an open source effort aiming at providing high-end simulation codes for strongly correlated quantum mechanical systems as well as C++ libraries for simplifying the development of such code. ALPS strives to increase software reuse in the physics community."""

alpsLicense = """ALPS LIBRARY LICENSE version 1.1
Copyright (C) 2003-2005 Ian McCulloch.  Everyone is permitted to copy and distribute this license document.

This License applies to any software containing a notice placed by the copyright holder saying that it may be distributed under the terms of the ALPS Library License version 1.1.  Such software is herein referred to as the "Library". This license grants permission to use, reproduce, display, distribute, execute and transmit the Library, and to prepare derivative works of the Library, and to permit others to do so for non-commercial academic use, all subject to the following conditions:

1. In any scientific publication based wholly or in part on the Library, the use of the Library must be acknowledged and the publications listed in the accompanying CITATIONS.txt document must be cited.

2. You may copy and distribute verbatim copies of the Library in the form that you received it, as long as all copyright notices and references to this license and warranty disclaimer are kept intact, and all recipients also receive a copy of this license, warranty disclaimer and CITATIONS.txt document.

3. You may modify your copy or copies of the Library, thus forming a work based on the Library, and use, copy or distribute such modified works under the terms of sections 1 and 2 above, provided that you also meet all of these conditions:

a. You must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.

b. All citations listed in the CITATIONS.txt document that refer to sections of the Library that exist in the modified work must be preserved irrespective of the extent of the modification.

c. You must cause any work that you distribute or publish, that in whole or in part contains or is derived from the Library or any part thereof, to be licensed as a whole at no charge to all third parties under terms compatible with this License.

4. This Software, or modifications under section 3 above, may be distributed in object code or executable form, provided that you meet all of these conditions:

a. This complete License, warranty disclaimer and accompanying CITATIONS.txt document is included.

b. The executable program is accompanied with the complete machine-readable source code to the Library as used in the executable, which must be distributed under the terms of sections 2 and 3 above. Alternatively, you may provide instructions for obtaining the source code at no cost (for example, a hyper-text link).

5. A program that contains no derivative of any portion of the Library, but is designed to work with the Library by being compiled or linked with it, is not a derivative work of the Library, and therefore falls outside the scope of this License.

However, linking such a work with the Library creates an executable that is a derivative of the Library (because it contains portions of the Library).  The executable is therefore covered by this License.  Section 4 states terms for distribution of such executables.

6. You must cause executable programs that utilize this Library to print or display, when started in the most basic way, a prominent announcement including a copyright notice and citation requirements as listed in the accompanying CITATIONS.txt document.  If the executable program utilizes the Library in a modified form (under section 3 above), then the announcement must state this. Exception: if the announcement would not normally be visible to the user, or the announcement would interfere with normal operations of the executable application, then the executable program is not required to print an announcement.

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT.  IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."""

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
