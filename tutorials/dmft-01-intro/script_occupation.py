# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2012 by Jakub Imriska <jimriska@phys.ethz.ch> 
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


import pyalps

print "DESCRIPTION: This script prints out the density for each flavor, for single site problems; the results are loaded from the h5-files in the current directory (recursively)."
print

result_files = pyalps.getResultFiles()

#print "OCCUPATION (in the last iteration) taken directly from the result files :"
#occ = pyalps.loadMeasurements(result_files,['n'])
#for o in occ:
#  if len(o)>0:
#    print "  from file : ",o[0].props['filename']
#    for f in range(0,len(o[0].x)):
#      print '    flavor ',f,' : ',o[0].y[f]
#print
      
print "DENSITY (in the last iteration) obtained from the Green's function in imaginary time as -G(beta^-):"
ll=pyalps.load.Hdf5Loader()
for a in result_files:
  res_file = [a]
  obs = pyalps.loadObservableList(res_file)
  if obs[0].count('G_tau')>0:
    flavors = int(pyalps.loadProperties(res_file)[0]["FLAVORS"])
    #print "  from file : ",a
    listobs=[]
    for f in range(0,flavors):
      listobs.append(str(f))  # previous format: "Green_"+str(f)
    data_G_tau = ll.ReadMeasurementFromFile(res_file, respath='/simulation/results/G_tau', measurements=listobs, verbose=True)  # [result_files(here only 1)][measurements(here: flavors)]
    for f in range(0,flavors):
      print '    flavor ',f,' : ',-data_G_tau[0][f].y[len(data_G_tau[0][f].x)-1]

