*******************************
API Reference Manual for Python
*******************************

.. toctree::
   :maxdepth: 4
   
   .

Run Applications
================
.. autosummary::
  pyalps.writeParameterFile
  pyalps.writeInputFiles

  pyalps.runApplication
  pyalps.runDMFT
  pyalps.evaluateLoop
  pyalps.evaluateSpinMC
  pyalps.evaluateQWL
  pyalps.evaluateFulldiagVersusT
  pyalps.evaluateFulldiagVersusH

.. autofunction:: pyalps.writeParameterFile
.. autofunction:: pyalps.writeInputFiles

.. autofunction:: pyalps.runApplication
.. autofunction:: pyalps.runDMFT
.. autofunction:: pyalps.evaluateLoop
.. autofunction:: pyalps.evaluateSpinMC
.. autofunction:: pyalps.collectXY
.. autofunction:: pyalps.evaluateQWL
.. autofunction:: pyalps.evaluateFulldiagVersusT
.. autofunction:: pyalps.evaluateFulldiagVersusH


Loading data
============
.. autosummary::
  pyalps.getResultFiles

  pyalps.loadMeasurements
  pyalps.loadBinningAnalysis
  pyalps.loadEigenstateMeasurements
  pyalps.loadSpectra
  pyalps.loadDMFTIterations
  pyalps.loadProperties
  pyalps.loadObservableList

.. autofunction:: pyalps.getResultFiles

.. autofunction:: pyalps.loadMeasurements
.. autofunction:: pyalps.loadBinningAnalysis
.. autofunction:: pyalps.loadEigenstateMeasurements
.. autofunction:: pyalps.loadSpectra
.. autofunction:: pyalps.loadDMFTIterations
.. autofunction:: pyalps.loadProperties
.. autofunction:: pyalps.loadObservableList


Evaluation
==========
.. autosummary::
  pyalps.DataSet
  
  pyalps.collectXY
  pyalps.groupSets
  pyalps.select
  pyalps.select_by_property
  pyalps.mergeDataSets
  pyalps.mergeMeasurements
  pyalps.select
  
  pyalps.fit_wrapper.Parameter
  pyalps.fit_wrapper.fit
  
  pyalps.SetLabels
  pyalps.CycleColors
  pyalps.CycleMarkers
  pyalps.collectXY

.. autoclass:: pyalps.DataSet
  :members:


Tools
-----
.. autofunction:: pyalps.collectXY
.. autofunction:: pyalps.groupSets
.. autofunction:: pyalps.select
.. autofunction:: pyalps.select_by_property
.. autofunction:: pyalps.mergeDataSets
.. autofunction:: pyalps.mergeMeasurements

Fit wrapper
-----------
.. autofunction:: pyalps.fit_wrapper.Parameter
.. autofunction:: pyalps.fit_wrapper.fit

Plot
----
.. autofunction:: pyalps.SetLabels
.. autofunction:: pyalps.CycleColors
.. autofunction:: pyalps.CycleMarkers

