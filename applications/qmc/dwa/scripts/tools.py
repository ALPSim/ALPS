import numpy;
import scipy;
import pyalps;

def thermalized(h5_outfile, observables, tolerance=0.01, includeLog=False):
  if isinstance(observables, str):
    observables = [observables];

  results = [];
  for observable in observables:
    timeseries = pyalps.hdf5.iArchive(h5_outfile).read("/simulation/results/" + observable)['timeseries']['data']; 
    mean = timeseries.mean();

    index = scipy.linspace(0, timeseries.size-1, timeseries.size);
    timeseries = scipy.polyval(scipy.polyfit(index, timeseries, 1), index);  # timeseries get fitted

    percentage_increment = (timeseries[-1] - timeseries[0])/mean;
    result = abs(percentage_increment) < tolerance;

    if not includeLog:
      results.append(result);
    else:
      results.append({'observable': observable, 'percentage_increment' : percentage_increment, 'thermalized': result})

  return results;

