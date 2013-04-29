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


def converged(h5_outfile, observables, includeLog=False):
  if isinstance(observables, str):
    observables = [observables];

  results = [];
  for observable in observables:
    measurements = pyalps.hdf5.iArchive(h5_outfile).read("/simulation/results/" + observable);

    result = (measurements['mean']['error_convergence'] == 0);

    if not includeLog:
      results.append(result);
    else:
      mean  = measurements['mean']['value'];
      error = measurements['mean']['error'];
      tau   = measurements['tau']['value'];
      count = measurements['count'];
      results.append({'observable': observable, 'converged': result, 'mean': mean, 'error': error, 'tau': tau, 'count': count});

  return results;


def tau(h5_outfile, observables):
  if isinstance(observables, str):
    observables = [observables];

  results = [];
  for observable in observables:
    measurements = pyalps.hdf5.iArchive(h5_outfile).read("/simulation/results/" + observable);
    results.append(measurements['tau']['value']);

  return results;




