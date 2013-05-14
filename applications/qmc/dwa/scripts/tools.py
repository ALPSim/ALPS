import numpy;
import scipy;
import pyalps;

def thermalized(h5_outfile, observables, tolerance=0.01, simplified=False, includeLog=False):
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

  if includeLog or not simplified:
    return results;
  else:
    return reduce(lambda x,y: x*y, results);



def converged(h5_outfile, observables, simplified=False, includeLog=False):
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

  if includeLog or not simplified:
    return results;
  else:
    return reduce(lambda x,y: x*y, results);


def tau(h5_outfile, observables):
  if isinstance(observables, str):
    observables = [observables];

  results = [];
  for observable in observables:
    measurements = pyalps.hdf5.iArchive(h5_outfile).read("/simulation/results/" + observable);
    results.append(measurements['tau']['value']);

  return results;


def write_status(filename, status):
  ar = pyalps.hdf5.h5ar(filename, 'w');
  if ar.is_group('/status'):  
    ar['/status/' + str(len(ar.list_children('/status')))] = status; 
  else:
    ar['/status/0'] = status;

def status(filename, includeLog=False):
  ar = pyalps.hdf5.h5ar(filename);
  if not ar.is_group('/status'):
    return None;
  else:
    if not includeLog:
      return ar['/status/' + str(len(ar.list_children('/status'))-1)];
    else:
      return ar['/status'] 


def recursiveRun(cmd, cmd_lang='command_line', follow_up_script=None, n=None, break_if=None, break_elseif=None, write_status=None, loc=None):
  ### 
  ### Either recursively run cmd for n times, or until the break_if condition holds true.
  ###
  ###  Note:
  ###    1. cmd              : command to be recursively run (quoted as a python str)
  ###    2. cmd_lang         : language of cmd, either "command_line" (default), or "python".
  ###    3. n                : number of recursions 
  ###    4. break_if         : condition to break recursion loop (quoted as a python str, interpreted as python command)
  ###    5. break_elseif     : further condition to break recursion loop (""")     
  ###    5. follow_up_script : script to be run after command (""")
  ###    6. loc              : Python dict of local variables 
  ###

  if loc != None:
    locals().update(loc);

  if cmd_lang == 'command_line':
    pyalps.executeCommand(cmd.split());
  elif cmd_lang == 'python':
    eval(cmd);
  else:
    print "Error: The options for cmd_lang are 1) 'command_line' (default), or 2) 'python'."
    return;

  if follow_up_script != None:  
    eval(follow_up_script);

  if write_status != None:
    eval(write_status);

  if n != None:             # if n exists
    if isinstance(n, int):  # if n is a python integer
      if n <= 1:
        return;
      else:
        return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, n=n-1, write_status=write_status, loc=locals()); 

  elif break_if != None:    # otherwise, if break_if exists
    if eval(break_if):
      return;
    else:
      if break_elseif != None:   # otherotherwise, if break_elseif exists
        if eval(break_elseif):
          return;
        else:
          return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, break_if=break_if, break_elseif=break_elseif, write_status=write_status, loc=locals());
      else:
        return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, break_if=break_if, write_status=write_status, loc=locals());


  else:                     # otherwise, recursiveRun only runs once
    return;

