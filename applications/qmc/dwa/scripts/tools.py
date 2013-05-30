import os;
import numpy;
import scipy;
import matplotlib;
import matplotlib.pyplot;
import pyalps;
from pyalps.dwa.worldlines_c import worldlines;

def format_string(string, loc):
  result = []; 
  string = string.replace("(", "( ");
  string = string.replace(",", " , ");
  string = string.replace(")", " )");
  string = string.replace("=", " = ");
  for item in string.split():
    try:
      item = "'" + ("{" + item + "}").format(**loc) + "'";
    except Exception:
      pass;
    result.append(item);
  string = " ".join(result);
  string = string.replace("( ", "(");
  string = string.replace(" , ", ", ");
  string = string.replace(" )", ")");  
  return string;

def str_quote(item):
  return ('"' + str(item) + '"');

def cyclic_shift(array, by):
  # array: numpy.array
  if not isinstance(array,numpy.ndarray):
    raise Exception("array must be of type <numpy.ndarray>");

  for i in range(array.ndim):
    array = numpy.roll(array, by[i], axis=i);

  return array;

def addToObservable(h5_outfile, RealObservable=None, RealVectorObservable=None, measurement=None):
  if measurement == None:
    return;
  if RealObservable == None and RealVectorObservable == None:
    return;
  if RealObservable != None and RealVectorObservable != None:
    return;

  ar = pyalps.ngs.archive(h5_outfile, 'w');
  measurements = pyalps.ngs.observables();
  measurements.load(ar, '/simulation/results');
  
  ar.set_context('/simulation/results');
  if RealObservable != None:
    if not ar.is_group(RealObservable):
      measurements.createRealObservable(RealObservable);
    measurements[RealObservable] << measurement;
  elif RealVectorObservable != None: 
    if not ar.is_group(RealVectorObservable):
      measurements.createRealVectorObservable(RealVectorObservable);
    measurements[RealVectorObservable] << measurement;
      
  measurements.save(ar);

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

def extract_worldlines(infile, outfile=None):
  wl = worldlines()
  wl.load(infile);

  if outfile != None:
    wl.save(outfile);
    return;
  else:
    return wl;

def show_worldlines(wl=None, fig_outfile=None, reshape=None, at=None, Nmax=20):
  if wl == None:
    return;

  [wl_siteindicator, wl_time, wl_state] = [wl.worldlines_siteindicator(), wl.worldlines_time(), wl.worldlines_state()];

  return wl_siteindicator;

  wl_coordinates = [];
  for idx1 in range(wl.num_sites()):
    for idx2 in range(1,wl.num_kinks(idx1)):
      wl_coordinates.append([idx1,wl_time[idx1][idx2]]);
  [wl_coordinates_site, wl_coordinates_time] = numpy.array(wl_coordinates).transpose();

  wl_state_segments = [];
  for idx1 in range(wl.num_sites()):
    for idx2 in range(wl.num_kinks(idx1)-1):
      wl_state_segments.append([[idx1,wl_time[idx1][idx2]], [idx1,wl_time[idx1][idx2+1]], wl_state[idx1][idx2]])
    wl_state_segments.append([[idx1,wl_time[idx1][wl.num_kinks(idx1)-1]], [idx1,1.0], wl_state[idx1][wl.num_kinks(idx1)-1]]);
  
  wl_n_state_segments = [];
  for n in range(Nmax+1):
    wl_n_state_segments.append([wl_state_segment[0:2] for wl_state_segment in wl_state_segments if wl_state_segment[2] == n]);  

  wl_vertex_segments = [];
  for idx1 in range(wl.num_sites()):
    for idx2 in range(1,wl.num_kinks(idx1)):
      if   (wl_siteindicator[idx1][idx2] - idx1 == 1):
        wl_vertex_segments.append([[idx1,wl_time[idx1][idx2]],[idx1+0.5,wl_time[idx1][idx2]]]);
      elif (wl_siteindicator[idx1][idx2] - idx1 == -1):
        wl_vertex_segments.append([[idx1-0.5,wl_time[idx1][idx2]],[idx1,wl_time[idx1][idx2]]]);
      elif (idx1 == 0):
        wl_vertex_segments.append([[idx1-0.5,wl_time[idx1][idx2]],[idx1,wl_time[idx1][idx2]]]);
      else:
        wl_vertex_segments.append([[idx1-0.5,wl_time[idx1][idx2]],[idx1+0.5,wl_time[idx1][idx2]]]);

  if fig_outfile != None:
    matplotlib.use('Agg');

  matplotlib.pyplot.figure(frameon=False);
  matplotlib.pyplot.xticks([]);
  matplotlib.pyplot.yticks([0,1]);
  

'''
def garbage():  
  ### Plot worldlines
  
  if (len(sys.argv) == 2):
    matplotlib.pyplot.xlim(-0.5,wl.num_sites()+0.5);
  else:
    xlim = eval(sys.argv[2]);
    matplotlib.pyplot.xlim(xlim[0]-0.5,xlim[1]-0.5);
  
  matplotlib.pyplot.ylim(-0.05,1.05);
  
  # scatter layout
  #matplotlib.pyplot.scatter(wl_coordinates_site, wl_coordinates_time);
  
  for wl_n_state_segment in wl_n_state_segments[0]:
    [segment_site, segment_time] = numpy.array(wl_n_state_segment).transpose();  
    matplotlib.pyplot.plot(segment_site, segment_time, '--k', linewidth=2);
  
  for wl_n_state_segment in wl_n_state_segments[1]:
    [segment_site, segment_time] = numpy.array(wl_n_state_segment).transpose();
    matplotlib.pyplot.plot(segment_site, segment_time, '-k', linewidth=2);
  
  wl_n_state_segment_line_spacing = 0.2;  
  
  for n in range(2,21):
    for wl_n_state_segment in wl_n_state_segments[n]:
      [segment_site, segment_time] = numpy.array(wl_n_state_segment).transpose();
      for m in range(n):
        matplotlib.pyplot.plot(segment_site - (m - (n-1)/2.)*wl_n_state_segment_line_spacing, segment_time, '-k', linewidth=2);
  
  for wl_vertex_segment in wl_vertex_segments:
    [segment_site, segment_time] = numpy.array(wl_vertex_segment).transpose();
    matplotlib.pyplot.plot(segment_site, segment_time, '-k', linewidth=2);
  
  matplotlib.pyplot.show();
  
  ### OUTFILE
  #outfile = infile[0:-3] + ".pdf";
  #matplotlib.pyplot.savefig(outfile, transparent=True);
'''



def recursiveRun(cmd, cmd_lang='command_line', follow_up_script=None, n=None, break_if=None, break_elseif=None, write_status=None, loc=None, loc0=None, batch_submit=False, batch_cmd_prefix=None, batch_run_directory=None, batch_run_script='run.script', batch_next_run_script=None, batch_run_now=False, batch_noRun=False):
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

  # set absolute path for current path
  if batch_run_directory == None:
    batch_run_directory = os.getcwd();

  # Format string 
  if loc != None:
    locals().update(loc);
    cmd = format_string(cmd, loc);
    if follow_up_script != None:
      follow_up_script = format_string(follow_up_script, loc);
    if break_if != None:
      break_if = format_string(break_if, loc);
    if break_elseif != None:
      break_elseif = format_string(break_elseif, loc);
    if write_status != None:
      write_status = format_string(write_status, loc);
    batch_run_script = format_string(batch_run_script, loc);
    if batch_next_run_script != None:
      batch_next_run_script = format_string(batch_next_run_script, loc);
    return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, n=n, break_if=break_if, break_elseif=break_elseif, write_status=write_status, loc0=loc, batch_submit=batch_submit, batch_cmd_prefix=batch_cmd_prefix, batch_run_directory=batch_run_directory, batch_run_script=batch_run_script, batch_next_run_script=batch_next_run_script, batch_run_now=batch_run_now, batch_noRun=batch_noRun);

  if loc0 != None:
    locals().update(loc0);

  # preparing batch run script 
  if batch_submit:
    if not batch_run_now:
      batch_cmd = '';
      batch_cmd += 'cd ' + str(batch_run_directory) + '\n\n';
      batch_cmd += 'python <<@@\n';
      batch_cmd +=   'import pyalps;\n';
      batch_cmd +=   'import pyalps.dwa\n\n';
      batch_cmd +=   'pyalps.dwa.recursiveRun(' + str_quote(cmd);
      if cmd_lang != None: 
        batch_cmd += ', cmd_lang = ' + str_quote(cmd_lang);
      if follow_up_script != None:
        batch_cmd += ', \n\tfollow_up_script = ' + str_quote(follow_up_script);
      if n != None:
        batch_cmd += ', \n\tn = ' + str(n);
      if break_if != None:
        batch_cmd += ', \n\tbreak_if = ' + str_quote(break_if);
      if break_elseif != None:
        batch_cmd += ', \n\tbreak_elseif = ' + str_quote(break_elseif);
      if write_status != None:
        batch_cmd += ', \n\twrite_status = ' + str_quote(write_status);
      batch_cmd += ', \n\tbatch_submit = ' + str(batch_submit);
      if batch_cmd_prefix != None:
        batch_cmd += ', \n\tbatch_cmd_prefix = ' + str_quote(batch_cmd_prefix);
      batch_cmd += ', \n\tbatch_run_directory = ' + str_quote(batch_run_directory);
      batch_cmd += ', \n\tbatch_run_script = ' + str_quote(batch_run_script);
      if batch_next_run_script != None:
        batch_cmd += ', \n\tbatch_next_run_script = ' + str_quote(batch_next_run_script);
      batch_cmd += ', \n\tbatch_run_now = True';
      batch_cmd +=   ');\n\n';
      batch_cmd += '@@';

      pyalps.executeCommand(['echo', batch_cmd, '>', batch_run_script]);
      pyalps.executeCommand(['chmod','755',batch_run_script]);

      if batch_noRun:
        return;

      command = [];
      if batch_cmd_prefix != None:
        command += batch_cmd_prefix.split();
      command += ['./' + batch_run_script];
      return pyalps.executeCommand(command);

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
        if batch_next_run_script != None:
          command = [];
          if batch_cmd_prefix != None:
            command += batch_cmd_prefix.split();
          command += ['./' + batch_next_run_script];
          return pyalps.executeCommand(command);
        else:
          return;
      else:
        return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, n=n-1, write_status=write_status, loc0=loc0, batch_submit=batch_submit, batch_cmd_prefix=batch_cmd_prefix, batch_run_directory=batch_run_directory, batch_run_script=batch_run_script, batch_next_run_script=batch_next_run_script, batch_run_now=False); 

  elif break_if != None:    # otherwise, if break_if exists
    if eval(break_if):
      if batch_next_run_script != None:
        command = [];
        if batch_cmd_prefix != None:
          command += batch_cmd_prefix.split();
        command += ['./' + batch_next_run_script];
        return pyalps.executeCommand(command);
      else:
        return;
    else:
      if break_elseif != None:   # otherotherwise, if break_elseif exists
        if eval(break_elseif):
          if batch_next_run_script != None:
            command = [];
            if batch_cmd_prefix != None:
              command += batch_cmd_prefix.split();
            command += ['./' + batch_next_run_script];
            return pyalps.executeCommand(command);
          else:
            return;
        else:
          return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, break_if=break_if, break_elseif=break_elseif, write_status=write_status, loc0=loc0, batch_submit=batch_submit, batch_cmd_prefix=batch_cmd_prefix, batch_run_directory=batch_run_directory, batch_run_script=batch_run_script, batch_next_run_script=batch_next_run_script, batch_run_now=False);
      else:
        return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, break_if=break_if, write_status=write_status, loc0=loc0, batch_submit=batch_submit, batch_cmd_prefix=batch_cmd_prefix, batch_run_script=batch_run_script, batch_run_directory=batch_run_directory, batch_next_run_script=batch_next_run_script, batch_run_now=False);

  else:                     # otherwise, recursiveRun only runs once
    if batch_next_run_script != None:
      command = [];
      if batch_cmd_prefix != None:
        command += batch_cmd_prefix.split();
      command += ['./' + batch_next_run_script];
      return pyalps.executeCommand(command);
    else:
      return;


def startRunScript(batch_run_script, batch_cmd_prefix=None, loc=None):
  # Format string 
  if loc != None:
    locals().update(loc);
    batch_run_script = format_string(batch_run_script, loc);
    return startRunScript(batch_run_script, batch_cmd_prefix=batch_cmd_prefix);

  command = [];
  if batch_cmd_prefix != None:
    command += batch_cmd_prefix.split();
  command += ['./' + batch_run_script];
  
  return pyalps.executeCommand(command);


