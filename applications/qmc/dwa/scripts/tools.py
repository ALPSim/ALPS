import numpy;
import scipy;
import pyalps;
import inspect;

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

def extract_worldlines(infile, outfile):
  wl = pyalps.dwa.worldlines()
  wl.load(infile);
  wl.save(outfile);
  return;

def recursiveRun(cmd, cmd_lang='command_line', follow_up_script=None, n=None, break_if=None, break_elseif=None, write_status=None, loc=None, batch_submit=False, batch_cmd_prefix=None, batch_run_script='run.script', batch_next_run_script=None, batch_run_now=False, batch_noRun=False):
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
    return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, n=n, break_if=break_if, break_elseif=break_elseif, write_status=write_status, batch_submit=batch_submit, batch_cmd_prefix=batch_cmd_prefix, batch_run_script=batch_run_script, batch_next_run_script=batch_next_run_script, batch_run_now=batch_run_now, batch_noRun=batch_noRun);

  # preparing batch run script 
  if batch_submit:
    if not batch_run_now:
      batch_cmd = '';
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
        return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, n=n-1, write_status=write_status, batch_submit=batch_submit, batch_cmd_prefix=batch_cmd_prefix, batch_run_script=batch_run_script, batch_next_run_script=batch_next_run_script, batch_run_now=False); 

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
          return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, break_if=break_if, break_elseif=break_elseif, write_status=write_status, batch_submit=batch_submit, batch_cmd_prefix=batch_cmd_prefix, batch_run_script=batch_run_script, batch_next_run_script=batch_next_run_script, batch_run_now=False);
      else:
        return recursiveRun(cmd, cmd_lang=cmd_lang, follow_up_script=follow_up_script, break_if=break_if, write_status=write_status, batch_submit=batch_submit, batch_cmd_prefix=batch_cmd_prefix, batch_run_script=batch_run_script, batch_next_run_script=batch_next_run_script, batch_run_now=False);

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


