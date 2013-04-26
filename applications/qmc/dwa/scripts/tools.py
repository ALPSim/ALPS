def thermalized(parmfiles, observables):
  if isinstance(parmfiles, str):
    parmfiles = [parmfiles];
  if isinstance(observables, str):
    observables = [observables];

  return [parmfiles, observables];
