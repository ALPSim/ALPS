# How to stop a simulation      {#stop_simulation}

[TOC]

## Number of sweeps
The default behavior is to run until the desired number of sweeps is reached:

    SWEEPS = X

However sometimes the application might be killed by the scheduling system or other stopping criteria might be desired. See the following sections for more options.

## Time limit
The MPS applications have a simple command line parameter to specify the maximum run time. If the limit is reached, the application will stop and safely write down the current state of the simulation such that it can be restarted. The limit is checked at the end of each iteration.

The time limit is specified with the option ```-T X``` with *X* the number of seconds it is allowed to run.

When using Python to run the applications, the time limit parameter is an argument of the ```runApplication()``` function. For example:

~~~~~~~~~~~~~{.py}
res = pyalps.runApplication('mps_optim',input_file,T=300)
~~~~~~~~~~~~~

## Energy threshold
Often you don't know how many sweeps are needed for convergence. Here the best is to submit the submit the calculation with a large numebr of sweeps (now playing the role of maximum number of sweeps) and enable the *relative energy threshold*. At each sweep the application computes the relative energy change, if this is below the desired threshold, the calculation terminates and the current MPS wave function is writte to disk as a checkpoint file.

To enable this feature, just use for example:

    SWEEPS = 20
    rel_en_thresh    = 1e-8
    rel_en_thresh_at = "half"

The value of ```rel_en_thresh_at``` specifies at which site to check for the energy threshold. Possible values are "half" (in the middle of the system) and "end" (at the right end of the system).

*Note:* the energy convergence is not a universal indicator for convergence. It happens that the calculation get stuck in unphysical local minima which might have only 1e-4 relative energy difference to the actual ground state. Usually, looking at local observables it is pretty straight forward to see that these don't fulfill the expected symmetry.
