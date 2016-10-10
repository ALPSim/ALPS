# Simulation files      {#simulation_files}



## Input files
The MPS applications expect a set of XML files describing all the parameter sets to be computed. To generate the input files you may rely on the ALPS standard tools, i.e. using the ```paramter2xml``` application or generate them using Python scripts.

Assuming a basename *parms*, the list of input files is:

    .
    ├── parms.in.xml        # list of all parameter sets and links to their input files
    ├── parms.task1.in.xml  # parameters for the first simulation
    ├── parms.task2.in.xml  # parameters for the second simulation
    ├── ...
    └── parms.taskN.in.xml  # parameters for the n-th simulation

## Output files
After running a simulation you will find a new set of files being generated. Here follows a description.

    .
    ├── parms.out.xml           # list of all parameter sets and links to their output files
    ├── parms.task1.out.chkp/   # the MPS wavefunction. This can be used to continue the calculation, or to initialize new simulations.
    ├── parms.task1.out.h5      # output in HDF5 format. All iteration and final measurements.
    ├── parms.task1.out.xml     # output in XML format. Observables will be present only if running with the `--write-xml` option
    └── ...

Note that to re-run the simulation you have to delete all output files, in particular the *out.chkp* files.

If you enable the calculation of more eigenstates (with the the *NUMBER_EIGENVALUES* parameter), there will be one set of files for each eigenstate:

    .
    ├── parms.out.xml           # list of all parameter sets and links to their output files
    ├── parms.task1.out.0.chkp/ # the MPS wavefunction of the ground state
    ├── parms.task1.out.0.h5    # iteration measurements for computing the ground state
    ├── parms.task1.out.1.chkp/ # the MPS wave function of the first excited state
    ├── parms.task1.out.1.h5    # iteration measurements for computing the first excited states
    ├── parms.task1.out.h5      # final measurements of all eigenstates in HDF5 format
    ├── parms.task1.out.xml     # final measurements of all eigenstates in XML format
    └── ...

