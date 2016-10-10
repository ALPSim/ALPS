# Documentation      {#mainpage}

## Using the code
Example simulations are available in the `tutorials/mps-XX-YY` folders. In order to run them you need a complete installation of ALPS, i.e. including the Python tools. The python scripts must be executed with the `alpspython` program, such that the ALPS libraries are imported correclty.

For the first tutorial:

    $ cd tutorials/mps-01-optim
    $ alpspython spin_one_half.py

### Parameters
In the [Parameters section](@ref parameters) one can find a detailed description of the simulation parameters for all the ALPS %MPS tools.

### Guides
* [Simulation files](@ref simulation_files)
* [How to stop a simulation](@ref stop_simulation)

## File list
A short description of the source files is available in the [File list section](files.html).

## Installation


### Dependencies

 * CMake >= 2.8.x
 * Boost C++ libraries >= 1.52.0
 * HDF5 >= 1.8.2
 * BLAS and LAPACK libraries
 * Python 2.6.x or 2.7.x with Numpy, SciPy, Matplotlib.

Example package names for the above dependencies are available on the [ALPS website][alps-web].


### Build

The CMake configuration system will find automatically most 
dependencies. In case of configuration problems, please refer to the 
[ALPS website][alps-web] and the [ALPS-Users forum][alps-mail].

For a Unix or Darwin system, the following command are enough to build 
the package. For Windows it is advised to use the CMake GUI utility or 
to download a binary installer (>=2.2.0b1) from the [ALPS website][alps-web].

#### Configuration

    $ cd build
    $ cmake ../alps

#### Compile and install

    $ make # -jN for compilation in parallel using N processes
    $ make install


## License
See alps/LICENSE.txt and alps/LICENSE-applications.txt. Details on the [ALPS license website](https://alps.comp-phys.org/mediawiki/index.php/Licensing).

Use of `mps_optim`, `mps_tevol`, `mps_meas` or `mps_overlap` requires 
citation of the ALPS %MPS paper \cite dolfi2014. Use of any ALPS program requires 
citation of the ALPS paper \cite alps-paper.

## References
See [Bibliography](citelist.html)  page.


[alps-web]: http://alps.comp-phys.org
[alps-mail]: comp-phys-alps-users@lists.comp-phys.org
