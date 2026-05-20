Installation
============
ALPS comes as a single package containing both libraries and
applications. Unless stated otherwise, the ALPS Libraries are
distributed under the ALPS library license for libraries and ALPS
application license for applications. 


You can install ALPS 2 by simply downloading and installing a binary
package for your system or build ALPS from source yourself. To install
ALPS you can either download and install the binary package for your
system or you may build ALPS from the source yourself.


How to install ALPS
-------------------
Using the LiveALPS Linux distribution

To make your first steps with ALPS, without installing it to your system, the easiest way to proceed is using the LiveALPS Linux distribution.


Installing the binary releases on MacOSX and Windows
----------------------------------------------------

Installing the ALPS applications and binaries
+++++++++++++++++++++++++++++++++++++++++++++
To install the binary releases for MacOS X and Windows just download the
appropriate installer from the download page and install it on your
system. The installation location is "/opt/alps" on MacOS and
"C:\Program Files\ALPS" on Windows. 

If it wants to install into "C:\Program Files (x86)\ALPS", then you have
chosen the 32 bit version of ALPS for a 64 bit version of Windows. While
this works it is not recommended. 

Note: on Windows make sure to select the option to add ALPS to the PATH
when asked. If you forget to do so, this can still be done
manually. Choose "Control Panel" --> "System" --> "Advanced" -->
"Environment Variables" and add the path to the ALPS programs
(e.g. C:\Program Files\ALPS\bin) to the PATH variable. 


Non-default install locations for the binary packages
+++++++++++++++++++++++++++++++++++++++++++++++++++++
If you install it in a non-default location you will have to 

  * set the environment variable ALPS_XML_DIR to point the the directory
    containing the ALPS XML and XSL files (the lib/xml subdirectory of the
    ALPS installation 

  * On MacOS X set the environment variable DYLD_LIBRARY_PATH to include
    the lib directory of the ALPS installation so that the dynamic
    libraries can be found. 

How to build ALPS
-----------------


In case of problems
-------------------
Should you have any questions or inquiries about the ALPS package,
please don't hesitate to consult our mailing list:
comp-phys-alps-users@lists.comp-phys.org
