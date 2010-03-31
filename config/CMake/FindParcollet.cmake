#  Copyright Olivier Parcollet and Matthias Troyer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
#  Python settings : 
#
#  This module checks that : 
#  - the python interpreter is working and version >= 2.5
#  - it has modules : distutils, numpy, tables
# 
#  This module defines the variables
#  - PYTHON_INTERPRETER : name of the python interpreter
#  - PYTHON_INCLUDE_DIR : include for compilation
#  - PYTHON_INCLUDE_NUMPY : include for compilation with numpy
#  - PYTHON_LIBRARY : link flags 
#  - PYTHON_SITE_PKG : path to the standard packages of the python interpreter
#  - PYTHON_EXTRA_LIBS :  libraries which must be linked in when embedding
#  - PYTHON_LINK_FOR_SHARED :  linking flags needed when building a shared lib for external modules

set(PYTHON_INTERPRETER "python" CACHE STRING "Command for the Python interpreter" )
set(PYTHON_MINIMAL_VERSION 2.5)

#
# The function EXEC_PYTHON_SCRIPT executes the_script in  python interpreter
# and set the variable of output_var_name in the calling scope
#
FUNCTION ( EXEC_PYTHON_SCRIPT the_script output_var_name)
  EXECUTE_PROCESS(COMMAND ${PYTHON_INTERPRETER} -c "${the_script}" 
    OUTPUT_VARIABLE res RESULT_VARIABLE returncode OUTPUT_STRIP_TRAILING_WHITESPACE)
  IF (NOT returncode EQUAL 0)
    MESSAGE(FATAL_ERROR "The script : ${the_script} \n did not run properly in the Python interpreter. Check your python installation.") 
  ENDIF (NOT returncode EQUAL 0)
  SET( ${output_var_name} ${res} PARENT_SCOPE)
ENDFUNCTION (EXEC_PYTHON_SCRIPT)

#
# Check the interpreter and its version
#
EXEC_PYTHON_SCRIPT ("import sys, string; print sys.version.split()[0]" PYTHON_VERSION)
STRING(COMPARE GREATER ${PYTHON_MINIMAL_VERSION} ${PYTHON_VERSION} PYTHON_VERSION_NOT_OK)
IF (PYTHON_VERSION_NOT_OK)
MESSAGE(FATAL_ERROR "Python intepreter version is ${PYTHON_VERSION} will it should be >= ${PYTHON_MINIMAL_VERSION}")
ENDIF (PYTHON_VERSION_NOT_OK)
EXEC_PYTHON_SCRIPT ("import distutils " nulle) # check that distutils is there...
EXEC_PYTHON_SCRIPT ("import numpy" nulle) # check that numpy is there...
EXEC_PYTHON_SCRIPT ("import tables" nulle) # check that tables is there...
MESSAGE(STATUS "Python interpreter ok : version ${PYTHON_VERSION}" )

#
# Check for Python include path
#
EXEC_PYTHON_SCRIPT ("import distutils ; from distutils.sysconfig import * ; print distutils.sysconfig.get_python_inc()"  PYTHON_INCLUDE_DIR )
message(status "PYTHON_INCLUDE_DIR =  ${PYTHON_INCLUDE_DIR}" )
mark_as_advanced(PYTHON_INCLUDE_DIR)

#
# include files for numpy
#
EXEC_PYTHON_SCRIPT ("import numpy;print numpy.get_include()" PYTHON_INCLUDE_NUMPY)
MESSAGE(STATUS "PYTHON_INCLUDE_NUMPY = ${PYTHON_INCLUDE_NUMPY}" )
mark_as_advanced(PYTHON_INCLUDE_NUMPY)

#
# Check for Python library path
#
#EXEC_PYTHON_SCRIPT ("import string; from distutils.sysconfig import * ;print string.join(get_config_vars('VERSION'))"  PYTHON_VERSION_MAJOR_MINOR)         
EXEC_PYTHON_SCRIPT ("import string; from distutils.sysconfig import * ;print '-L%s/config -lpython%s'%(get_python_lib(0,1),string.join(get_config_vars('VERSION')))" PYTHON_LIBRARY)
MESSAGE(STATUS "PYTHON_LIBRARY = ${PYTHON_LIBRARY}" )
mark_as_advanced(PYTHON_LIBRARY)

#
# Check for site packages
#
EXEC_PYTHON_SCRIPT ("from distutils.sysconfig import * ;print get_python_lib(0,0)"
		    PYTHON_SITE_PKG)
MESSAGE(STATUS "Python site packages found it in: ${PYTHON_SITE_PKG}" )
mark_as_advanced(PYTHON_SITE_PKG)

#
# libraries which must be linked in when embedding
#
EXEC_PYTHON_SCRIPT ("from distutils.sysconfig import * ;print (str(get_config_var('LOCALMODLIBS')) + ' ' + str(get_config_var('LIBS'))).strip()"
		    PYTHON_EXTRA_LIBS)
MESSAGE(STATUS "PYTHON_EXTRA_LIBS =${PYTHON_EXTRA_LIBS}" )
mark_as_advanced(PYTHON_EXTRA_LIBS)

#
# linking flags needed when embedding (building a shared lib)
# To BE RETESTED
#
EXEC_PYTHON_SCRIPT ("from distutils.sysconfig import *;print get_config_var('LINKFORSHARED')"
		    PYTHON_LINK_FOR_SHARED)
MESSAGE(STATUS "PYTHON_LINK_FOR_SHARED =  ${PYTHON_LINK_FOR_SHARED}" )
mark_as_advanced(PYTHON_LINK_FOR_SHARED)

# Correction on Mac
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
SET (PYTHON_LINK_FOR_SHARED -u _PyMac_Error -framework Python)
SET (PYTHON_LINK_MODULE -bundle -undefined dynamic_lookup)
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
SET (PYTHON_LINK_MODULE -shared)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

#
# This function writes down a script to compile f2py modules
# indeed, one needs to use the f2py of the correct numpy module.
#
FUNCTION( WriteScriptToBuildF2pyModule filename fcompiler_desc modulename module_pyf_name filelist )
  # Copy all the files
  EXECUTE_PROCESS(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${module_pyf_name} ${CMAKE_CURRENT_BINARY_DIR} )
  FOREACH( f ${filelist})
    EXECUTE_PROCESS(COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR} )
  ENDFOREACH(f)
  # write the script that will build the f2py extension
  SET(filename ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
  FILE(WRITE ${filename} "import sys\n")
  FILE(APPEND ${filename} "from numpy.f2py import main\n")
  FILE(APPEND ${filename} "sys.argv = [''] +'-c --fcompiler=${fcompiler_desc} -m ${modulename} ${modulename}.pyf ${filelist} -llapack'.split()\n")
  FILE(APPEND ${filename} "main()\n")
ENDFUNCTION(WriteScriptToBuildF2pyModule)



