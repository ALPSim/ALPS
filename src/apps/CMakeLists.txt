#  Copyright Matthias Troyer 2009 - 2010.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

######################################################################
# applications
######################################################################

# The components built with Eigen (ALPS_EIGEN_COMPONENT_NAMES, which have their
# own exports) that a target links, directly or through other targets.
function(alps_linked_eigen_components target out)
  get_property(eigen_components GLOBAL PROPERTY ALPS_EIGEN_COMPONENT_NAMES)
  set(found "")
  set(seen "")
  set(queue ${target})
  while(queue)
    list(POP_FRONT queue current)
    # Generator expressions and external libraries are not ALPS components.
    if(NOT TARGET "${current}" OR current IN_LIST seen)
      continue()
    endif()
    list(APPEND seen ${current})
    get_target_property(aliased ${current} ALIASED_TARGET)
    if(aliased)
      list(APPEND queue ${aliased})
      continue()
    endif()
    get_target_property(imported ${current} IMPORTED)
    if(imported)
      continue()
    endif()
    get_target_property(name ${current} EXPORT_NAME)
    if(name AND name IN_LIST eigen_components)
      list(APPEND found ${name})
    endif()
    foreach(property IN ITEMS LINK_LIBRARIES INTERFACE_LINK_LIBRARIES)
      get_target_property(libraries ${current} ${property})
      if(libraries)
        list(APPEND queue ${libraries})
      endif()
    endforeach()
  endwhile()
  list(REMOVE_DUPLICATES found)
  set(${out} ${found} PARENT_SCOPE)
endfunction()

function(alps_install_applications)
  install(
    TARGETS ${ARGN}
    EXPORT ALPSApplicationTargets
    RUNTIME_DEPENDENCY_SET alps_runtime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT applications)
  # ALPSConfig.cmake provides these components' libraries with the applications,
  # so that consumers bundling the programs (pyalps) need not know about them.
  foreach(application IN LISTS ARGN)
    alps_linked_eigen_components(${application} components)
    set_property(GLOBAL APPEND PROPERTY ALPS_APPLICATION_EIGEN_COMPONENTS ${components})
  endforeach()
  list(TRANSFORM ARGN PREPEND "ALPS::")
  set_property(TARGET alps APPEND PROPERTY ALPS_APPLICATION_TARGETS ${ARGN})
endfunction()

add_subdirectory(mc)
add_subdirectory(qmc)
add_subdirectory(diag)
add_subdirectory(dmrg)
add_subdirectory(dmft)
add_subdirectory(maxent)

install(
  EXPORT ALPSApplicationTargets
  NAMESPACE ALPS::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  COMPONENT applications)
