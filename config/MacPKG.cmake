## Prepare Distribution.xml file and make_package.sh script used for the Mac Package system

include(${CMAKE_BINARY_DIR}/CPackConfig.cmake)

set(PACK_IDENTIFIER "org.comp-phys.alps.${CPACK_PACKAGE_NAME}")
set(PACK_LOGO "${PROJECT_SOURCE_DIR}/config/alps_logo.png")

## helper macro creating javascript code for dependency
macro(build_selected_code RET COMP)
  string(TOUPPER ${COMP} _COMP_UNAME)
  
  set(_ret "my.choice.selected")
  foreach(_DEP_COMP ${CPACK_COMPONENT_${_COMP_UNAME}_DEPENDS})
    set(_ret "${_ret} &amp;&amp; choices['${_DEP_COMP}Choice'].selected")
  endforeach()
  foreach(_DEP_COMP ${COMP_${_COMP_UNAME}_DEPENDANTS})
    set(_ret "${_ret} || choices['${_DEP_COMP}Choice'].selected")
  endforeach()
  set(${RET} " selected=\"${_ret}\"")
endmacro(build_selected_code)


set(GROUPS_ALL)
foreach(COMP ${CPACK_COMPONENTS_ALL})
    # message(STATUS "Found component ${COMP}")
    
    string(TOUPPER ${COMP} _CPACK_COMP_UNAME)
    
    if(CPACK_COMPONENT_${_CPACK_COMP_UNAME}_GROUP)
        # message(STATUS "-- In group: ${CPACK_COMPONENT_${_CPACK_COMP_UNAME}_GROUP}")
        list(APPEND GROUPS_ALL ${CPACK_COMPONENT_${_CPACK_COMP_UNAME}_GROUP})
        
        set(COMP_GROUP ${CPACK_COMPONENT_${_CPACK_COMP_UNAME}_GROUP})
        string(TOUPPER ${COMP_GROUP} _COMP_GROUP_UNAME)
        list(APPEND GROUP_${_COMP_GROUP_UNAME}_COMPONENTS ${COMP})
    endif()
    
    ## Setup dependants lookup
    foreach(_DEP_COMP ${CPACK_COMPONENT_${_CPACK_COMP_UNAME}_DEPENDS})
      string(TOUPPER ${_DEP_COMP} _DEP_COMP_UNAME)
      list(APPEND COMP_${_DEP_COMP_UNAME}_DEPENDANTS ${COMP})
      # message(STATUS "Dependants for ${_DEP_COMP_UNAME}: ${COMP_${_DEP_COMP_UNAME}_DEPENDANTS}")
    endforeach()
    
endforeach()

if(GROUPS_ALL)
  list(REMOVE_DUPLICATES GROUPS_ALL)
  # message(STATUS "All groups: ${GROUPS_ALL}")
endif(GROUPS_ALL)

set(CHOICES_OUTLINE "")
set(CHOICES "")
set(PKGREFS "")
foreach(COMP_GROUP ${GROUPS_ALL})
    string(TOUPPER ${COMP_GROUP} _COMP_GROUP_UNAME)
    set(CHOICES "${CHOICES}\n<choice id=\"${COMP_GROUP}Choice\" title=\"${CPACK_COMPONENT_GROUP_${_COMP_GROUP_UNAME}_DISPLAY_NAME}\" description=\"${CPACK_COMPONENT_GROUP_${_COMP_GROUP_UNAME}_DESCRIPTION}\" auth=\"Admin\" />")
    
    set(CHOICES_OUTLINE "${CHOICES_OUTLINE}\n<line choice=\"${COMP_GROUP}Choice\">")
    foreach(COMP ${GROUP_${_COMP_GROUP_UNAME}_COMPONENTS})
        string(TOUPPER ${COMP} _COMP_UNAME)
        set(_comp_id "${PACK_IDENTIFIER}.${COMP}")
        set(_comp_hidden "")
        if(CPACK_COMPONENT_${_COMP_UNAME}_HIDDEN)
            set(_comp_hidden " visible=\"false\"")
        endif()
        build_selected_code(_select_code ${COMP})
        
        set(CHOICES_OUTLINE "${CHOICES_OUTLINE}\n  <line choice=\"${COMP}Choice\" />")
        set(CHOICES "${CHOICES}\n  <choice id=\"${COMP}Choice\" title=\"${CPACK_COMPONENT_${_COMP_UNAME}_DISPLAY_NAME}\" description=\"${CPACK_COMPONENT_${_COMP_UNAME}_DESCRIPTION}\" ${_select_code} ${_comp_hidden} auth=\"Admin\">")
        set(CHOICES "${CHOICES}\n    <pkg-ref id=\"${_comp_id}\" />")
        set(CHOICES "${CHOICES}\n  </choice>")
        set(PKGREFS "${PKGREFS}\n  <pkg-ref id=\"${_comp_id}\" version=\"${CPACK_PACKAGE_VERSION}\" auth=\"Admin\" onConclusion=\"None\">${CPACK_PACKAGE_FILE_NAME}-${COMP}.pkg</pkg-ref>")
    endforeach()
    
    set(CHOICES_OUTLINE "${CHOICES_OUTLINE}\n</line>")
endforeach()


foreach(COMP ${CPACK_COMPONENTS_ALL})
    string(TOUPPER ${COMP} _COMP_UNAME)
    if(NOT CPACK_COMPONENT_${_COMP_UNAME}_GROUP)
        set(_comp_id "${PACK_IDENTIFIER}.${COMP}")
        set(_comp_hidden "")
        if(CPACK_COMPONENT_${_COMP_UNAME}_HIDDEN)
            set(_comp_hidden " visible=\"false\"")
        endif()
        build_selected_code(_select_code ${COMP})
        
        set(CHOICES_OUTLINE "${CHOICES_OUTLINE}\n<line choice=\"${COMP}Choice\" />")
        set(CHOICES "${CHOICES}\n<choice id=\"${COMP}Choice\" title=\"${CPACK_COMPONENT_${_COMP_UNAME}_DISPLAY_NAME}\" description=\"${CPACK_COMPONENT_${_COMP_UNAME}_DESCRIPTION}\" ${_select_code} ${_comp_hidden} auth=\"Admin\">")
        set(CHOICES "${CHOICES}\n  <pkg-ref id=\"${_comp_id}\" />")
        set(CHOICES "${CHOICES}\n</choice>")
        set(PKGREFS "${PKGREFS}\n<pkg-ref id=\"${_comp_id}\" version=\"${CPACK_PACKAGE_VERSION}\" auth=\"Admin\" onConclusion=\"None\">${CPACK_PACKAGE_FILE_NAME}-${COMP}.pkg</pkg-ref>")
    endif(NOT CPACK_COMPONENT_${_COMP_UNAME}_GROUP)
endforeach()

set(PACK_CHOICE_OPTIONS "allow")
if("${CPACK_COMPONENTS_ALL}" STREQUAL "")
  set(PACK_CHOICE_OPTIONS "never")
  
  set(COMP "default")
  set(_comp_id "${PACK_IDENTIFIER}")
  
  set(CHOICES_OUTLINE "${CHOICES_OUTLINE}\n<line choice=\"${COMP}Choice\" />")
  set(CHOICES "${CHOICES}\n<choice id=\"${COMP}Choice\" visible=\"false\" auth=\"Admin\">")
  set(CHOICES "${CHOICES}\n  <pkg-ref id=\"${_comp_id}\" />")
  set(CHOICES "${CHOICES}\n</choice>")
  set(PKGREFS "${PKGREFS}\n<pkg-ref id=\"${_comp_id}\" version=\"${CPACK_PACKAGE_VERSION}\" auth=\"Admin\" onConclusion=\"None\">${CPACK_PACKAGE_FILE_NAME}-default.pkg</pkg-ref>")
endif("${CPACK_COMPONENTS_ALL}" STREQUAL "")

string (REPLACE ";" " " CPACK_COMPONENTS_ALL_STR "${CPACK_COMPONENTS_ALL}")
set(PACK_LOGO "${PROJECT_SOURCE_DIR}/config/alps_logo.png")
configure_file(${PROJECT_SOURCE_DIR}/config/make_package.sh.in ${CMAKE_BINARY_DIR}/make_package.sh)
configure_file(${PROJECT_SOURCE_DIR}/config/Distribution.xml.in ${CMAKE_BINARY_DIR}/Distribution.xml)


add_custom_target(mac_package COMMAND "/bin/bash" "${CMAKE_BINARY_DIR}/make_package.sh")
