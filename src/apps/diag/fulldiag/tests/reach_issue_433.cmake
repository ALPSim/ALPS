# Issue #433: fulldiag_evaluate given the task's .out.h5 file (written next
# to the .out.xml) must either evaluate it like the .out.xml, or stop with a
# usage error (exit status 1) and write no plots. It must not die in the
# catch-all exception handler (exit status 251).
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}/xml" "${test_dir}/h5")
set(ENV{ALPS_XML_PATH} "${xml_dir}")

function(run dir status_var)
  execute_process(COMMAND ${ARGN}
    WORKING_DIRECTORY "${test_dir}/${dir}" OUTPUT_VARIABLE out ERROR_VARIABLE err
    RESULT_VARIABLE status)
  set(${status_var} "${status}" PARENT_SCOPE)
  set(output "${out}${err}" PARENT_SCOPE)
endfunction()

foreach(dir xml h5)
  configure_file("${input}" "${test_dir}/${dir}/fd" COPYONLY)
  run(${dir} s "${parameter2xml}" fd)
  if(NOT s EQUAL 0)
    message(FATAL_ERROR "parameter2xml failed (${s}):\n${output}")
  endif()
  run(${dir} s "${fulldiag}" fd.in.xml)
  if(NOT s EQUAL 0)
    message(FATAL_ERROR "fulldiag failed (${s}):\n${output}")
  endif()
endforeach()

set(plots energy entropy free_energy magnetization specific_heat uniform_susceptibility)

# The .out.xml file evaluates.
run(xml s "${fulldiag_evaluate}" --T_MIN 0.5 --T_MAX 1.0 --DELTA_T 0.5 fd.task1.out.xml)
if(NOT s EQUAL 0)
  message(FATAL_ERROR "fulldiag_evaluate on fd.task1.out.xml failed (${s}):\n${output}")
endif()
foreach(p ${plots})
  if(NOT EXISTS "${test_dir}/xml/fd.task1.plot.${p}.xml")
    message(FATAL_ERROR "no ${p} plot from fd.task1.out.xml:\n${output}")
  endif()
endforeach()

# The .out.h5 file.
run(h5 s "${fulldiag_evaluate}" --T_MIN 0.5 --T_MAX 1.0 --DELTA_T 0.5 fd.task1.out.h5)
if(s EQUAL 0)
  foreach(p ${plots})
    set(f "fd.task1.plot.${p}.xml")
    if(NOT EXISTS "${test_dir}/h5/${f}")
      message(FATAL_ERROR "no ${p} plot from fd.task1.out.h5:\n${output}")
    endif()
    file(READ "${test_dir}/xml/${f}" a)
    file(READ "${test_dir}/h5/${f}" b)
    if(NOT a STREQUAL b)
      message(FATAL_ERROR "${f} from fd.task1.out.h5 differs from the one from fd.task1.out.xml")
    endif()
  endforeach()
elseif(s EQUAL 1)
  file(GLOB written "${test_dir}/h5/*.plot.*")
  if(written)
    message(FATAL_ERROR "fd.task1.out.h5 was refused but plots were written: ${written}")
  endif()
else()
  message(FATAL_ERROR "fulldiag_evaluate on fd.task1.out.h5 neither evaluated it "
    "nor refused it as a usage error (exit status ${s}, expected 0 or 1):\n${output}")
endif()
