# FCIDUMP="file" with several tasks: tasks with their own FCIDUMP each write
# theirs; a second task naming the same FCIDUMP fails the run before any task
# runs, instead of overwriting the first one's file.
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")
set(common [=[
LATTICE="chain lattice"
L=2
MODEL="fermion Hubbard"
Nup_total=1
Ndown_total=1
]=])

function(run_fulldiag name)
  execute_process(COMMAND "${parameter2xml}" ${name}
    WORKING_DIRECTORY "${test_dir}" COMMAND_ERROR_IS_FATAL ANY)
  execute_process(COMMAND "${fulldiag}" ${name}.in.xml
    WORKING_DIRECTORY "${test_dir}" OUTPUT_VARIABLE out ERROR_VARIABLE err
    RESULT_VARIABLE status)
  set(output "${out}${err}" PARENT_SCOPE)
  set(status "${status}" PARENT_SCOPE)
endfunction()

function(require_u file u)
  file(STRINGS "${test_dir}/${file}" line REGEX " 1 +1 +1 +1$")
  if(NOT line MATCHES "^ *${u}\\.0+e\\+00 ")
    message(FATAL_ERROR "${file}: expected (11|11) = ${u}, found '${line}'")
  endif()
endfunction()

file(WRITE "${test_dir}/own" "${common}"
  "{ U=2; FCIDUMP=\"u2.fcidump\" }\n{ U=8; FCIDUMP=\"u8.fcidump\" }\n")
run_fulldiag(own)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "fulldiag failed with its own FCIDUMP per task:\n${output}")
endif()
require_u(u2.fcidump 2)
require_u(u8.fcidump 8)

file(WRITE "${test_dir}/shared" "${common}" "FCIDUMP=\"shared.fcidump\"\n{ U=2 }\n{ U=8 }\n")
run_fulldiag(shared)
if(status EQUAL 0)
  message(FATAL_ERROR "fulldiag succeeded with a shared FCIDUMP:\n${output}")
endif()
if(NOT output MATCHES "shared.fcidump is also the FCIDUMP of [^\n]*shared.task1.in.xml; give each task its own FCIDUMP")
  message(FATAL_ERROR "the second task's FCIDUMP was not refused:\n${output}")
endif()
string(REGEX MATCHALL "Wrote FCIDUMP" writes "${output}")
list(LENGTH writes writes)
if(NOT writes EQUAL 0 OR EXISTS "${test_dir}/shared.fcidump")
  message(FATAL_ERROR "a task ran although the FCIDUMP was shared:\n${output}")
endif()
