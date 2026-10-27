# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#
# Helpers for CMake-script tests (cmake -P) that run a Monte Carlo
# application and compare its averages with exactly known values.
# CMake's math() is integer-only, so numbers are compared in units of 1e-9.

# Run a command in ${test_dir}; fail the test with its output if it fails.
# The combined output is left in ${output} in the caller's scope.
function(alps_run name)
  execute_process(COMMAND ${ARGN}
    WORKING_DIRECTORY "${test_dir}" OUTPUT_VARIABLE out ERROR_VARIABLE err
    RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${name} failed (${status}):\n${out}${err}")
  endif()
  set(output "${out}${err}" PARENT_SCOPE)
endfunction()

# Convert a decimal number, such as -0.3454579 or 4.5e-05, to an integer
# number of 1e-9 units (truncated).
function(alps_to_nano var text)
  if(NOT text MATCHES "^([-+]?)([0-9]*)\\.?([0-9]*)([eE]([-+]?[0-9]+))?$")
    message(FATAL_ERROR "not a number: '${text}'")
  endif()
  set(sign "${CMAKE_MATCH_1}")
  set(digits "${CMAKE_MATCH_2}${CMAKE_MATCH_3}")
  set(exponent "${CMAKE_MATCH_5}")
  if(exponent STREQUAL "")
    set(exponent 0)
  endif()
  string(REGEX REPLACE "^\\+" "" exponent "${exponent}")
  string(LENGTH "${CMAKE_MATCH_2}" point)
  math(EXPR keep "${point} + ${exponent} + 9")
  if(keep LESS_EQUAL 0)
    set(${var} 0 PARENT_SCOPE)
    return()
  endif()
  string(LENGTH "${digits}" have)
  while(have LESS keep)
    string(APPEND digits "0")
    math(EXPR have "${have} + 1")
  endwhile()
  string(SUBSTRING "${digits}" 0 ${keep} digits)
  string(REGEX REPLACE "^0+" "" digits "${digits}")
  if(digits STREQUAL "")
    set(digits 0)
  endif()
  if(sign STREQUAL "-")
    set(digits "-${digits}")
  endif()
  set(${var} "${digits}" PARENT_SCOPE)
endfunction()

# Read <SCALAR_AVERAGE name="..."> ... <MEAN>m</MEAN> ... <ERROR>e</ERROR>
# from an ALPS result file into ${mean} and ${error}.
function(alps_read_average file name)
  file(READ "${file}" xml)
  string(REGEX MATCH "<SCALAR_AVERAGE name=\"${name}\"[^>]*>.*</SCALAR_AVERAGE>" block "${xml}")
  if(NOT block)
    message(FATAL_ERROR "${file} has no average '${name}':\n${xml}")
  endif()
  # Cut at the end of the first SCALAR_AVERAGE element.
  string(REGEX REPLACE "</SCALAR_AVERAGE>.*" "" block "${block}")
  if(NOT block MATCHES "<MEAN[^>]*>([^<]*)</MEAN>")
    message(FATAL_ERROR "average '${name}' in ${file} has no mean:\n${block}")
  endif()
  set(mean "${CMAKE_MATCH_1}" PARENT_SCOPE)
  if(NOT block MATCHES "<ERROR[^>]*>([^<]*)</ERROR>")
    message(FATAL_ERROR "average '${name}' in ${file} has no error:\n${block}")
  endif()
  set(error "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

# Fail unless |mean - exact| <= sigmas * error, with error > 0.
function(alps_check_within_error what mean error exact sigmas)
  alps_to_nano(m "${mean}")
  alps_to_nano(e "${error}")
  alps_to_nano(x "${exact}")
  if(e LESS_EQUAL 0)
    message(FATAL_ERROR "${what}: error bar ${error} is not positive")
  endif()
  math(EXPR deviation "${m} - ${x}")
  if(deviation LESS 0)
    math(EXPR deviation "-(${deviation})")
  endif()
  math(EXPR allowed "${sigmas} * ${e}")
  if(deviation GREATER allowed)
    message(FATAL_ERROR
      "${what}: ${mean} +/- ${error} is more than ${sigmas} error bars from the exact ${exact}")
  endif()
  message(STATUS "${what}: ${mean} +/- ${error}, exact ${exact}")
endfunction()

# Fail unless |value - exact| <= percent/100 * |exact|.
function(alps_check_within_percent what value exact percent)
  alps_to_nano(v "${value}")
  alps_to_nano(x "${exact}")
  math(EXPR deviation "${v} - ${x}")
  if(deviation LESS 0)
    math(EXPR deviation "-(${deviation})")
  endif()
  if(x LESS 0)
    math(EXPR allowed "-(${x}) * ${percent} / 100")
  else()
    math(EXPR allowed "${x} * ${percent} / 100")
  endif()
  if(deviation GREATER allowed)
    message(FATAL_ERROR
      "${what}: ${value} deviates by more than ${percent}% from the exact ${exact}")
  endif()
  message(STATUS "${what}: ${value}, exact ${exact}")
endfunction()

# Read the <point><x>..</x><y>..</y></point> list of an ALPS plot file into
# ${xs} and ${ys}.
function(alps_read_plot file)
  file(READ "${file}" xml)
  string(REGEX MATCHALL "<point><x>[^<]*</x><y>[^<]*</y></point>" points "${xml}")
  set(xs)
  set(ys)
  foreach(point IN LISTS points)
    string(REGEX REPLACE "<point><x>([^<]*)</x><y>([^<]*)</y></point>" "\\1" x "${point}")
    string(REGEX REPLACE "<point><x>([^<]*)</x><y>([^<]*)</y></point>" "\\2" y "${point}")
    list(APPEND xs "${x}")
    list(APPEND ys "${y}")
  endforeach()
  set(xs "${xs}" PARENT_SCOPE)
  set(ys "${ys}" PARENT_SCOPE)
endfunction()
