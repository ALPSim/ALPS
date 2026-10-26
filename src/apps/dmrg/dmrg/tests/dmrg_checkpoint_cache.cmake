# Failed runs keep a checkpoint that a rerun resumes from; stale entries are
# collected. ALPS_DMRG_CHECKPOINT_TEST simulates failures at a given commit.
# An 8-site chain with 2 sweeps commits 16 times (8 iterations per sweep).
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}/xml")
file(COPY "${xml_dir}/lattices.xml" "${xml_dir}/models.xml"
  "${xml_dir}/ALPS.xsl" DESTINATION "${test_dir}/xml")
set(ENV{ALPS_XML_PATH} "${test_dir}/xml")
cmake_host_system_information(RESULT host QUERY HOSTNAME)

set(physics [=[
LATTICE="open chain lattice"
MODEL="spin"
CONSERVED_QUANTUMNUMBERS="N,Sz"
Sz_total=0
L=8
J=1
SWEEPS=2
NUMBER_EIGENVALUES=1
MAXSTATES=20
SEED=12345
]=])

# Run dmrg in ${test_dir}/<name> with extra parameter lines; sets <name>_exit,
# <name>_output and <name>_energy. With the same SEED, a resumed run gives the
# same energy as an uninterrupted one to the last digit.
function(run_dmrg name extra hook)
  set(dir "${test_dir}/${name}")
  if(NOT EXISTS "${dir}/parameters.in.xml")
    file(MAKE_DIRECTORY "${dir}")
    file(WRITE "${dir}/parameters" "${physics}${extra}\n")
    execute_process(COMMAND "${parameter2xml}" parameters
      WORKING_DIRECTORY "${dir}" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
  endif()
  if(hook)
    set(ENV{ALPS_DMRG_CHECKPOINT_TEST} "${hook}")
  else()
    unset(ENV{ALPS_DMRG_CHECKPOINT_TEST})
  endif()
  execute_process(COMMAND "${dmrg}" --write-xml parameters.in.xml
    WORKING_DIRECTORY "${dir}" RESULT_VARIABLE result
    OUTPUT_VARIABLE out ERROR_VARIABLE err)
  unset(ENV{ALPS_DMRG_CHECKPOINT_TEST})
  set(${name}_exit "${result}" PARENT_SCOPE)
  set(${name}_output "${out}${err}" PARENT_SCOPE)
  set(energy "")
  if(EXISTS "${dir}/parameters.task1.out.xml")
    file(READ "${dir}/parameters.task1.out.xml" xml)
    if(xml MATCHES "name=\"Energy\"><MEAN>([^<]+)")
      set(energy "${CMAKE_MATCH_1}")
    endif()
  endif()
  set(${name}_energy "${energy}" PARENT_SCOPE)
endfunction()

function(expect condition message)
  cmake_language(EVAL CODE "if(NOT (${condition}))\n  message(FATAL_ERROR [==[${message}]==])\nendif()")
endfunction()

function(entries cache out)
  file(GLOB found LIST_DIRECTORIES true "${cache}/*")
  set(${out} "${found}" PARENT_SCOPE)
endfunction()

function(fake_entry dir pid)
  file(MAKE_DIRECTORY "${dir}")
  file(WRITE "${dir}/meta" "ALPS DMRG checkpoint 1\nhost ${host}\npid ${pid}\nfingerprint\nfake\n")
endfunction()

# Reference without checkpoints; nothing may be written to a cache.
run_dmrg(reference "CHECKPOINT=0\nCHECKPOINT_DIRECTORY=\"${test_dir}/unused\"" "")
expect("reference_exit EQUAL 0" "reference run failed:\n${reference_output}")
expect("NOT reference_energy STREQUAL \"\"" "reference energy missing")
expect("NOT EXISTS \"${test_dir}/unused\"" "CHECKPOINT=0 created a checkpoint directory")

# A successful run removes its entry.
run_dmrg(success "CHECKPOINT_DIRECTORY=\"${test_dir}/cache_success\"" "")
expect("success_exit EQUAL 0" "checkpointed run failed:\n${success_output}")
entries("${test_dir}/cache_success" left)
expect("left STREQUAL \"\"" "successful run left its checkpoint: ${left}")
expect("success_energy STREQUAL reference_energy" "checkpointed run changed the energy")

# Each failure leaves an entry that the rerun resumes from and then removes;
# kill-staged dies with uncommitted files, kill:16 after the last commit.
foreach(hook throw:3 kill:6 kill-staged:5 term:10 kill:16)
  string(REPLACE ":" "_" name "${hook}")
  string(REPLACE "-" "_" name "${name}")
  set(cache "${test_dir}/cache_${name}")
  run_dmrg(${name} "CHECKPOINT_DIRECTORY=\"${cache}\"" "${hook}")
  expect("NOT ${name}_exit EQUAL 0" "${hook}: the failing run succeeded")
  entries("${cache}" kept)
  list(LENGTH kept count)
  expect("count EQUAL 1" "${hook}: expected one checkpoint, found: ${kept}")
  expect("EXISTS \"${kept}/progress\"" "${hook}: checkpoint has no progress record")
  if(hook MATCHES "^(throw|term)")
    expect("${name}_output MATCHES \"Checkpoint kept at\"" "${hook}: no failure message:\n${${name}_output}")
  endif()
  set(failed_output "${${name}_output}")
  run_dmrg(${name} "" "")
  expect("${name}_exit EQUAL 0" "${hook}: rerun failed:\n${${name}_output}")
  expect("${name}_output MATCHES \"resuming from checkpoint\"" "${hook}: rerun did not resume")
  expect("${name}_energy STREQUAL reference_energy"
    "${hook}: resumed energy ${${name}_energy} differs from ${reference_energy}")
  entries("${cache}" left)
  expect("left STREQUAL \"\"" "${hook}: resumed run left its checkpoint")
endforeach()

# A run with other parameters starts fresh and leaves the other entry alone.
set(cache "${test_dir}/cache_params")
run_dmrg(params_a "CHECKPOINT_DIRECTORY=\"${cache}\"" "kill:4")
run_dmrg(params_b "h=0.1\nCHECKPOINT_DIRECTORY=\"${cache}\"" "")
expect("params_b_exit EQUAL 0" "run with other parameters failed:\n${params_b_output}")
expect("NOT params_b_output MATCHES \"resuming\"" "resumed from a checkpoint of other parameters")
entries("${cache}" kept)
list(LENGTH kept count)
expect("count EQUAL 1" "checkpoint of other parameters was removed")

# An entry from another format version is discarded, not resumed.
set(cache "${test_dir}/cache_version")
run_dmrg(version "CHECKPOINT_DIRECTORY=\"${cache}\"" "kill:4")
entries("${cache}" kept)
file(READ "${kept}/meta" meta)
string(REPLACE "ALPS DMRG checkpoint 1" "ALPS DMRG checkpoint 0" meta "${meta}")
file(WRITE "${kept}/meta" "${meta}")
run_dmrg(version "" "")
expect("version_exit EQUAL 0" "run after version change failed:\n${version_output}")
expect("version_output MATCHES \"discarding checkpoint\"" "old-version checkpoint not discarded")
expect("NOT version_output MATCHES \"resuming\"" "resumed from an old-version checkpoint")
expect("version_energy STREQUAL reference_energy" "fresh start after version change changed the energy")

# Garbage collection: expired entries go; fresh, live and foreign directories stay.
set(cache "${test_dir}/cache_gc")
fake_entry("${cache}/expired" 999999)
fake_entry("${cache}/fresh" 999999)
fake_entry("${cache}/live" 1)  # pid 1 always exists
file(MAKE_DIRECTORY "${cache}/not_an_entry")
file(WRITE "${cache}/not_an_entry/meta" "something else\n")
foreach(old expired live not_an_entry)
  execute_process(COMMAND touch -t 202001010000 "${cache}/${old}/meta" COMMAND_ERROR_IS_FATAL ANY)
endforeach()
run_dmrg(gc "CHECKPOINT_DIRECTORY=\"${cache}\"" "")
expect("gc_exit EQUAL 0" "run with stale entries failed:\n${gc_output}")
expect("NOT EXISTS \"${cache}/expired\"" "expired checkpoint was not removed")
foreach(kept fresh live not_an_entry)
  expect("EXISTS \"${cache}/${kept}\"" "garbage collection removed ${kept}")
endforeach()

# Over CHECKPOINT_MAX_SIZE the oldest entry goes first.
set(cache "${test_dir}/cache_size")
string(REPEAT "x" 600000 payload)
fake_entry("${cache}/older" 999999)
file(WRITE "${cache}/older/data" "${payload}")
execute_process(COMMAND "${CMAKE_COMMAND}" -E sleep 1)
fake_entry("${cache}/newer" 999999)
file(WRITE "${cache}/newer/data" "${payload}")
run_dmrg(size "CHECKPOINT_DIRECTORY=\"${cache}\"\nCHECKPOINT_MAX_SIZE=0.001" "")
expect("size_exit EQUAL 0" "run over the size cap failed:\n${size_output}")
expect("NOT EXISTS \"${cache}/older\"" "oldest checkpoint not evicted over the size cap")
expect("EXISTS \"${cache}/newer\"" "size cap evicted more than needed")

# An unusable directory never fails the run; without a terminal there is no prompt.
run_dmrg(unusable "CHECKPOINT_DIRECTORY=\"/dev/null/checkpoints\"" "")
expect("unusable_exit EQUAL 0" "unusable checkpoint directory failed the run:\n${unusable_output}")
expect("unusable_output MATCHES \"continuing without checkpoints\"" "no warning for an unusable directory")
expect("unusable_energy STREQUAL reference_energy" "run without checkpoints changed the energy")

# A relative CHECKPOINT_DIRECTORY is resolved against the working directory.
run_dmrg(relative "CHECKPOINT_DIRECTORY=\"here\"" "kill:4")
entries("${test_dir}/relative/here" kept)
list(LENGTH kept count)
expect("count EQUAL 1" "relative CHECKPOINT_DIRECTORY not used: ${relative_output}")
