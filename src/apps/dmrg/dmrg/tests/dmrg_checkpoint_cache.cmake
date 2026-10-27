# Failed runs keep a checkpoint next to their task files that a rerun resumes
# from. ALPS_DMRG_CHECKPOINT_TEST simulates failures at a given commit.
# An 8-site chain with 2 sweeps commits 16 times (8 iterations per sweep).
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}/xml")
file(COPY "${xml_dir}/lattices.xml" "${xml_dir}/models.xml"
  "${xml_dir}/ALPS.xsl" DESTINATION "${test_dir}/xml")
set(ENV{ALPS_XML_PATH} "${test_dir}/xml")

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
# <name>_output, <name>_energy and <name>_checkpoint (its checkpoint
# directory). With the same SEED, a resumed run gives the same energy as an
# uninterrupted one to the last digit.
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
  set(${name}_checkpoint "${dir}/parameters.task1.checkpoint" PARENT_SCOPE)
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

# Rewrite the parameters of an existing run and regenerate its input.
function(change_parameters name parameters)
  set(dir "${test_dir}/${name}")
  file(REMOVE "${dir}/parameters.in.xml")
  file(WRITE "${dir}/parameters" "${parameters}\n")
  execute_process(COMMAND "${parameter2xml}" parameters
    WORKING_DIRECTORY "${dir}" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
endfunction()

# Reference without checkpoints; nothing may be written next to the task.
run_dmrg(reference "CHECKPOINT=0" "")
expect("reference_exit EQUAL 0" "reference run failed:\n${reference_output}")
expect("NOT reference_energy STREQUAL \"\"" "reference energy missing")
expect("NOT EXISTS \"${reference_checkpoint}\"" "CHECKPOINT=0 created a checkpoint")

# A successful run removes its checkpoint.
run_dmrg(success "" "")
expect("success_exit EQUAL 0" "checkpointed run failed:\n${success_output}")
expect("success_output MATCHES \"parameters.task1.checkpoint\"" "checkpoint not next to the task files:\n${success_output}")
expect("NOT EXISTS \"${success_checkpoint}\"" "successful run left its checkpoint")
expect("success_energy STREQUAL reference_energy" "checkpointed run changed the energy")

# Each failure leaves a checkpoint that the rerun resumes from and then
# removes; kill-staged dies with uncommitted files, kill:16 after the last commit.
foreach(hook throw:3 kill:6 kill-staged:5 term:10 kill:16)
  string(REPLACE ":" "_" name "${hook}")
  string(REPLACE "-" "_" name "${name}")
  run_dmrg(${name} "" "${hook}")
  expect("NOT ${name}_exit EQUAL 0" "${hook}: the failing run succeeded")
  expect("EXISTS \"${${name}_checkpoint}/progress\"" "${hook}: no checkpoint kept")
  if(hook MATCHES "^(throw|term)")
    expect("${name}_output MATCHES \"Checkpoint kept in\"" "${hook}: no failure message:\n${${name}_output}")
  endif()
  run_dmrg(${name} "" "")
  expect("${name}_exit EQUAL 0" "${hook}: rerun failed:\n${${name}_output}")
  expect("${name}_output MATCHES \"resuming from checkpoint\"" "${hook}: rerun did not resume")
  expect("${name}_energy STREQUAL reference_energy"
    "${hook}: resumed energy ${${name}_energy} differs from ${reference_energy}")
  expect("NOT EXISTS \"${${name}_checkpoint}\"" "${hook}: resumed run left its checkpoint")
endforeach()

# SIGUSR1 and SIGUSR2 ask ALPS programs for a checkpoint; the run goes on.
run_dmrg(usr1 "" "usr1:5")
expect("usr1_exit EQUAL 0" "SIGUSR1 ended the run:\n${usr1_output}")
expect("usr1_energy STREQUAL reference_energy" "SIGUSR1 changed the energy")

# Fail twice in a row: the second run resumes and fails again, the third finishes.
run_dmrg(twice "" "kill:6")
run_dmrg(twice "" "kill:3")
expect("twice_output MATCHES \"resuming from checkpoint\"" "second run did not resume")
run_dmrg(twice "" "")
expect("twice_exit EQUAL 0" "third run failed:\n${twice_output}")
expect("twice_output MATCHES \"resuming from checkpoint after sweep 2/2\"" "third run did not resume from the second")
expect("twice_energy STREQUAL reference_energy" "twice resumed energy ${twice_energy} differs from ${reference_energy}")

# A regenerated input with a new SEED, as parameter2xml and pyalps.writeInputFiles
# write it, still resumes.
run_dmrg(reseed "" "kill:6")
string(REPLACE "SEED=12345" "SEED=777" reseeded "${physics}")
change_parameters(reseed "${reseeded}")
run_dmrg(reseed "" "")
expect("reseed_exit EQUAL 0" "rerun with a new SEED failed:\n${reseed_output}")
expect("reseed_output MATCHES \"resuming from checkpoint\"" "rerun with a new SEED did not resume")

# Changed parameters start over.
run_dmrg(params "" "kill:4")
change_parameters(params "${physics}h=0.1")
run_dmrg(params "" "")
expect("params_exit EQUAL 0" "run with other parameters failed:\n${params_output}")
expect("params_output MATCHES \"it belongs to other parameters\"" "checkpoint of other parameters not discarded")
expect("NOT params_output MATCHES \"resuming\"" "resumed from a checkpoint of other parameters")
expect("NOT EXISTS \"${params_checkpoint}\"" "run with other parameters left a checkpoint")

# A checkpoint of another format version is discarded, not resumed.
run_dmrg(version "" "kill:4")
file(READ "${version_checkpoint}/meta" meta)
string(REPLACE "ALPS DMRG checkpoint 1" "ALPS DMRG checkpoint 0" meta "${meta}")
file(WRITE "${version_checkpoint}/meta" "${meta}")
run_dmrg(version "" "")
expect("version_exit EQUAL 0" "run after version change failed:\n${version_output}")
expect("version_output MATCHES \"discarding checkpoint\"" "old-version checkpoint not discarded")
expect("NOT version_output MATCHES \"resuming\"" "resumed from an old-version checkpoint")
expect("version_energy STREQUAL reference_energy" "fresh start after version change changed the energy")

# A checkpoint file that is missing or truncated is never resumed from: DMTK
# would go on with an empty block and give wrong results.
foreach(damage missing truncated)
  run_dmrg(${damage} "" "kill:6")
  file(GLOB blocks "${${damage}_checkpoint}/block_*")
  list(GET blocks 0 victim)
  if(damage STREQUAL "missing")
    file(REMOVE "${victim}")
  else()
    file(WRITE "${victim}" "short")
  endif()
  run_dmrg(${damage} "" "")
  expect("${damage}_exit EQUAL 0" "${damage} file: rerun failed:\n${${damage}_output}")
  expect("${damage}_output MATCHES \"discarding checkpoint\"" "${damage} file: checkpoint not discarded:\n${${damage}_output}")
  expect("NOT ${damage}_output MATCHES \"resuming\"" "${damage} file: resumed from a damaged checkpoint")
  expect("${damage}_energy STREQUAL reference_energy" "${damage} file: energy ${${damage}_energy} differs")
endforeach()

# A checkpoint that every resumed run fails on is dropped after two attempts.
run_dmrg(attempts "" "kill:6")
foreach(attempt 1 2)
  run_dmrg(attempts "" "kill-staged:1")
  expect("attempts_output MATCHES \"resuming\"" "attempt ${attempt} did not resume")
endforeach()
run_dmrg(attempts "" "")
expect("attempts_exit EQUAL 0" "run after failed resumes failed:\n${attempts_output}")
expect("attempts_output MATCHES \"failed before the next checkpoint\"" "checkpoint kept after failed resumes:\n${attempts_output}")
expect("attempts_energy STREQUAL reference_energy" "fresh start after failed resumes changed the energy")

# A run given TEMP_DIRECTORY keeps its files there, unless CHECKPOINT is set.
file(MAKE_DIRECTORY "${test_dir}/scratch")
run_dmrg(temp "TEMP_DIRECTORY=\"${test_dir}/scratch\"" "kill:4")
expect("NOT EXISTS \"${temp_checkpoint}\"" "TEMP_DIRECTORY run wrote a checkpoint")
run_dmrg(temp_on "TEMP_DIRECTORY=\"${test_dir}/scratch\"\nCHECKPOINT=1" "kill:4")
expect("EXISTS \"${temp_on_checkpoint}/progress\"" "CHECKPOINT=1 with TEMP_DIRECTORY kept no checkpoint")

# A checkpoint location that cannot be used never fails the run.
file(MAKE_DIRECTORY "${test_dir}/blocked")
file(WRITE "${test_dir}/blocked/parameters.task1.checkpoint" "a file in the way\n")
run_dmrg(blocked "" "")
expect("blocked_exit EQUAL 0" "unusable checkpoint location failed the run:\n${blocked_output}")
expect("blocked_output MATCHES \"continues without checkpoints\"" "no warning for an unusable location:\n${blocked_output}")
expect("blocked_energy STREQUAL reference_energy" "run without checkpoints changed the energy")

# Two runs of the same task at once (execute_process starts its commands
# together): the one that does not hold the checkpoint goes on without it.
# The holder keeps it for 3 seconds at its first commit.
file(MAKE_DIRECTORY "${test_dir}/concurrent")
file(WRITE "${test_dir}/concurrent/parameters" "${physics}\n")
execute_process(COMMAND "${parameter2xml}" parameters
  WORKING_DIRECTORY "${test_dir}/concurrent" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
set(ENV{ALPS_DMRG_CHECKPOINT_TEST} "sleep:1")
execute_process(
  COMMAND "${dmrg}" parameters.in.xml
  COMMAND "${dmrg}" parameters.in.xml
  WORKING_DIRECTORY "${test_dir}/concurrent" OUTPUT_VARIABLE out ERROR_VARIABLE err)
unset(ENV{ALPS_DMRG_CHECKPOINT_TEST})
expect("err MATCHES \"another run of this task is using\"" "concurrent runs shared the checkpoint:\n${out}${err}")
