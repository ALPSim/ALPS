# Failed runs keep a checkpoint next to their task files that a rerun with
# RESUME=1 resumes from; a rerun without it refuses to start.
# ALPS_DMRG_CHECKPOINT_TEST simulates failures at a given commit.
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

# Rerun an existing run with RESUME=1 (and the given extra parameter lines).
function(set_resume name)
  change_parameters(${name} "${physics}RESUME=1\n${ARGN}")
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

# Each failure leaves a checkpoint. A rerun without RESUME refuses to start and
# leaves it alone; one with RESUME=1 resumes from it and then removes it.
# kill-staged dies with uncommitted files, kill:16 after the last commit.
foreach(hook throw:3 kill:6 kill-staged:5 term:10 kill:16)
  string(REPLACE ":" "_" name "${hook}")
  string(REPLACE "-" "_" name "${name}")
  run_dmrg(${name} "" "${hook}")
  expect("NOT ${name}_exit EQUAL 0" "${hook}: the failing run succeeded")
  expect("EXISTS \"${${name}_checkpoint}/progress\"" "${hook}: no checkpoint kept")
  if(hook MATCHES "^(throw|term)")
    expect("${name}_output MATCHES \"Checkpoint kept in\"" "${hook}: no failure message:\n${${name}_output}")
  endif()
  file(READ "${${name}_checkpoint}/progress" progress_before)
  run_dmrg(${name} "" "")
  expect("NOT ${name}_exit EQUAL 0" "${hook}: rerun without RESUME started")
  expect("${name}_output MATCHES \"Set RESUME=1\"" "${hook}: no RESUME hint:\n${${name}_output}")
  file(READ "${${name}_checkpoint}/progress" progress_after)
  expect("progress_after STREQUAL progress_before" "${hook}: rerun without RESUME changed the checkpoint")
  set_resume(${name})
  run_dmrg(${name} "" "")
  expect("${name}_exit EQUAL 0" "${hook}: rerun failed:\n${${name}_output}")
  expect("${name}_output MATCHES \"resuming from checkpoint\"" "${hook}: rerun did not resume")
  expect("${name}_energy STREQUAL reference_energy"
    "${hook}: resumed energy ${${name}_energy} differs from ${reference_energy}")
  expect("NOT EXISTS \"${${name}_checkpoint}\"" "${hook}: resumed run left its checkpoint")
endforeach()

# SIGUSR1 and SIGUSR2 ask ALPS programs for a checkpoint; the run goes on,
# also after a second one (the ALPS handlers are one-shot).
run_dmrg(usr1 "" "usr1:5")
expect("usr1_exit EQUAL 0" "SIGUSR1 ended the run:\n${usr1_output}")
expect("usr1_energy STREQUAL reference_energy" "SIGUSR1 changed the energy")
run_dmrg(usr1_twice "" "usr1-twice:5")
expect("usr1_twice_exit EQUAL 0" "a second SIGUSR1 ended the run:\n${usr1_twice_output}")
expect("usr1_twice_energy STREQUAL reference_energy" "two SIGUSR1 changed the energy")

# Fail twice in a row: the second run resumes and fails again, the third finishes.
run_dmrg(twice "" "kill:6")
set_resume(twice)
run_dmrg(twice "" "kill:3")
expect("twice_output MATCHES \"resuming from checkpoint\"" "second run did not resume")
run_dmrg(twice "" "")
expect("twice_exit EQUAL 0" "third run failed:\n${twice_output}")
expect("twice_output MATCHES \"resuming from checkpoint after sweep 2/2\"" "third run did not resume from the second")
expect("twice_energy STREQUAL reference_energy" "twice resumed energy ${twice_energy} differs from ${reference_energy}")

# A regenerated input with a new SEED, as parameter2xml and pyalps.writeInputFiles
# write it, still resumes. The run continues with the SEED it started with,
# so it gives the reference results and its output lists that SEED.
run_dmrg(reseed "" "kill:6")
string(REPLACE "SEED=12345" "SEED=777" reseeded "${physics}")
change_parameters(reseed "${reseeded}RESUME=1")
run_dmrg(reseed "" "")
expect("reseed_exit EQUAL 0" "rerun with a new SEED failed:\n${reseed_output}")
expect("reseed_output MATCHES \"resuming from checkpoint\"" "rerun with a new SEED did not resume")
expect("reseed_output MATCHES \"started with SEED=12345\"" "original SEED not reported:\n${reseed_output}")
expect("reseed_energy STREQUAL reference_energy" "reseeded energy ${reseed_energy} differs from ${reference_energy}")
file(READ "${test_dir}/reseed/parameters.task1.out.xml" reseed_xml)
expect("reseed_xml MATCHES \"name=\\\"SEED\\\">12345<\"" "output does not list the SEED the results follow")

# RESUME=1 never starts over: a checkpoint it cannot resume fails the run and
# is left for the user to delete.
macro(expect_refused name what pattern)
  expect("NOT ${name}_exit EQUAL 0" "${what}: the run started")
  expect("${name}_output MATCHES \"${pattern}\"" "${what}: wrong message:\n${${name}_output}")
  expect("NOT ${name}_output MATCHES \"resuming\"" "${what}: resumed")
endmacro()

# RESUME=1 can stay in the parameters: a task without a checkpoint starts.
run_dmrg(nothing "RESUME=1" "")
expect("nothing_exit EQUAL 0" "RESUME=1 without a checkpoint failed:\n${nothing_output}")
expect("nothing_output MATCHES \"it starts from the beginning\"" "no notice for RESUME=1 without a checkpoint")
expect("nothing_energy STREQUAL reference_energy" "RESUME=1 without a checkpoint changed the energy")
expect("NOT EXISTS \"${nothing_checkpoint}\"" "RESUME=1 without a checkpoint left a directory")

# A job of several tasks killed in task 2 finishes with RESUME=1 set once for
# the job: tasks without a checkpoint start, task 2 resumes.
file(MAKE_DIRECTORY "${test_dir}/multitask")
string(REPLACE "L=8\n" "" job "${physics}")
set(tasks "{ L=6 }\n{ L=8 }\n{ L=6, J=0.5 }\n")
change_parameters(multitask "${job}${tasks}")
run_dmrg(multitask "" "kill:12")
expect("EXISTS \"${test_dir}/multitask/parameters.task2.checkpoint/progress\"" "multi-task job: no checkpoint for task 2")
file(REMOVE "${test_dir}/multitask/parameters.out.xml")
# a parameter after the task blocks would apply to no task
change_parameters(multitask "${job}RESUME=1\n${tasks}")
run_dmrg(multitask "" "")
expect("multitask_exit EQUAL 0" "multi-task job with RESUME=1 failed:\n${multitask_output}")
expect("multitask_output MATCHES \"resuming from checkpoint\"" "task 2 did not resume")
foreach(task 1 2 3)
  expect("EXISTS \"${test_dir}/multitask/parameters.task${task}.out.xml\"" "multi-task job: task ${task} has no results")
endforeach()

# Changed parameters: not resumed, with or without RESUME, and kept.
run_dmrg(params "" "kill:4")
change_parameters(params "${physics}h=0.1")
run_dmrg(params "" "")
expect_refused(params "other parameters without RESUME" "Set RESUME=1")
change_parameters(params "${physics}h=0.1\nRESUME=1")
run_dmrg(params "" "")
expect_refused(params "other parameters" "belongs to other parameters")
expect("EXISTS \"${params_checkpoint}/progress\"" "checkpoint of other parameters was removed")

# A checkpoint of another format version. Version 1 is the format before
# BMatrix files lost their unused hash members.
run_dmrg(version "" "kill:4")
file(READ "${version_checkpoint}/meta" meta)
string(REPLACE "ALPS DMRG checkpoint 2" "ALPS DMRG checkpoint 1" meta "${meta}")
file(WRITE "${version_checkpoint}/meta" "${meta}")
set_resume(version)
run_dmrg(version "" "")
expect_refused(version "old-version checkpoint" "unreadable or from another version")
expect("EXISTS \"${version_checkpoint}/progress\"" "old-version checkpoint was removed")

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
  set_resume(${damage})
  run_dmrg(${damage} "" "")
  expect_refused(${damage} "${damage} file" "it is damaged")
  expect("EXISTS \"${${damage}_checkpoint}/progress\"" "${damage} file: checkpoint was removed")
endforeach()

# A checkpoint that every resumed run fails on before completing an iteration
# is refused after two attempts.
run_dmrg(attempts "" "kill:6")
set_resume(attempts)
foreach(attempt 1 2)
  run_dmrg(attempts "" "kill-iter:1")
  expect("attempts_output MATCHES \"resuming\"" "attempt ${attempt} did not resume")
endforeach()
run_dmrg(attempts "" "")
expect_refused(attempts "after failed resumes" "failed before completing an iteration")

# Runs killed after completing an iteration are not counted (#371): three
# kills in the final sweep, which is never checkpointed, still resume.
run_dmrg(final_kills "" "kill:16")
set_resume(final_kills)
foreach(attempt 1 2 3)
  run_dmrg(final_kills "" "kill-iter:2")
  expect("NOT final_kills_exit EQUAL 0" "kill ${attempt} in the final sweep did not happen")
endforeach()
run_dmrg(final_kills "" "")
expect("final_kills_exit EQUAL 0" "resume after kills in the final sweep failed:\n${final_kills_output}")
expect("final_kills_energy STREQUAL reference_energy" "final_kills energy ${final_kills_energy} differs")

# RESUME needs checkpoints and excludes START_SWEEP.
run_dmrg(resume_off "RESUME=1\nCHECKPOINT=0" "")
expect_refused(resume_off "RESUME=1 with CHECKPOINT=0" "RESUME=1 needs checkpoints")
run_dmrg(resume_start "RESUME=1\nSTART_SWEEP=1" "")
expect_refused(resume_start "RESUME=1 with START_SWEEP" "RESUME and START_SWEEP cannot be used together")

# SCRATCH_IN_MEMORY keeps the files in memory and writes them to the checkpoint
# only on SIGTERM, SIGINT, SIGQUIT (and then ends the run), SIGUSR1 or SIGUSR2.
run_dmrg(memory "SCRATCH_IN_MEMORY=1" "memory-term:5")
expect("NOT memory_exit EQUAL 0" "SIGTERM did not end the in-memory run")
expect("EXISTS \"${memory_checkpoint}/progress\"" "SIGTERM: in-memory files not written to the checkpoint")
foreach(lanczos iter_ALPS.dat vectors_ALPS.dat)
  expect("NOT EXISTS \"${memory_checkpoint}/${lanczos}\"" "in-memory run wrote ${lanczos} to the checkpoint")
endforeach()
change_parameters(memory "${physics}SCRATCH_IN_MEMORY=1\nRESUME=1")
run_dmrg(memory "" "")
expect("memory_exit EQUAL 0" "in-memory resume failed:\n${memory_output}")
expect("memory_output MATCHES \"resuming from checkpoint\"" "in-memory run did not resume")
expect("memory_energy STREQUAL reference_energy" "in-memory resumed energy ${memory_energy} differs")
expect("NOT EXISTS \"${memory_checkpoint}\"" "in-memory resumed run left its checkpoint")
run_dmrg(memory_quiet "SCRATCH_IN_MEMORY=1" "")
expect("memory_quiet_exit EQUAL 0" "in-memory run failed:\n${memory_quiet_output}")
expect("NOT memory_quiet_output MATCHES \"changed in-memory files to the checkpoint\"" "in-memory run wrote a checkpoint unasked")
run_dmrg(memory_usr1 "SCRATCH_IN_MEMORY=1" "memory-usr1:3")
expect("memory_usr1_exit EQUAL 0" "SIGUSR1 ended the in-memory run:\n${memory_usr1_output}")
expect("memory_usr1_output MATCHES \"changed in-memory files to the checkpoint\"" "SIGUSR1: no checkpoint written")
expect("memory_usr1_energy STREQUAL reference_energy" "SIGUSR1 changed the in-memory energy")
# The checkpoint of an in-memory run resumes on disk, and the other way round.
run_dmrg(memory_to_disk "SCRATCH_IN_MEMORY=1" "memory-term:5")
set_resume(memory_to_disk)
run_dmrg(memory_to_disk "" "")
expect("memory_to_disk_energy STREQUAL reference_energy" "in-memory checkpoint resumed on disk: ${memory_to_disk_energy}")
run_dmrg(disk_to_memory "" "throw:5")
set_resume(disk_to_memory "SCRATCH_IN_MEMORY=1")
run_dmrg(disk_to_memory "" "")
expect("disk_to_memory_energy STREQUAL reference_energy" "disk checkpoint resumed in memory: ${disk_to_memory_energy}")

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

# A directory that is in the way but is not a checkpoint is left alone.
file(MAKE_DIRECTORY "${test_dir}/foreign/parameters.task1.checkpoint")
file(WRITE "${test_dir}/foreign/parameters.task1.checkpoint/notes.txt" "keep me\n")
run_dmrg(foreign "" "")
expect("foreign_exit EQUAL 0" "foreign directory failed the run:\n${foreign_output}")
expect("foreign_output MATCHES \"is not a DMRG checkpoint\"" "no warning for a foreign directory:\n${foreign_output}")
expect("EXISTS \"${foreign_checkpoint}/notes.txt\"" "contents of a foreign directory were deleted")
expect("foreign_energy STREQUAL reference_energy" "run next to a foreign directory changed the energy")
# parameter2xml does not overwrite the inputs of a finished job
file(REMOVE "${test_dir}/foreign/parameters.out.xml")
change_parameters(foreign "${physics}RESUME=1")
run_dmrg(foreign "" "")
expect("foreign_exit EQUAL 0" "RESUME=1 next to a foreign directory failed:\n${foreign_output}")
expect("foreign_output MATCHES \"is not a DMRG checkpoint\"" "RESUME=1: no warning for a foreign directory")
expect("EXISTS \"${foreign_checkpoint}/notes.txt\"" "RESUME=1 deleted the contents of a foreign directory")

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

# KEEP_CHECKPOINT keeps the checkpoint of a successful run, and RESUME=1 with a
# larger SWEEPS continues it. With a constant number of states, one sweep
# continued by a second gives the energy of a two-sweep run.
string(REPLACE "MAXSTATES=20" "NUMSTATES=20" constant "${physics}")
string(REPLACE "SWEEPS=2" "SWEEPS=1" one_sweep "${constant}")
file(MAKE_DIRECTORY "${test_dir}/keep_reference")
change_parameters(keep_reference "${constant}CHECKPOINT=0")
run_dmrg(keep_reference "" "")
expect("keep_reference_exit EQUAL 0" "two-sweep reference failed:\n${keep_reference_output}")
foreach(mode disk memory)
  if(mode STREQUAL "memory")
    set(scratch "SCRATCH_IN_MEMORY=1\n")
  else()
    set(scratch "")
  endif()
  set(name keep_${mode})
  file(MAKE_DIRECTORY "${test_dir}/${name}")
  change_parameters(${name} "${one_sweep}${scratch}KEEP_CHECKPOINT=1")
  run_dmrg(${name} "" "")
  expect("${name}_exit EQUAL 0" "${mode}: KEEP_CHECKPOINT run failed:\n${${name}_output}")
  expect("${name}_output MATCHES \"checkpoint kept in\"" "${mode}: no notice of the kept checkpoint")
  expect("EXISTS \"${${name}_checkpoint}/progress\"" "${mode}: KEEP_CHECKPOINT kept no checkpoint")
  file(REMOVE "${test_dir}/${name}/parameters.out.xml")
  change_parameters(${name} "${constant}${scratch}RESUME=1")
  run_dmrg(${name} "" "")
  expect("${name}_exit EQUAL 0" "${mode}: continuation failed:\n${${name}_output}")
  expect("${name}_output MATCHES \"resuming from checkpoint after sweep 1/1,\"" "${mode}: continuation did not resume:\n${${name}_output}")
  expect("${name}_energy STREQUAL keep_reference_energy"
    "${mode}: continued energy ${${name}_energy} differs from ${keep_reference_energy}")
  expect("NOT EXISTS \"${${name}_checkpoint}\"" "${mode}: continuation without KEEP_CHECKPOINT left its checkpoint")
endforeach()
