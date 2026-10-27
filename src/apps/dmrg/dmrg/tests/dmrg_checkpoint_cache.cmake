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

# SIGUSR1 and SIGUSR2 ask ALPS programs for a checkpoint; the run goes on.
run_dmrg(usr1 "CHECKPOINT_DIRECTORY=\"${test_dir}/cache_usr1\"" "usr1:5")
expect("usr1_exit EQUAL 0" "SIGUSR1 ended the run:\n${usr1_output}")
expect("usr1_energy STREQUAL reference_energy" "SIGUSR1 changed the energy")

# Fail twice in a row: the second run resumes and fails again, the third finishes.
set(cache "${test_dir}/cache_twice")
run_dmrg(twice "CHECKPOINT_DIRECTORY=\"${cache}\"" "kill:6")
run_dmrg(twice "" "kill:3")
expect("twice_output MATCHES \"resuming from checkpoint\"" "second run did not resume")
run_dmrg(twice "" "")
expect("twice_exit EQUAL 0" "third run failed:\n${twice_output}")
expect("twice_output MATCHES \"resuming from checkpoint after sweep 2/2\"" "third run did not resume from the second")
expect("twice_energy STREQUAL reference_energy" "twice resumed energy ${twice_energy} differs from ${reference_energy}")

# A regenerated input with a new SEED, as parameter2xml and pyalps.writeInputFiles
# write it, still resumes.
set(cache "${test_dir}/cache_reseed")
run_dmrg(reseed "CHECKPOINT_DIRECTORY=\"${cache}\"" "kill:6")
file(REMOVE "${test_dir}/reseed/parameters.in.xml")
string(REPLACE "SEED=12345" "SEED=777" reseeded "${physics}")
file(WRITE "${test_dir}/reseed/parameters" "${reseeded}CHECKPOINT_DIRECTORY=\"${cache}\"\n")
execute_process(COMMAND "${parameter2xml}" parameters
  WORKING_DIRECTORY "${test_dir}/reseed" OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
run_dmrg(reseed "" "")
expect("reseed_exit EQUAL 0" "rerun with a new SEED failed:\n${reseed_output}")
expect("reseed_output MATCHES \"resuming from checkpoint\"" "rerun with a new SEED did not resume")

# A checkpoint file that is missing or truncated is never resumed from: DMTK
# would go on with an empty block and give wrong results.
foreach(damage missing truncated)
  set(cache "${test_dir}/cache_${damage}")
  run_dmrg(${damage} "CHECKPOINT_DIRECTORY=\"${cache}\"" "kill:6")
  file(GLOB blocks "${cache}/*/block_*")
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
set(cache "${test_dir}/cache_attempts")
run_dmrg(attempts "CHECKPOINT_DIRECTORY=\"${cache}\"" "kill:6")
foreach(attempt 1 2)
  run_dmrg(attempts "" "kill-staged:1")
  expect("attempts_output MATCHES \"resuming\"" "attempt ${attempt} did not resume")
endforeach()
run_dmrg(attempts "" "")
expect("attempts_exit EQUAL 0" "run after failed resumes failed:\n${attempts_output}")
expect("attempts_output MATCHES \"failed before the next checkpoint\"" "checkpoint kept after failed resumes:\n${attempts_output}")
expect("attempts_energy STREQUAL reference_energy" "fresh start after failed resumes changed the energy")

# ~ and environment variables in CHECKPOINT_DIRECTORY are expanded; an unset
# variable warns and the run goes on without checkpoints. A parameter file
# only passes ${NAME} on; XML input can also use $NAME.
set(ENV{HOME} "${test_dir}/home")
set(ENV{ALPS_DMRG_TEST_SCRATCH} "${test_dir}/scratch")
run_dmrg(tilde "CHECKPOINT_DIRECTORY=\"~/checkpoints\"" "kill:4")
entries("${test_dir}/home/checkpoints" kept)
list(LENGTH kept count)
expect("count EQUAL 1" "~ not expanded: ${tilde_output}")
run_dmrg(variable "CHECKPOINT_DIRECTORY=\"\${ALPS_DMRG_TEST_SCRATCH}/x\"" "kill:4")
entries("${test_dir}/scratch/x" kept)
list(LENGTH kept count)
expect("count EQUAL 1" "environment variables not expanded: ${variable_output}")
run_dmrg(unset "CHECKPOINT_DIRECTORY=\"\${ALPS_DMRG_TEST_UNSET}/x\"" "")
expect("unset_exit EQUAL 0" "unset variable failed the run:\n${unset_output}")
expect("unset_output MATCHES \"ALPS_DMRG_TEST_UNSET is not set\"" "no warning for an unset variable:\n${unset_output}")

# An entry of a run on another host stays while that run updates it, even
# over CHECKPOINT_MAX_SIZE, since a lock on a network filesystem may not be
# seen from here.
set(cache "${test_dir}/cache_remote")
file(MAKE_DIRECTORY "${cache}/remote")
file(WRITE "${cache}/remote/meta" "ALPS DMRG checkpoint 1\nhost elsewhere:0\npid 1\nfingerprint\nfake\n")
file(WRITE "${cache}/remote/data" "${payload}")
fake_entry("${cache}/local" 999999)
file(WRITE "${cache}/local/data" "${payload}")
run_dmrg(remote "CHECKPOINT_DIRECTORY=\"${cache}\"\nCHECKPOINT_MAX_SIZE=0.001" "")
expect("remote_exit EQUAL 0" "run next to a remote entry failed:\n${remote_output}")
expect("EXISTS \"${cache}/remote\"" "entry of a live run on another host was removed")
expect("NOT EXISTS \"${cache}/local\"" "dead local entry was kept over the size cap")
