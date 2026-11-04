# convert2xml starts the scheduler without a job file to read a Monte Carlo
# task file. std::filesystem::absolute("") throws, so the empty job file
# must not be made absolute (56868297c stopped every such conversion with
# "cannot make absolute path").
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")
file(WRITE "${test_dir}/task.xml" [=[
<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>
<SIMULATION xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://xml.comp-phys.org/2002/10/ALPS.xsd">
  <PARAMETERS>
    <PARAMETER name="SWEEPS">10</PARAMETER>
  </PARAMETERS>
  <AVERAGES>
  </AVERAGES>
</SIMULATION>
]=])

execute_process(COMMAND "${convert2xml}" task.xml
  WORKING_DIRECTORY "${test_dir}" OUTPUT_VARIABLE out ERROR_VARIABLE err
  RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "convert2xml failed (${status}):\n${out}${err}")
endif()
file(READ "${test_dir}/task.xml" task)
if(NOT task MATCHES "name=\"SWEEPS\">10<")
  message(FATAL_ERROR "convert2xml lost the parameters:\n${task}")
endif()
