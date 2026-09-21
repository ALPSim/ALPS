#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_executable(ising1 main.C ising.C)
add_executable(ising2 main2.C ising2.C)
add_executable(ising3 main3.C ising.C)

target_link_libraries(ising1 PRIVATE ALPS::alps)
target_link_libraries(ising2 PRIVATE ALPS::alps)
target_link_libraries(ising3 PRIVATE ALPS::alps)
#ADD_TEST(ising1 ising1 ./parm.in.xml)
#ADD_TEST(ising2 ising2 ./parm2.in.xml)
#ADD_TEST(ising3 ising3 ./parm3.in.xml)
