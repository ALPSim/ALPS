#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_executable(sparsediag sparsediag.C factory.C)
target_link_libraries(sparsediag PRIVATE ALPS::alps)
alps_install_applications(sparsediag)
