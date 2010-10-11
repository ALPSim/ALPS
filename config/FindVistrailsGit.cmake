#  Copyright Matthias Troyer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_PATCH_VISTRAILS AND NOT VISTRAILS_GIT_DIR)
    set(SEARCH_DIRS ${PROJECT_SOURCE_DIR}/../../git/vistrails ${PROJECT_SOURCE_DIR}/../vistrails )
    find_path(VISTRAILS_GIT_DIR scripts ${SEARCH_DIRS} NO_DEFAULT_PATH)
    mark_as_advanced(VISTRAILS_GIT_DIR)
endif (ALPS_PATCH_VISTRAILS AND NOT VISTRAILS_GIT_DIR)
