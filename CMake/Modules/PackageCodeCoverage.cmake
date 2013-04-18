##############################################################################
#  Code coverage with gcov/lcov
##############################################################################

if ( NOT CODE_COVERAGE )
    set ( DEFAULT_CODE_COVERAGE FALSE )
endif ( NOT CODE_COVERAGE )

set ( CODE_COVERAGE ${DEFAULT_CODE_COVERAGE} CACHE BOOL "Enable / Disable Code Coverage" )

if ( CODE_COVERAGE )
    set ( COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage" )
# TODO: move to copmiler flags!
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}" )
endif ( CODE_COVERAGE )