##############################################################################
#  Code coverage with gcov/lcov
##############################################################################

set ( LAMA_USE_CODE_COVERAGE FALSE CACHE BOOL "Enable / Disable use of Code Coverage" )

if ( LAMA_USE_CODE_COVERAGE )
    set ( COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage" )
# TODO: move to copmiler flags!
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}" )
endif ( LAMA_USE_CODE_COVERAGE )