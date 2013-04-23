##############################################################################
#  Code coverage with gcov/lcov
##############################################################################

setAndCheckCache ( CODE_COVERAGE "Code Coverage" )

if ( CODE_COVERAGE_FOUND AND LAMA_USE_CODE_COVERAGE )
    set ( COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage" )
# TODO: move to copmiler flags!
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}" )
endif ( CODE_COVERAGE_FOUND AND LAMA_USE_CODE_COVERAGE )