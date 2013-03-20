cmake_minimum_required (VERSION 2.8)

project (FastGramSchmidt)

# make sure that the default is a RELEASE

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Package OpenMP works only for C

add_library(QR QR.cpp)

set ( OPENMP_FLAG "-fopenmp" )
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_FLAG}" )

add_executable(main main.cpp)
add_executable(bench bench.cpp)


set ( MKL_ROOT_PATH /home/brandes/intel/mkl )

include_directories( ${MKL_ROOT_PATH}/include )

target_link_libraries( main 
                       -L${MKL_ROOT_PATH}/lib/intel64
                       -lmkl_intel_lp64
                       -lmkl_sequential
                       -lmkl_core
                      )

target_link_libraries( bench
                       -L${MKL_ROOT_PATH}/lib/intel64
                       -lmkl_intel_lp64
                       -lmkl_sequential
                       -lmkl_core
                      )

target_link_libraries( main QR )
target_link_libraries( bench QR )
