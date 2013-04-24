### C Compiler
if(CMAKE_COMPILER_IS_GNUCC)
    execute_process( COMMAND ${CMAKE_C_COMPILER} --version OUTPUT_VARIABLE _compiler_output )
    string( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" GNUCC_COMPILER_VERSION ${_compiler_output})
endif (CMAKE_COMPILER_IS_GNUCC )

### CXX Compiler
if(CMAKE_COMPILER_IS_GNUCXX)
    execute_process( COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE _compiler_output )
    string( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" GNUCXX_COMPILER_VERSION ${_compiler_output})
endif (CMAKE_COMPILER_IS_GNUCXX )