###
 # @file scai/benchmark/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief CMake configuration file for SCAI package benchmark
 # @author Thomas Brandes
 # @date 13.09.2017
###

## Include required scai macros

include ( scai_macro/scai_module )
include ( scai_macro/scai_project )
include ( scai_macro/scai_subdirectories )
include ( scai_macro/scai_library )
include ( scai_macro/scai_test_scripts )
include ( scai_function/relative_install )
 
## Benchmark framework uses only common classes and logging
 
scai_module (

    MODULE_NAME   benchmark
    INTERNAL_DEPS common logging
) 

scai_project (

    CLASSES                  # .cpp, .hpp

        Benchmark
        BenchmarkPrinter
        BFError
        BFException
        Config
        InputSet
    )

## Dummy library required to resolve dependencies

scai_library ( PREFIX  ${SCAI_LIBRARY_PREFIX}
               TYPE    ${SCAI_LIBRARY_TYPE}
               VERSION ${SCAI_VERSION}
               ${CXX_SOURCES}
              )

### add executables ###

foreach ( executable RunBenchmark RegisteredInputSets RegisteredBenchmarks )

    add_executable ( ${executable} ${executable}.cpp )

    target_link_libraries ( ${executable} ${MODULE_LIBRARY} )

    install( TARGETS ${executable} DESTINATION bin )

endforeach ( executable )

install ( PROGRAMS BenchmarkRunner DESTINATION bin )

install( DIRECTORY pybench DESTINATION bin )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 

include ( CustomCommands/SphinxDoc )

