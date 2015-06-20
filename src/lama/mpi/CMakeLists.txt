###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief lama/mpi/CMakeLists.txt
 # @author Jan Ecker
 # @date 02.04.2013
 # @since 1.0.0
###

### set this dir as lama source dir ###
lama_set_source_dir ()

### Add classes and header files ###
set ( CLASSES
        MPICommunicator
        MPICommunicatorManager
        MPIException
        MPISyncToken
    )
    
set ( HEADERS
        MPIUtils
    )

lama_classes ( ${CLASSES} )
lama_headers ( ${HEADERS} )

### MPI source files need path to find "mpi.h"

include_directories ( ${MPI_INCLUDE_PATH} )

################################################################################
#TODO: RECHECK!
if ( CMAKE_COMPILER_IS_GNUCXX )
    set_source_files_properties ( ${CXX_SOURCES} PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter" )
endif ( CMAKE_COMPILER_IS_GNUCXX )

################################################################################

### add library ###
add_library ( amampi ${LAMA_LIBRARY_TYPE} ${CXX_SOURCES} )

### set target link libraries ###  

target_link_libraries ( amampi ${MPI_LIBRARIES} )

### add definitions ###
if ( WIN32 )
    add_definitions ( -DLAMA_COMPILING_DLL )
endif ( WIN32 )

### install ####
install ( TARGETS amampi DESTINATION lib )
