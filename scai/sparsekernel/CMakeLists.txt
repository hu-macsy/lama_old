###
 # @file scai/sparsekernel/CMakeLists.txt
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
 # @brief scai/sparsekernel/CMakeLists.txt
 # @author Eric Schricker
 # @date 18.02.2016
###

## Define name of the SCAI module project

scai_module(

    MODULE_NAME   sparsekernel
    INTERNAL_DEPS common logging tracing tasking hmemo kregistry blaskernel utilskernel 
    EXTERNAL_DEPS OpenMP CUDA MIC SCAI_BLAS 
)

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp


    HEADERS                  # .hpp only

        COOKernelTrait
        CSRKernelTrait
        DenseKernelTrait
        DIAKernelTrait
        ELLKernelTrait
        JDSKernelTrait
        StencilKernelTrait
        
)

## add subdirectories

scai_subdirectories( cuda examples external mic openmp test )

### add library ###

scai_library()

if ( USE_MIC )
    set ( MKLROOT $ENV{MKLROOT} )
    set ( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} 
          "-offload-option,mic,compiler,\"-L${MKLROOT}/lib/mic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group\"" )
    #ToDo: set_target_properties not working, any better way?
endif ( USE_MIC )

### install ###

## install headers
relative_install( FILES ${CXX_HEADERS} DESTINATION include/scai/${MODULE_NAME} )

## install library.hpp
install ( FILES ../${MODULE_NAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 

include ( CustomCommands/SphinxDoc )
