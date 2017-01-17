###
 # @file Dependencies/NeededCacheVariables.cmake
 #
 # @license
 # Copyright (c) 2009-2016
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
 # @brief Whitelist ( Cache Variables that should be passed to subprojects )
 # @author Lauretta Schubert
 # @date 17.08.2015
###

set ( WHITELIST
        # FLAGS
        ADDITIONAL_CXX_FLAGS
        ADDITIONAL_CXX_FLAGS_CODE_COVERAGE
        ADDITIONAL_CXX_FLAGS_DEBUG
        ADDITIONAL_CXX_FLAGS_NO_OFFLOAD
        ADDITIONAL_CXX_FLAGS_OPENMP
        ADDITIONAL_CXX_FLAGS_RELEASE
        ADDITIONAL_CXX_FLAGS_STATIC
        ADDITIONAL_LINKER_FLAGS
        ADDITIONAL_WARNING_FLAGS
        # EXAMPLES
        BUILD_EXAMPLES
        # TEST
        BUILD_TEST
        # DOC
        BUILD_DOC
        SCAI_DOC_TYPE
        Sphinx-apidoc_EXECUTABLE
        Sphinx-build_EXECUTABLE
        Sphinx_DIR
        PYTHON_EXECUTABLE
        # C++11
        CXX_SUPPORTS_C11
        CXX11_SUPPORTED_FEATURE_LIST
        CXX11_UNSUPPORTED_FEATURE_LIST
        # CMAKE
        CMAKE_BUILD_TYPE
        CMAKE_INSTALL_PREFIX
        CMAKE_CXX_COMPILER
        CMAKE_C_COMPILER # may be unused
        CMAKE_SKIP_INSTALL_RPATH
        CMAKE_SKIP_RPATH
        # SCAI
        SCAI_ASSERT_LEVEL
        SCAI_CMAKE_VERBOSE
        SCAI_HOST_TYPES_LIST
        SCAI_INDEX_TYPE        
        SCAI_LIBRARY_TYPE
        SCAI_LOGGING_LEVEL
        SCAI_TRACING
        # USE
        USE_CODE_COVERAGE
        USE_OPENMP             # put this here, because OpenMP is used with 'common' in each subproject
        USE_JAVA               # make it availabe in sub project
    )
    
set ( BOOST_WHITELIST
        BOOST_ROOT
        Boost_INCLUDE_DIR
        Boost_LIBRARY_DIRS
        Boost_NO_BOOST_CMAKE
        Boost_UNIT_TEST_FRAMEWORK_LIBRARY
        Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG
        Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE
    )

set ( CUDA_WHITELIST
        ADDITIONAL_NVCC_FLAGS
        ADDITIONAL_NVCC_FLAGS_DEBUG
        ADDITIONAL_NVCC_RELEASE_FLAGS
        CUDA_64_BIT_DEVICE_CODE
        #CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE
        #CUDA_BUILD_CUBIN
        #CUDA_BUILD_EMULATION
        CUDA_COMPUTE_CAPABILITY
        CUDA_CUDART_LIBRARY
        CUDA_CUDA_LIBRARY
        CUDA_GENERATE_CODE
        CUDA_HAVE_GPU
        CUDA_HOST_COMPILATION_CPP
        CUDA_HOST_COMPILER
        CUDA_NVCC_EXECUTABLE
        #CUDA_NVCC_FLAGS
        #CUDA_NVCC_FLAGS_DEBUG
        #CUDA_NVCC_FLAGS_RELEASE
        CUDA_PROPAGATE_HOST_FLAGS
        CUDA_SEPARABLE_COMPILATION
        #CUDA_TARGET_CPU_ARCH
        CUDA_TOOLKIT_INCLUDE
        CUDA_TOOLKIT_ROOT_DIR
        #CUDA_TOOLKIT_TARGET_DIR
        CUDA_VERBOSE_BUILD
        CUDA_VERSION
        CUDA_cublas_LIBRARY
        #CUDA_cufft_LIBRARY
        #CUDA_cupti_LIBRARY
        #CUDA_curand_LIBRARY
        CUDA_cusparse_LIBRARY
        #CUDA_nppc_LIBRARY
        #CUDA_nppi_LIBRARY
        #CUDA_npps_LIBRARY
        USE_CUDA
        USE_CUSPARSE
    )

set ( GPI_WHITELIST
        USE_GPI
        GPI2_ROOT
        IBVERBS_ROOT
    )

set ( GRAPHPARTITIONING_WHITELIST
        GRAPHPARTITIONING_INCLUDE_DIR
        METIS_INCLUDE_DIR
        METIS_LIBRARY
        METIS_ROOT
        PARMETIS_INCLUDE_DIR
        PARMETIS_LIBRARY
        PARMETIS_ROOT
        SCAI_GRAPHPARTITIONING_LIBRARIES
        USE_GRAPHPARTIONING
    )
    
set ( MIC_WHITELIST 
        USE_MIC
    )

set ( MPI_WHITELIST 
        MPIEXEC
        MPIEXEC_MAX_NUMPROCS
        MPIEXEC_NUMPROC_FLAG
        MPI_CXX_COMPILER
        MPI_CXX_INCLUDE_PATH
        MPI_CXX_LIBRARIES
        MPI_CXX_LINK_FLAGS
        MPI_C_COMPILER
        MPI_C_INCLUDE_PATH
        MPI_C_LIBRARIES
        MPI_C_LINK_FLAGS
        MPI_EXTRA_LIBRARY
        MPI_LIBRARY
        USE_MPI
    )

set ( OPENMP_WHITELIST
        OpenMP_CXX_FLAGS
    )
    
set ( SCAI_BLAS_WHITELIST
        BLAS_blas_LIBRARY
        BLAS_goto2_LIBRARY
        #GFORTRAN_LIBRARY
        MKL_Is64
        MKL_LIBRARY_CORE
        MKL_LIBRARY_GNU
        MKL_LIBRARY_INTEL
        MKL_LIBRARY_LP64
        LAPACK_goto2_LIBRARY
        LAPACK_lapack_LIBRARY
        SCAI_BLAS_LIBRARY
        USE_BLAS
        USE_LAPACK
        USE_MKL
        USE_SCALAPACK
    )

set ( JAVA_WHITELIST
		Java_JAVAC_EXECUTABLE
	)
    
#set ( THREAD_WHITELIST
#    )
