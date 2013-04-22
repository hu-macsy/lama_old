# Find required packages
set ( REQUIRED_PACKAGES_TO_FIND
        OpenMP
        LAMA_BLAS
        #add required packages here
    )
    
# Find optional packages
set ( OPTIONAL_PACKAGES_TO_FIND
        Doxygen
        Threads
        #add optional packages here
    )
    
#CUDA Only works with GCC on Linux
#TODO: This needs to be checked on windows
if ( CMAKE_COMPILER_IS_GNUCC )
    set ( OPTIONAL_PACKAGES_TO_FIND
          ${OPTIONAL_PACKAGES_TO_FIND}
          CUDA
    )
endif ( CMAKE_COMPILER_IS_GNUCC )