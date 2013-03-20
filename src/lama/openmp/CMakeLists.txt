### set include directories ### 

include_directories ( . )


### set source files ###

set ( OPENMP_SOURCES

      openmp/OpenMPInterface.cpp
      
      openmp/BLASHelper.cpp

      openmp/OpenMPBLAS1.cpp
      openmp/OpenMPBLAS2.cpp
      openmp/OpenMPBLAS3.cpp
      openmp/OpenMPLAPACK.cpp

      openmp/OpenMPUtils.cpp
      openmp/OpenMPCSRUtils.cpp
      openmp/OpenMPELLUtils.cpp
      openmp/OpenMPCOOUtils.cpp
      openmp/OpenMPDenseUtils.cpp
      openmp/OpenMPJDSUtils.cpp
      openmp/OpenMPDIAUtils.cpp
    )

if ( SCALAPACK_FOUND AND MPI_FOUND )

# can only be compiled if MKL ScaLAPACK is availabe

set ( OPENMP_SOURCES ${OPENMP_SOURCES}
      openmp/OpenMPSCALAPACK.cpp
    )

endif()

set ( CXXSOURCES 
      ${CXXSOURCES} 
      ${OPENMP_SOURCES} 

      PARENT_SCOPE
    )

### install ###

file ( GLOB includefiles *.hpp )
install ( FILES ${includefiles} DESTINATION include/lama/openmp )
