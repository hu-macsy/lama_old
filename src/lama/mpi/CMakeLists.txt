# This directory contains the (optional) LAMA-MPI interface

### set include directories ###

include_directories ( . )

# This directory contains all relevant MPI stuff and is completely optional

set ( MPISOURCES
       MPICommunicator.cpp
       MPICommunicatorManager.cpp
       MPISyncToken.cpp
       MPIException.cpp
     )

if ( CMAKE_COMPILER_IS_GNUCXX )
    set_source_files_properties ( ${MPISOURCES} PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter" )
endif ( CMAKE_COMPILER_IS_GNUCXX )

if ( LAMA_TRACE_LEVEL STREQUAL "VT" )
    find_package ( VampirTrace REQUIRED)
    include_directories ( ${VAMPIRTRACE_INCLUDE_DIR} )
endif ( LAMA_TRACE_LEVEL STREQUAL "VT" )

add_library ( lamampi SHARED ${MPISOURCES} )

### set target link libraries ###  

if ( LAMA_TRACE_LEVEL STREQUAL "VT" )
    target_link_libraries ( lamampi ${VAMPIRTRACE_LIBRARIES} )
endif ( LAMA_TRACE_LEVEL STREQUAL "VT" )

target_link_libraries ( lamampi ${MPI_LIBRARIES} )
    
if ( NOT HAVE_INTEL_MPI )
    if ( METIS_FOUND )
        set ( PARTION_LIBRARIES ${PARTION_LIBRARIES} ${METIS_LIBRARIES} )
    endif ( METIS_FOUND )
    
    if ( PARMETIS_FOUND )
        set ( PARTION_LIBRARIES ${PARTION_LIBRARIES} ${PARMETIS_LIBRARIES} )
    endif ( PARMETIS_FOUND )

    target_link_libraries ( lamampi ${PARTION_LIBRARIES} )
endif ( NOT HAVE_INTEL_MPI ) 

### add definitions ###

if ( WIN32 )
    add_definitions ( -DLAMA_COMPILING_DLL )
endif ( WIN32 )

# Install of the LAMA-MPI library

install ( TARGETS lamampi DESTINATION lib )
