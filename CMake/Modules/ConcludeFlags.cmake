#### concluding all defined compiler flags to CMAKE_..._FLAGS ####

## scai common adds OpenMP_CXX_FLAGS to SCAI_COMMON_FLAGS if found
if    ( SCAI_COMMON_FOUND ) 
	set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} ${SCAI_COMMON_FLAGS}" )
else  ( SCAI_COMMON_FOUND )
	set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif ( SCAI_COMMON_FOUND ) 

## add variables to cache with new names so they can be modified by the user via CCMAKE

set ( ADDITIONAL_CXX_FLAGS "${LAMA_CXX_FLAGS}" CACHE STRING "additional flags for cxx compile and link" )
set ( ADDITIONAL_WARNING_FLAGS "${SCAI_WARNING_FLAGS}" CACHE STRING "compilation flags concerning warnings" )
set ( ADDITIONAL_CXX_FLAGS_RELEASE "${LAMA_CXX_FLAGS_RELEASE}" CACHE STRING "addtional cxx compiler flags for release optimizations" )
set ( ADDITIONAL_LINKER_FLAGS "${LAMA_LINKER_FLAGS}" CACHE STRING "additional linker flags" )

mark_as_advanced ( ADDITIONAL_CXX_FLAGS ADDITIONAL_WARNING_FLAGS ADDITIONAL_CXX_FLAGS_RELEASE ADDITIONAL_LINKER_FLAGS )

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS}")
set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ADDITIONAL_CXX_FLAGS_RELEASE} " )
set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )
set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )