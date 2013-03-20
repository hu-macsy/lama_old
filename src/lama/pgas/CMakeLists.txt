# This directory contains the optional LAMA-PGAS interface
# It builds as target the shared library 'lamapgas'

include_directories ( . )

### set source files ###

set ( PGASSOURCES
        LAMAInstance.cpp
        LAMAArguments.cpp
        PGASInterface.cpp
        PGASCommunicator.cpp
        PGASCommunicatorManager.cpp
        PGASSyncToken.cpp
        PGASContext.cpp
        PGASFunctor.cpp
        PGASContextManager.cpp
        PGASNoInterface.cpp
        PGASCheckInterface.cpp
        PGASAllocator.cpp
        PGASSimpleAllocator.cpp
    )
     
#TODO change to PGASFound !!!!!!!!!!!

if ( GASNET_FOUND )
    if ( OSHMEM_FOUND )
        message( STATUS "Including OpenSHMEMInterface" )
        set ( PGASSOURCES ${PGASSOURCES} pgas/OpenShMemInterface.cpp )
    endif ( OSHMEM_FOUND )
endif ( GASNET_FOUND )

if ( GPI_FOUND )
    message( STATUS "Including GPIInterface" )
    set ( PGASSOURCES
          ${PGASSOURCES} 
          pgas/GPIInterface.cpp
    )
endif ( GPI_FOUND )

set ( PGASSOURCES
      ${PGASSOURCES} 
      PARENT_SCOPE
)
if ( GPI_FOUND )
    include_directories ( ${GPI_INCLUDE_DIR} )
endif ( GPI_FOUND )

if ( GASNET_FOUND AND OSHMEM_FOUND )
    include_directories ( ${OSHMEM_INCLUDE_DIR} )
endif (  GASNET_FOUND AND OSHMEM_FOUND )

### add library ###

add_library ( lamapgas  SHARED ${PGASSOURCES} )

### set target link libraries ###  

#TODO change to PGASFound !!!!!!!!!!!
if ( OSHMEM_FOUND )
    if ( GASNET_FOUND )
        target_link_libraries ( lamapgas ${OSHMEM_LIBRARIES} )
        target_link_libraries ( lamapgas ${GASNET_LIBRARIES} )
        target_link_libraries ( lamapgas ibverbs )
    else ( GASNET_FOUND )
        message ( STATUS "OpenShMem without GASNET cannot be used" )
    endif ( GASNET_FOUND )
endif ( OSHMEM_FOUND )
if ( GPI_FOUND )
    target_link_libraries ( lamapgas ${GPI_LIBRARIES} )
    target_link_libraries ( lamapgas ibverbs )
endif ( GPI_FOUND )

# Install of the LAMA-PGAS interface library

install ( TARGETS lamapgas DESTINATION lib )

