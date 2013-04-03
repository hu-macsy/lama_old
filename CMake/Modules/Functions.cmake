# Function to set FOUND_XXX variables depending on the corresponding LAMA_USE_XXX variable
function ( setAndCheckCache PACKAGE_NAME )
    # Create variable names with LAMA_USE_XXX and FOUND_XXX
    set ( CACHE_VARIABLE_NAME LAMA_USE_${PACKAGE_NAME} )
    set ( FOUND_VARIABLE_NAME ${PACKAGE_NAME}_FOUND )

    # Check if cache variable is already set
    if ( DEFINED ${CACHE_VARIABLE_NAME} )
        # if use of package is enabled
        if ( ${CACHE_VARIABLE_NAME} )
            if ( NOT ${FOUND_VARIABLE_NAME} )
                # if package is enabled, but not found: ERROR!
                message ( FATAL_ERROR "${PACKAGE_NAME} enabled, but not found!" )
            endif ( NOT ${FOUND_VARIABLE_NAME} )
        
        # if use of package is disabled
        else ( ${CACHE_VARIABLE_NAME} )
            # disable the package
            set ( ${FOUND_VARIABLE_NAME} FALSE PARENT_SCOPE )
        endif ( ${CACHE_VARIABLE_NAME} )
    
    # if cache variable is NOT set
    else ( DEFINED ${CACHE_VARIABLE_NAME} )
        # Check if package was found
        if ( ${FOUND_VARIABLE_NAME} )
            set ( USE_PACKAGE TRUE )
        else ( ${FOUND_VARIABLE_NAME} )
            set ( USE_PACKAGE FALSE )
        endif ( ${FOUND_VARIABLE_NAME} )
        
        # Set cache variable
        set ( ${CACHE_VARIABLE_NAME} ${USE_PACKAGE} CACHE BOOL "Enable / Disable use of ${PACKAGE_NAME}" )
    endif ( DEFINED ${CACHE_VARIABLE_NAME} )
endfunction ( setAndCheckCache )


function ( get_relative_path RELATIVE_PATH )
    # get relative path
    string ( LENGTH "${CMAKE_SOURCE_DIR}lama/" CMAKE_SOURCE_DIR_LENGTH )
    string ( LENGTH ${CMAKE_CURRENT_SOURCE_DIR} CMAKE_CURRENT_SOURCE_DIR_LENGTH )
    math ( EXPR PATH_LENGTH ${CMAKE_CURRENT_SOURCE_DIR_LENGTH}-${CMAKE_SOURCE_DIR_LENGTH} )
    string ( SUBSTRING ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR_LENGTH} ${PATH_LENGTH} PATH )
    
    # set return parameter (via PARENT_SCOPE)
    set ( ${RELATIVE_PATH} "${PATH}/" PARENT_SCOPE )
endfunction ( get_relative_path )

function ( lama_get_relative_path RELATIVE_PATH )
    # get relative path
    string ( LENGTH "${LAMA_SOURCE_DIR}" LAMA_SOURCE_DIR_LENGTH )
    string ( LENGTH ${CMAKE_CURRENT_SOURCE_DIR} CMAKE_CURRENT_SOURCE_DIR_LENGTH )
   
    if ( ${LAMA_SOURCE_DIR_LENGTH} LESS ${CMAKE_CURRENT_SOURCE_DIR_LENGTH} )
        math ( EXPR LAMA_SOURCE_DIR_LENGTH ${LAMA_SOURCE_DIR_LENGTH}+1 )
        set ( PATH_SUFFIX / )
    endif ()

    math ( EXPR PATH_LENGTH ${CMAKE_CURRENT_SOURCE_DIR_LENGTH}-${LAMA_SOURCE_DIR_LENGTH} )
    string ( SUBSTRING ${CMAKE_CURRENT_SOURCE_DIR} ${LAMA_SOURCE_DIR_LENGTH} ${PATH_LENGTH} PATH )
    
    # set return parameter (via PARENT_SCOPE)
    set ( ${RELATIVE_PATH} ${PATH}${PATH_SUFFIX} PARENT_SCOPE )
endfunction ( lama_get_relative_path )






function ( install_header_files )
    # get actual subdir
    string ( LENGTH ${CMAKE_SOURCE_DIR} CMAKE_SOURCE_DIR_LENGTH )
    string ( LENGTH ${CMAKE_CURRENT_SOURCE_DIR} CMAKE_CURRENT_SOURCE_DIR_LENGTH )
    math ( EXPR PATH_LENGTH ${CMAKE_CURRENT_SOURCE_DIR_LENGTH}-${CMAKE_SOURCE_DIR_LENGTH} )
    string ( SUBSTRING ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR_LENGTH} ${PATH_LENGTH} SUB_PATH )
    
    if ( DEFINED ARGV0 )
        set ( SUB_PATH ${SUB_PATH}/${ARGV0} )
        file ( GLOB INCLUDE_FILES ${ARGV0}/*.hpp )
    else ()
        file ( GLOB INCLUDE_FILES *.hpp )
    endif ( DEFINED ARGV0 )
    
    # find all *.hpp files and copy them in the correct subdir
    install ( FILES ${INCLUDE_FILES} DESTINATION "include${SUB_PATH}" )
endfunction ( install_header_files )




macro ( lama_set_source_dir )
    set ( LAMA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
    set ( CXX_SOURCES "" )
endmacro ( lama_set_source_dir )


# Needs to be macros not functions, because modifications of the parent scope
macro ( lama_classes )
    lama_sources ( ${ARGN} )
    lama_headers ( ${ARGN} )
endmacro ( lama_classes )

macro ( lama_sources )    
    lama_get_relative_path ( LAMA_RELATIVE_PATH )
    
    foreach(SOURCE_FILE ${ARGN})
        set ( CXX_SOURCES ${CXX_SOURCES} "${LAMA_RELATIVE_PATH}${SOURCE_FILE}.cpp" )
    endforeach()
endmacro ( lama_sources )

macro ( lama_headers )
    lama_get_relative_path ( LAMA_RELATIVE_PATH )
    
    # clear CXX_HEADERS
    set ( CXX_HEADERS "" )
    
    foreach(SOURCE_FILE ${ARGN})
        set ( CXX_SOURCES ${CXX_SOURCES} "${LAMA_RELATIVE_PATH}${SOURCE_FILE}.hpp" )
        set ( CXX_HEADERS ${CXX_HEADERS} "${SOURCE_FILE}.hpp" )
    endforeach()
    
    # install CXX_HEADERS
    install ( FILES ${CXX_HEADERS} DESTINATION "include/lama/${LAMA_RELATIVE_PATH}" )
endmacro ( lama_headers )

macro ( lama_add )
    # publish CXX_SOURCES in parent scope
    set ( CXX_SOURCES ${CXX_SOURCES} PARENT_SCOPE )
endmacro ( lama_add )















