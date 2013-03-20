### set include directories ### 

include_directories( . )


### set source files ###

set ( CXXSOURCES
        SourceLocation.cpp
        Level.cpp
        GenLogger.cpp
        GenLoggerCreator.cpp
        LoggerProvider.cpp
        Logger.cpp
    )


### add definitions ###

# GenLogger uses also logging
add_definitions ( -DLAMA_LOG_LEVEL_DEBUG )

if ( WIN32 )
	add_definitions ( -DLOG4LAMA_COMPILING_DLL )
endif ( WIN32 )


### Glob for headers to make them available in Visual Studio Projects ###
file ( GLOB_RECURSE HEADERS *.hpp )


### add library ###

add_library ( log4lama SHARED ${CXXSOURCES} ${HEADERS} )

### Installation: IncludeFiles + Library ###

file ( GLOB IncludeFiles *.hpp )

install ( FILES ${IncludeFiles} DESTINATION include/logging )
install ( TARGETS log4lama DESTINATION lib )
