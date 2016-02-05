## adding packages


macro    ( addInternalPackages )
	message ( STATUS "add internal Packages" )
	foreach    ( PACKAGE_TO_FIND ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
    	find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
    	message ( STATUS "Package: ${PACKAGE_TO_FIND}")
	endforeach ( PACKAGE_TO_FIND ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
endmacro ( addInternalPackages )

macro    ( addExternalPackages )
	message ( STATUS "add external Packages" )
	foreach    ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
    	include( Package/${module} )
    	message ( STATUS "Package: ${module}")
	endforeach ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
endmacro ( addExternalPackages )

macro    ( addInternalAndExternalPackages )
	addInternalPackages()
	addExternalPackages()
endmacro ( addInternalAndExternalPackages )

## adding includes dirs from packages

macro    ( addInternalIncludes )
	foreach    ( module ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
		string ( TOUPPER ${module} upper_module )
    	include_directories( ${${upper_module}_INCLUDE_DIR} )
	endforeach ( module ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
endmacro ( addInternalIncludes )

macro    ( addExternalIncludes )
	foreach    ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
		string ( TOUPPER ${module} upper_module )
    	include_directories( ${${upper_module}_INCLUDE_DIR} )
	endforeach ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
endmacro ( addExternalIncludes )

macro    ( addInternalAndExternalIncludes )
	foreach    ( module ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
		string ( TOUPPER ${module} upper_module )
    	include_directories( ${${upper_module}_INCLUDE_DIR} )
	endforeach ( module ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
endmacro ( addInternalAndExternalIncludes )

## adding link libraries of packages

macro    ( addInternalLinkLibraries )
	set ( REVERT_LIST ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} ) # because list does not accept variable recursion
	list ( REVERSE REVERT_LIST )
	foreach    ( module ${REVERT_LIST} )
		string ( TOUPPER ${module} upper_module )
		set ( ${UPPER_PROJECT_NAME}_LINK_LIBRARIES ${${UPPER_PROJECT_NAME}_LINK_LIBRARIES} ${${upper_module}_LIBRARY} )
	endforeach ( module ${REVERT_LIST} )
	target_link_libraries ( ${PROJECT_NAME} ${SCAI_START_LINK_LIBRARIES} ${${UPPER_PROJECT_NAME}_LINK_LIBRARIES} ${SCAI_END_LINK_LIBRARIES} )
endmacro ( addInternalLinkLibraries )

macro    ( addExternalLinkLibraries )
	foreach    ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
		string ( TOUPPER ${module} upper_module )
    	target_link_libraries ( ${PROJECT_NAME} ${SCAI_${upper_module}_LIBRARIES} )
	endforeach ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
endmacro ( addExternalLinkLibraries )

macro    ( addInternalAndExternalLinkLibraries )
	addInternalLinkLibraries()
	addExternalLinkLibraries()
endmacro ( addInternalAndExternalLinkLibraries )
