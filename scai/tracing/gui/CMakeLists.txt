
# Enable Java Compilation
# this will set at least CMAKE_Java_COMPILER CMAKE_Java_ARCHIVE 

find_package( Java )

# Note: Java_FOUND not set by some versions, so ask for JAVAC executable

if ( Java_JAVAC_EXECUTABLE )

	include( UseJava )
	
	# directory where to find the needed jar files of other Java code
	
	set( LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib )
	
	# Java GUI uses Java Graph Package Grappa
	find_jar( JAR1 grappa1_2 PATHS ${LIBDIR} )
	
	find_jar( JAR2 jfreechart-1.0.8 PATHS ${LIBDIR} )
	find_jar( JAR3 jcommon-1.0.12 PATHS ${LIBDIR} )
	find_jar( JAR4 log4j-1.2.13 PATHS ${LIBDIR} )
	
	mark_as_advanced( JAR1 JAR2 JAR3 JAR4 )
	
	file( GLOB_RECURSE JAVA_SRC "*.java" )
	
	set(VERSION 1.0.0)
	
	set(GUI_JAR "trace_gui-${VERSION}")
	
	add_jar( ${GUI_JAR} 
	         SOURCES ${JAVA_SRC}
	         INCLUDE_JARS ${JAR1} ${JAR2} ${JAR3} ${JAR4} )
	
	install_jar( ${GUI_JAR} bin )
	
	install( FILES ${JAR1} DESTINATION bin )
	install( FILES ${JAR2} DESTINATION bin )
	install( FILES ${JAR3} DESTINATION bin )
	install( FILES ${JAR4} DESTINATION bin )
	
	configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/bin/TraceCT.in" "${PROJECT_BINARY_DIR}/TraceCT" )
	
	install( FILES ${PROJECT_BINARY_DIR}/TraceCT 
	         PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ
	         DESTINATION bin )

endif ( Java_JAVAC_EXECUTABLE )
