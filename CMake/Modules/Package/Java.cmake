###
 # @file Java.cmake
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief findPackage and configuration of Java
 # @author Jan Ecker
 # @date 25.04.2013
###

##  Apple not supported here

if ( APPLE )
    return ()
endif ()

if ( CMAKE_VERSION VERSION_GREATER 2.8.11 )

    # Enable Java Compilation
    #  Java_JAVA_EXECUTABLE    = the full path to the Java runtime
    #  Java_JAVAC_EXECUTABLE   = the full path to the Java compiler
    #  Java_VERSION_STRING     = Version of the package found (java version), eg. 1.6.0_12
    #  Java_FOUND              - TRUE if all components are found, not helpful here

    find_package( Java COMPONENTS Development )

    if ( Java_JAVAC_EXECUTABLE )
        set ( JAVA_FOUND True )
    else ()
        set ( JAVA_FOUND False )
    endif ()

    if ( Java_JAR_EXECUTABLE )
        # that is okay, nothing to do
    else ()
        set ( JAVA_FOUND False )
        message( STATUS "JAR exectuable not found" )
    endif ()

    scai_build_variable ( NAME      USE_JAVA
                          BOOL 
                          DEFAULT   ${JAVA_FOUND}
                          DOCSTRING "use java for compiling tracing GUI" )

    # LAMA irrelevant entries will be removed from cmake GUI completely

    set ( Java_JAVADOC_EXECUTABLE "${Java_JAVADOC_EXECUTABLE}" CACHE INTERNAL "" )
    set ( Java_JAVAH_EXECUTABLE "${Java_JAVAH_EXECUTABLE}" CACHE INTERNAL "" )
    set ( Java_JAVA_EXECUTABLE "${Java_JAVA_EXECUTABLE}" CACHE INTERNAL "" )

    scai_summary_external ( NAME        Java
                            ENABLED     ${USE_JAVA}
                            FOUND       ${JAVA_FOUND} 
                            VERSION     ${Java_VERSION_STRING}
                            EXECUTABLE  "${Java_JAVAC_EXECUTABLE}, ${Java_JAR_EXECUTABLE}" )

endif ()
