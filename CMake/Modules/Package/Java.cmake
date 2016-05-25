###
 # @file Java.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief findPackage and configuration of Java
 # @author Jan Ecker
 # @date 25.04.2013
###

if    ( NOT APPLE )
    if    ( CMAKE_VERSION VERSION_GREATER 2.8.11 )
        # Enable Java Compilation
        # this will set at least CMAKE_Java_COMPILER CMAKE_Java_ARCHIVE 
        find_package( Java COMPONENTS Development )

        setAndCheckCache ( JAVA )
        set ( USE_JAVA ${USE_JAVA} CACHE BOOL "Enable / Disable use of JAVA" )

        # LAMA irrelevant entries will be removed from cmake GUI completely
        set ( Java_JAR_EXECUTABLE "${Java_JAR_EXECUTABLE}" CACHE INTERNAL "" )
        set ( Java_JAVADOC_EXECUTABLE "${Java_JAVADOC_EXECUTABLE}" CACHE INTERNAL "" )
        set ( Java_JAVAH_EXECUTABLE "${Java_JAVAH_EXECUTABLE}" CACHE INTERNAL "" )
        set ( Java_JAVA_EXECUTABLE "${Java_JAVA_EXECUTABLE}" CACHE INTERNAL "" )
    endif ( CMAKE_VERSION VERSION_GREATER 2.8.11 )
endif ( NOT APPLE )
