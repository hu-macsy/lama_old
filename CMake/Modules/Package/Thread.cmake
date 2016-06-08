###
 # @file Package/Thread.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief find Thread library
 # @author Lauretta Schubert
 # @date 20.08.2015
###

### SCAI_THREAD_LIBRARY - needed Thread library

enable_language ( C )

set ( CMAKE_THREAD_PREFER_PTHREAD 1 )
set ( THREADS_PREFER_PTHREAD_FLAG 1 )

find_package( Threads ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )

if     ( CMAKE_THREAD_LIBS_INIT )
    set ( SCAI_THREAD_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} CACHE PATH "Pthread library" )
elseif ( ${CMAKE_USE_PTHREADS_INIT} )
	set ( SCAI_THREAD_LIBRARIES pthread CACHE PATH "Pthread library" )
endif  ( )

mark_as_advanced( SCAI_THREAD_LIBRARIES )

if    ( APPLE )
	execute_process ( COMMAND /usr/bin/otool -L /usr/lib/libSystem.dylib OUTPUT_VARIABLE _pthread_output )
	string ( REGEX MATCH "current version ([0-9]+[.]*[0-9]*)" __pthread_output ${_pthread_output} )
	string ( REGEX MATCH "([0-9]+[.]*[0-9]*)" SCAI_THREAD_VERSION ${__pthread_output} )
elseif ( UNIX )
	## get pthread version
	execute_process ( COMMAND /usr/bin/getconf GNU_LIBPTHREAD_VERSION OUTPUT_VARIABLE _pthread_output )
	string ( REGEX MATCH "([0-9]+\\.[0-9]*)" SCAI_THREAD_VERSION ${_pthread_output} )
endif (  )

###  Here we use PThread library for threads
###  Note: FindThreads in CMake is available as Module, but is buggy, needs update of CheckIncludeFiles.cmake
#find_library ( PTHREADS_LIBRARY NAMES pthread pthreads )
#set ( SCAI_THREAD_LIBRARY ${PTHREADS_LIBRARY} )
