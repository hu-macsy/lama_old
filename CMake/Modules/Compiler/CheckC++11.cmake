###
 # @file CheckC++11.cmake
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
 # @brief Check if compiler supports C++11 features.
 # @author Thomas Brandes
 # @date 09.06.2015
###

#### Check for -std=c++11
#### if supported it is set to SCAI_LANG_FLAGS
#
#
#  If C++ compiler supports ISO C++11 standard we will use smart pointers (shared_ptr, unique_ptr)
#  and function/bind. If not, features of Boost libraries are employed. 
# 
#  Note: In C++ files the following macro can be used to see if the flag -std=c++11 has been specified:
#  if __cplusplus > 199711L  
#    <code using C++11 features>
#  else
#    <code using boost>
#  endif

# use CMake module that provides CHECK_CXX_COMPILER_FLAG 

if    ( NOT WIN32 )

include ( CheckCXXCompilerFlag )

if ( NOT DEFINED CXX_SUPPORTS_C11 )
    CHECK_CXX_COMPILER_FLAG( -std=c++11 CXX_SUPPORTS_C11 )
endif ( NOT DEFINED CXX_SUPPORTS_C11 )

if    ( CXX_SUPPORTS_C11 )
    set ( SCAI_LANG_FLAGS "-std=c++11" )
else  ( CXX_SUPPORTS_C11 )
    set ( SCAI_LANG_FLAGS "" )
endif ( CXX_SUPPORTS_C11 )

set ( CXX_SUPPORTS_C11 ${CXX_SUPPORTS_C11} CACHE STRING "Compiler supports CXX-11." )

else  ( NOT WIN32 )
	
	if    ( MSVC ) # if Visual Studio
		#message ( STATUS "MSVC_VERSION ${MSVC_VERSION} " )
		if    ( ${MSVC_VERSION} GREATER 1600 ) # 1600 = VS 10.0
			set ( CXX_SUPPORTS_C11 TRUE )
			set ( SCAI_LANG_FLAGS "" ) # implicitly compiled with c++11
		else  ( ${MSVC_VERSION} GREATER 1600 )
			set ( CXX_SUPPORTS_C11 FALSE )
			set ( SCAI_LANG_FLAGS "" )
		endif ( ${MSVC_VERSION} GREATER 1600 )
	else  ( MSVC )
		message ( FATAL_ERROR "NO Visual Studio compiler, can not determine if C++11 support given." )
	endif ( MSVC )
	
endif ( NOT WIN32 )

set ( ADDITIONAL_CXX_FLAGS_LANG "${SCAI_LANG_FLAGS}" CACHE STRING "Language flag for using C++11 (if compiler capable)" )
mark_as_advanced ( ADDITIONAL_CXX_FLAGS_LANG )
