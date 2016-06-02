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
		# check for needed C++11 features
		# only use C++11 if all needed features are found

		include ( Functions/checkFeature )

		#set ( CXX11_FEATURE_LIST    "auto" "nullptr" "lambda" "static_assert" "rvalue_references" "decltype" "cstdint" "long_long" "variadic_templates" "constexpr" "sizeof_member" "__func__" )
		#set ( CXX11_FEATURE_NUMBER  "2546" "2431"    "2927"   "1720"          "2118"              "2343"     ""        "1811"      "2555"               "2235"      "2253"          "2340"     )

		set ( CXX11_FEATURE_LIST    "auto" )
		set ( CXX11_FEATURE_NUMBER  "2546" )

		set ( COUNT 0 )
		foreach    ( FEATURE ${CXX11_FEATURE_LIST} )
			list ( GET CXX11_FEATURE_NUMBER ${COUNT} FEATURE_NUMBER )
			checkFeature ( ${FEATURE} "${FEATURE_NUMBER}" ${FEATURE}_BOOLVALUE "-std=c++11" )
			message ( STATUS "${COUNT}: item ${FEATURE} number ${FEATURE_NUMBER}: ${${FEATURE}_BOOLVALUE}" )
			if    ( NOT ${FEATURE}_BOOLVALUE )
				set ( CXX_SUPPORTS_C11 FALSE )
			endif ( NOT ${FEATURE}_BOOLVALUE )
			math ( EXPR COUNT "${COUNT}+1" )
		endforeach ( FEATURE ${CXX11_FEATURE_LIST} )
	endif ( CXX_SUPPORTS_C11 )

else  ( NOT WIN32 )
	
	if    ( MSVC ) # if Visual Studio
		#message ( STATUS "MSVC_VERSION ${MSVC_VERSION} " )
		if    ( ${MSVC_VERSION} GREATER 1600 ) # 1600 = VS 10.0
			set ( CXX_SUPPORTS_C11 TRUE )
		else  ( ${MSVC_VERSION} GREATER 1600 )
			set ( CXX_SUPPORTS_C11 FALSE )
		endif ( ${MSVC_VERSION} GREATER 1600 )
	else  ( MSVC )
		message ( FATAL_ERROR "NO Visual Studio compiler, can not determine if C++11 support given." )
	endif ( MSVC )
	
endif ( NOT WIN32 )

set ( CXX_SUPPORTS_C11 ${CXX_SUPPORTS_C11} CACHE STRING "Compiler supports CXX-11." )

if    ( CXX_SUPPORTS_C11 AND NOT WIN32 ) # WIN32: implicitly compiled with c++11 by MSVC
    set ( SCAI_LANG_FLAGS "-std=c++11" )
else  ( CXX_SUPPORTS_C11 AND NOT WIN32 )
    set ( SCAI_LANG_FLAGS "" )
endif ( CXX_SUPPORTS_C11 AND NOT WIN32 )

set ( ADDITIONAL_CXX_FLAGS_LANG "${SCAI_LANG_FLAGS}" CACHE STRING "Language flag for using C++11 (if compiler capable)" )
mark_as_advanced ( ADDITIONAL_CXX_FLAGS_LANG )
