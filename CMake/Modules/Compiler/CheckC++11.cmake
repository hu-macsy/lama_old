###
 # @file CheckC++11.cmake
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

if ( NOT WIN32 )

    set ( CXX11_COMPILE_FLAG "-std=c++11" )

    # if ( CMAKE_CXX_COMPILER_ID MATCHES Clang )
    #     set ( CXX11_COMPILE_FLAG "${CXX11_COMPILE_FLAG} -stdlib=libc++" )
    # endif ( CMAKE_CXX_COMPILER_ID MATCHES Clang )

    include ( CheckCXXCompilerFlag )

    if ( NOT DEFINED CXX_SUPPORTS_C11 )
        CHECK_CXX_COMPILER_FLAG( ${CXX11_COMPILE_FLAG} CXX_SUPPORTS_C11 )
    endif ( NOT DEFINED CXX_SUPPORTS_C11 )

    if    ( CXX_SUPPORTS_C11 )
        # check for needed C++11 features
        # only use C++11 if all needed features are found
        include ( scai_function/checkFeature )

        # set c++11 flag as CMAKE_CXX_FLAGS and restore old value afterwards
        # setting flag int try_compile/try_run does not work
        set ( CHECK_CXX11_OLD_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} )
        set ( CMAKE_CXX_FLAGS "${CXX11_COMPILE_FLAG}" )

        # other c++(11) features - not used (yet)
        # add to list below if used in LAMA ( take care of ordering of the three lists (!) )
        #set ( CXX11_FEATURE_LIST    "auto" "nullptr" "lambda" "static_assert" "rvalue_references" "decltype" "variadic_templates" "constexpr" "sizeof_member" "__func__" )
        #set ( CXX11_FEATURE_NUMBER  "2546" "2431"    "2927"   "1720"          "2118"              "2343"     "2555"               "2235"      "2253"          "2340"     )
        #set ( CXX11_LINK_LIBRARIES  ""     ""        ""       ""              ""                  ""         ""                   ""          ""              ""         )

        # needed c++11 features - replacing boost                                            # replacing pthread
        set ( CXX11_FEATURE_LIST    "bind" "function" "shared_ptr" "unique_ptr" "weak_ptr" ) # "thread_local" )
        set ( CXX11_FEATURE_NUMBER  ""     ""         ""           ""           ""         ) # ""             )
        set ( CXX11_LINK_LIBRARIES  ""     ""         ""           ""           ""         ) # "pthread"      ) #todo:: pthreads only for gnu !?

        set ( COUNT 0 )
        foreach    ( FEATURE ${CXX11_FEATURE_LIST} )
            list ( GET CXX11_FEATURE_NUMBER ${COUNT} FEATURE_NUMBER )
            list ( GET CXX11_LINK_LIBRARIES ${COUNT} LINK_LIB )
            checkFeature ( ${FEATURE} "${FEATURE_NUMBER}" ${FEATURE}_BOOLVALUE "${CXX11_COMPILE_FLAG}" "${LINK_LIB}" )
            message ( "${COUNT}: item ${FEATURE} number ${FEATURE_NUMBER}: ${${FEATURE}_BOOLVALUE}" )

            if    ( NOT ${FEATURE}_BOOLVALUE )
                set ( CXX_SUPPORTS_C11 FALSE )
            endif ( NOT ${FEATURE}_BOOLVALUE )
            math ( EXPR COUNT "${COUNT}+1" )
        endforeach ( FEATURE ${CXX11_FEATURE_LIST} )

        if    ( NOT CXX_SUPPORTS_C11 )
            message ( STATUS "Compiler does not support all needed c++11 features. Turn CXX_SUPPORTS_C11 off.\n Unsupported features are: ${CXX11_UNSUPPORTED_FEATURE_LIST}" )
        endif ( NOT CXX_SUPPORTS_C11 )

        # restore CMAKE_CXX_FLAGS
        set ( CMAKE_CXX_FLAGS ${CHECK_CXX11_OLD_CMAKE_CXX_FLAGS} )

    else ( CXX_SUPPORTS_C11 )
        message( STATUS "Compiler does not support flag: ${CXX11_COMPILE_FLAG}" )
    endif ( CXX_SUPPORTS_C11 )

else ( NOT WIN32 )
    
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
    set ( SCAI_LANG_FLAGS "${CXX11_COMPILE_FLAG}" )
else  ( CXX_SUPPORTS_C11 AND NOT WIN32 )
    set ( SCAI_LANG_FLAGS "" )
endif ( CXX_SUPPORTS_C11 AND NOT WIN32 )

set ( ADDITIONAL_CXX_FLAGS_LANG "${SCAI_LANG_FLAGS}" CACHE STRING "Language flag for using C++11 (if compiler capable)" )

mark_as_advanced ( ADDITIONAL_CXX_FLAGS_LANG )
