###
 # @file FindFFTW.cmake
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
 # @brief Find FFTW
 # @author Eric Schricker
 # @date 04.10.2016
###

#
# Find the FFTW includes and libraries
#
# FFTW_FOUND       - Do not attempt to use if "no" or undefined
# FFTW_INCLUDE_DIR - the FFTW include dir
# FFTW_LIBRARIES   - List of fully qualified libraries to link against
	
IF ( NOT DEFINED FFTW_INCLUDE_DIR )
    FIND_PATH(FFTW_INCLUDE_DIR fftw3.h
    	/usr/local/include
    	/usr/include
    	/usr/include/FFTW
    	$ENV{FFTW_INCLUDE_PATH}
    	$ENV{FFTW_ROOT}/include
    	${FFTW_ROOT}/include
    )
ENDIF( )

find_library ( FFTW_FLOAT_LIBRARY fftw3
	/usr/local/lib
	/usr/lib
	$ENV{FFTW_LIBRARY_PATH}
	$ENV{FFTW_ROOT}/lib
	${FFTW_ROOT}/lib
)

find_library ( FFTW_DOUBLE_LIBRARY fftw3f
    /usr/local/lib
    /usr/lib
    $ENV{FFTW_LIBRARY_PATH}
    $ENV{FFTW_ROOT}/lib
    ${FFTW_ROOT}/lib
)

find_library ( FFTW_LONG_DOUBLE_LIBRARY fftw3l
    /usr/local/lib
    /usr/lib
    $ENV{FFTW_LIBRARY_PATH}
    $ENV{FFTW_ROOT}/lib
    ${FFTW_ROOT}/lib
)

if ( FFTW_INCLUDE_DIR )

    set ( FFTW_LIBRARIES "" )
    set ( FFTW_TYPES_LIST "" )

	if ( FFTW_FLOAT_LIBRARY )
		list( APPEND FFTW_LIBRARIES ${FFTW_FLOAT_LIBRARY})
		set( FFTW_FOUND TRUE )
		list ( APPEND FFTW_TYPES_LIST "float" )
	endif ()
	
	if ( FFTW_DOUBLE_LIBRARY )
		list( APPEND FFTW_LIBRARIES ${FFTW_DOUBLE_LIBRARY})
		set( FFTW_FOUND TRUE )
		list ( APPEND FFTW_TYPES_LIST "double" )
	endif ()
	
	if ( FFTW_LONG_DOUBLE_LIBRARY )
		list( APPEND FFTW_LIBRARIES ${FFTW_LONG_DOUBLE_LIBRARY})
		set( FFTW_FOUND TRUE )
		list ( APPEND FFTW_TYPES_LIST "long double" )
	endif ()
	
endif ()

mark_as_advanced( FFTW_FOUND FFTW_INCLUDE_DIR FFTW_LIBRARIES
                  FFTW_FLOAT_LIBRARY FFTW_DOUBLE_LIBRARY FFTW_LONG_DOUBLE_LIBRARY )
