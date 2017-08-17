###
 # @file scai_pragma_once.cmake
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
 # @brief Macro to build unit test for SCAI modules
 # @author Thomas Brandes
 # @date 13.07.2017
###

## This macro provides a similiar functionality as the C preprocessor directive.
##
## Mostly it does not matter whether a cmake configuration file is called twice 
## but there might be some inconveniences:
##
##   - if the configuration appends variables, e.g. summaries
##   - if the configuration prints warnings 

macro ( scai_pragma_once )

    ## Identify the file where this macro has been called

    get_filename_component ( name ${CMAKE_CURRENT_LIST_FILE} NAME_WE )

    if ( SCAI_CHECK_DONE_${name} )

        ## return works for the calling file, not for the macro only

        return ()

    else ()

        ## set a variable with a unique name, pass it to the parent scope if available

        set ( SCAI_CHECK_DONE_${name} True )

        get_directory_property ( hasParent PARENT_DIRECTORY )

        if ( hasParent )
            set ( SCAI_CHECK_DONE_${name} True PARENT_SCOPE )
        endif ()

    endif ()

endmacro ()
