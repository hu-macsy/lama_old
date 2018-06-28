###
 # @file relative_install.cmake
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
 # @brief Own cmake install function that keepss relative path names
 # @author Thomas Brandes
 # @date 01.10.2015
###

# alternative version to install ( FILES f1 f2 ... DESTINATION <dest> )
# this version keeps relative pathnames of f1 f2 ...

function ( relative_install )

   set ( options )
   set ( oneValueArgs DESTINATION )
   set ( multiValueArgs FILES )

   cmake_parse_arguments ( relative_install "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

   # message( STATUS "relative_install, FILES = ${relative_install_FILES}" )
   # message( STATUS "relative_install, DESTINATION = ${relative_install_DESTINATION}" )

   foreach   ( SOURCE_FILE ${relative_install_FILES} )

       if    ( CMAKE_VERSION VERSION_LESS 2.8.12 )
           get_filename_component( FILE_DIR ${SOURCE_FILE} PATH )
       else  ( CMAKE_VERSION VERSION_LESS 2.8.12 )
           get_filename_component( FILE_DIR ${SOURCE_FILE} DIRECTORY )
       endif ( CMAKE_VERSION VERSION_LESS 2.8.12 )

       # message( STATUS "install ${SOURCE_FILE} in ${relative_install_DESTINATION}/${FILE_DIR}" )

       install( FILES ${SOURCE_FILE} DESTINATION ${relative_install_DESTINATION}/${FILE_DIR} )

   endforeach ( )

endfunction ( relative_install )
