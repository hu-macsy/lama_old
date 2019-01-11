###
 # @file Package/Doxygen.cmake
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
 # @brief findPackage and configuration of doxygen
 # @author Jan Ecker
 # @date 25.04.2013
###

### DOXYGEN_FOUND
### DOXYGEN_DOT_EXECUTABLE
### DOXYGEN_EXECUTABLE

find_package ( Doxygen ${SCAI_FIND_PACKAGE_FLAGS} )

scai_build_variable ( NAME      USE_DOXYGEN
                      BOOL 
                      DEFAULT   ${DOXYGEN_FOUND}
                      DOCSTRING "use of doxygen (for system documentation)" )

scai_summary_external ( NAME       Doyxgen
                        FOUND      ${DOXYGEN_FOUND}
                        ENABLED    ${USE_DOXYGEN} 
                        VERSION    ${DOXYGEN_VERSION} 
                        EXECUTABLE ${DOXYGEN_EXECUTABLE} )
