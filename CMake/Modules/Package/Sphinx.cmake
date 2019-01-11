###
 # @file Sphinx.cmake
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
 # @brief Package Sphinx (might be enabled or disabled)
 # @author Jan Ecker
 # @date 11.11.2015
###

### SPHINX_FOUND

include ( scai_macro/scai_pragma_once )
include ( scai_macro/scai_build_variable )
include ( scai_macro/scai_summary )

## run this configuration only once to avoid multiple summaries, messages

scai_pragma_once ()

find_package ( Sphinx ${SCAI_FIND_PACKAGE_FLAGS} )

mark_as_advanced( Sphinx_DIR )

scai_build_variable ( NAME      USE_SPHINX
                      BOOL 
                      DEFAULT   ${SPHINX_FOUND}
                      DOCSTRING "use of Sphinx (for user documentation)" )

scai_summary_external ( NAME       Sphinx
                        ENABLED    ${USE_SPHINX}
                        FOUND      ${SPHINX_FOUND}
                        VERSION    ${Sphinx_VERSION_STRING} 
                        EXECUTABLE ${Sphinx-build_EXECUTABLE} ${Sphinx-apidoc_EXECUTABLE} )
