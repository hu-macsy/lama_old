###
 # @file Package/dl.cmake
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
 # @brief findPackage für lib dl
 # @author Lauretta Schubert
 # @date 05.02.2016
###

### SCAI_DL_LIBRARIES   - needed dl libraries

if    ( NOT WIN32 )

	find_library( SCAI_DL_LIBRARIES dl )
	
else  ( NOT WIN32 )
	set ( SCAI_DL_LIBRARIES "" )
endif ( NOT WIN32 )

#mark_as_advanced ( SCAI_DL_LIBRARIES )
