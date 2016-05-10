###
 # @file Summaries/common.cmake
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
 # @brief SCAI common Configuration Summary
 # @author Lauretta Schubert
 # @date 25.08.2015
###

include ( Functions/scaiMessages )

emptyline()
message ( STATUS "=====================================" )
message ( STATUS "Summary of SCAI common Configuration:" )
message ( STATUS "=====================================" )

include ( Summaries/Modules/Compiler )

# common (core)
heading ( "Required core:" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR )
    set ( REQUIRED_FOUND TRUE )
endif ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR  )

heading2 ( "External Libraries" "REQUIRED_FOUND" )

    # pthreads
    found_message ( "pThreads" "SCAI_THREAD_LIBRARIES" "REQUIRED" "Version ${SCAI_THREAD_VERSION}" )
    # boost
    found_message ( "Boost" "SCAI_BOOST_INCLUDE_DIR" "REQUIRED" "Version ${BOOST_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

include ( Summaries/Modules/Build )

include ( Summaries/Modules/Configuration )
