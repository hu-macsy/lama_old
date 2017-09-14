/**
 * @file frame_stdlib.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief frame_stdlib.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
 */
/*
 * frame_stdlib.cpp
 *
 *  Created on: 01.02.2011
 *      Author: rrehrman
 */

#include <cstdlib>

#include <scai/benchmark/frame_stdlib.hpp>
#include <scai/benchmark/BFError.hpp>
#include <scai/benchmark/readDir.hpp>

namespace bf
{

void getSharedLibraries( std::vector<std::string>& files )
{
    std::string benchmark_library_path;
    char* val = getenv( "BENCHMARK_LIBRARY_PATH" );
    if( val != NULL )
    {
        benchmark_library_path = val;
    }
    if( benchmark_library_path.empty() )
    {
        throw BFError( "Set BENCHMARK_LIBRARY_PATH to the directory holding the "
                       "shared libraries of your benchmarks." );
    }

    // collect shared libraries from path.
    getFilesFromPath( benchmark_library_path, files );

    if( files.empty() )
    {
        throw BFError( "No shared libraries found in BENCHMARK_LIBRARY_PATH" );
    }
}

} // namespace bf
