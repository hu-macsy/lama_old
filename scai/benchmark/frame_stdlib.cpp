/**
 * @file frame_stdlib.cpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief frame_stdlib.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
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
