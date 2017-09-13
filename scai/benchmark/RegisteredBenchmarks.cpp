/**
 * @file RegisteredBenchmarks.cpp
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
 * @brief RegisteredBenchmarks.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
 */
/*
 * RegisteredBenchmarks.cpp
 *
 *  Created on: 03.02.2011
 *      Author: rrehrman
 */
#include <scai/benchmark.hpp>
#include <scai/benchmark/BenchmarkPrinter.hpp>
#include <scai/benchmark/frame_stdlib.hpp>

extern "C" bf::BenchmarkRegistry* getBenchmarkRegistry();

int main( void )
{
    std::vector<std::string> files;

    try
    {
        bf::getSharedLibraries( files );
    }
    catch( std::exception& e )
    {
        bf::BenchmarkPrinter::error( e.what() );
        return 1;
    }

    std::map<std::string,std::string> benchmarks;
    bf::BenchmarkRegistry* registry;

    for( unsigned int i = 0; i < files.size(); ++i )
    {
        typedef bf::BenchmarkRegistry* (*registry_t)();
        registry_t reg_handle = NULL;

        LAMA_LIB_HANDLE_TYPE handle;

        int error = bf::loadLibAndGetFunctionHandle( reg_handle, handle, files[i].c_str(), "getBenchmarkRegistry" );

        if( error != 0 )
        {
            continue;
        }

        typedef void (*free_registry_t)();
        free_registry_t reg_free_handle;

        error = bf::getFunctionHandle( reg_free_handle, handle, "releaseBenchmarkRegistry" );
        if( error != 0 )
        {
            std::stringstream message;
            message << files[i] << " does not define releaseBenchmarkRegistry skipping it.";
            bf::BenchmarkPrinter::warning( message.str() );
            //bf::freeLibHandle( handle );
            continue;
        }

        // getting registry from library.
        registry = reg_handle();

        bf::Benchmark* bench;

        // iterate over registry and save all BenchmarkIDs and Names, if those
        // Benchmarks need any unknown arguments, insert '<parametered>', instead.
        for( bf::BenchmarkRegistry::const_iterator it = registry->begin(); it != registry->end(); ++it )
        {
            try
            {
                bench = it->second->create();
            }
            catch( bf::BFError& )
            {
                benchmarks.insert( std::make_pair( it->first, "<parametered>" ) );
                continue;
            }
            // insert ID and Name.
            benchmarks.insert( std::make_pair( bench->getId(), bench->getName() ) );
            delete bench;
        }

        reg_free_handle();
        //bf::freeLibHandle( handle );
    }

    std::stringstream message;
    for( std::map<std::string,std::string>::const_iterator it = benchmarks.begin(); it != benchmarks.end(); ++it )
    {
        message << it->first << "%_:_%" << it->second << "%,";
    }

    std::string idsAndNames = message.str();

    // remove the last '%,'
    bf::BenchmarkPrinter::print( idsAndNames.substr( 0, idsAndNames.size() - 2 ) );
}
