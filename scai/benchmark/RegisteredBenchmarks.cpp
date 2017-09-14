/**
 * @file RegisteredBenchmarks.cpp
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
 * @brief RegisteredBenchmarks.cpp
 * @author Jiri Kraus
 * @date 06.04.2011
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

using namespace scai;

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
