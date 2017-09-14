/**
 * @file RunBenchmark.cpp
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
 * @brief RunBenchmark.cpp
 * @author Robin Rehrmann
 * @date 04.05.2010
 */

#include <vector>
#include <cstring>
#include <cstdlib>

#include <scai/benchmark.hpp>

#include <scai/benchmark/frame_stdlib.hpp>
#include <scai/benchmark/BenchmarkPrinter.hpp>

using namespace scai;

extern "C" bf::BenchmarkRegistry* getBenchmarkRegistry();
extern "C" bf::BaseInputSetRegistry* getInputSetRegistry();

int main( int argc, const char* argv[] )
{
    // For testing manually.
    if( argc != 7 )
    {
        std::ostringstream message;
        message << "Wrong number of parameters. <BenchmarkID> <InputSetID> <minTime> <numRepititions> <path> <tmp>"
                << std::endl;
        message << "Got:" << std::endl;
        for( int i = 0; i < argc; ++i )
        {
            message << "'" << argv[i] << "' " << std::flush;
        }
        bf::BenchmarkPrinter::error( message.str() );
        return 2;
    }

    std::string benchId = argv[1];
    std::string inputSetId = argv[2];
    float minTime = static_cast<float>( atof( argv[3] ) );
    int numRep = atoi( argv[4] );
    std::string path = argv[5];
    std::string tmp = argv[6];

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

    unsigned int i = 0;
    bool foundBenchmark = false;
    bf::BenchmarkRegistry* registry = NULL;
    LAMA_LIB_HANDLE_TYPE handle = NULL;
    typedef void (*free_registry_t)();
    free_registry_t reg_free_handle;

    // try dlopen on the shared libraries, until the registry, holding the
    // requested benchmark, is loaded.
    do
    {
        typedef bf::BenchmarkRegistry* (*registry_t)();

        registry_t reg_handle = NULL;

        int error = bf::loadLibAndGetFunctionHandle( reg_handle, handle, files[i].c_str(), "getBenchmarkRegistry" );
        if( error != 0 )
        {
            ++i;
            continue;
        }

        error = bf::getFunctionHandle( reg_free_handle, handle, "releaseBenchmarkLibraryResources" );
        if( error != 0 )
        {
            std::stringstream message;
            message << files[i] << " does not define releaseBenchmarkLibraryResources skipping it.";
            bf::BenchmarkPrinter::warning( message.str() );
            //bf::freeLibHandle( handle );
            ++i;
            continue;
        }

        // getting registry from library.
        registry = reg_handle();

        foundBenchmark = registry->has( benchId );
        ++i;

        if( !foundBenchmark )
        {
            reg_free_handle();
            //bf::freeLibHandle( handle );
            if( i >= files.size() )
            {
                std::stringstream message;
                message << '\'' << benchId << "' not found.";
                bf::BenchmarkPrinter::warning( message.str() );
                return 1;
            }
        }
    } while( !foundBenchmark );

    bf::Config& conf = bf::Config::getInstance();
    conf.setValueFor( "path", path );
    conf.setValueFor( "tmp", tmp );

    //! bench *must* be released _before_ handle is closed!
    //! Segmentation fault, otherwise!
    std::auto_ptr<bf::Benchmark> bench;
    try
    {
        bench = registry->createBenchmark( benchId );
        bf::BenchmarkPrinter::setDoOutput( bench->doOutput() );
        bench->setInputSetId( inputSetId );
        bench->setMinTime( minTime );
        bench->setNumRepitions( numRep );
    }
    catch( std::exception& e )
    {
        registry->destroyBenchmark( bench.release() );
        //bench.reset( 0 );
        bf::BenchmarkPrinter::warning( e.what() );
        reg_free_handle();
        //bf::freeLibHandle( handle );
        return 1;
    }
    try
    {
        bench->run( std::cout );
    }
    catch( bf::BFException& be )
    {
        registry->destroyBenchmark( bench.release() );
        //bench.reset( 0 );
        bf::BenchmarkPrinter::warning( be.what() );
        reg_free_handle();
        //bf::freeLibHandle( handle );
        return 1;
    }
    catch( std::exception& e )
    {
        registry->destroyBenchmark( bench.release() );
        //bench.reset( 0 );
        bf::BenchmarkPrinter::error( e.what() );
        reg_free_handle();
        //bf::freeLibHandle( handle );
        return 1;
    }

    std::stringstream message;
    message << bench->getName() << "%," << bench->getInputSetId() << "%," << bench->getGid() << "%,"
            << bench->getNumThreads() << "%," << bench->getValueTypeSize() << "%," << bench->getExecutionFlops()
            << "%," << bench->getExecutionBandwidth() << "%," << bench->getSetupTime() << "%,"
            << bench->getMinExecutionTime() << "%," << bench->getMaxExecutionTime() << "%,"
            << bench->getExecutionTime() << "%," << bench->getTearDownTime() << "%,"
            << bench->getTotalExecutionTime() << "%," << bench->getExecutedTime() << "%,"
            << bench->getNumActualRepititons();

    bf::BenchmarkPrinter::print( message.str() );

    registry->destroyBenchmark( bench.release() );
    //bench.reset( 0 );
    reg_free_handle();

    //bf::freeLibHandle( handle );
    return 0;
}
