/**
 * @file RunBenchmark.cpp
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
 * @brief RunBenchmark.cpp
 * @author robin
 * @date 06.04.2011
 * $Id$
 */
/**
 * @file RunBenchmark.cpp
 * @author robin
 * Created on: 04.05.2010
 */

#include <vector>
#include <cstring>
#include <cstdlib>

#include <scai/benchmark.hpp>

#include <scai/benchmark/frame_stdlib.hpp>
#include <scai/benchmark/BenchmarkPrinter.hpp>

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
