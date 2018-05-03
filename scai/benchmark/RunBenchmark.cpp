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
 * @author Thomas Brandes
 * @date 14.09.2017
 */

#include <vector>
#include <cstring>
#include <cstdlib>

#include <scai/benchmark.hpp>

#include <scai/benchmark/BenchmarkPrinter.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/LibModule.hpp>
#include <scai/common/OpenMP.hpp>

using namespace scai;
using namespace benchmark;

int main( int argc, const char* argv[] )
{
    // For testing manually.

    if ( argc != 7 )
    {
        std::ostringstream message;
        message << "Wrong number of parameters. <BenchmarkID> <InputSetID> <minTime> <numRepititions>"
                << std::endl;
        message << "Got:" << std::endl;
        for( int i = 0; i < argc; ++i )
        {
            message << "'" << argv[i] << "' " << std::flush;
        }

        BenchmarkPrinter::error( message.str() );
        return 2;
    }

    std::string benchId = argv[1];
    std::string inputSetId = argv[2];
    float minTime = static_cast<float>( atof( argv[3] ) );
    int numRep = atoi( argv[4] );

    // arg 5 is the input path used as prefix for input files with relative pathname

    scai::common::Settings::putEnvironment( "SCAI_INPUT_PATH", argv[5], true );

    std::string benchLibPath;

    if ( !common::Settings::getEnvironment( benchLibPath, "BENCHMARK_LIBRARY_PATH" ) )
    {
        // try the official one if SCAI_ROOT has been set

        if ( common::Settings::getEnvironment( benchLibPath, "SCAI_ROOT" ) )
        {
            benchLibPath += "/benchmark";
        }
        else
        {
            throw BFError( "Set BENCHMARK_LIBRARY_PATH to the directory holding the "
                           "shared libraries of your benchmarks." );
        }
    }

    common::LibModule::loadLibsInDir( benchLibPath.c_str() );

    std::unique_ptr<Benchmark> bench;

    try
    {
        bench.reset( Benchmark::createWithArgument( benchId ) );
        benchmark::BenchmarkPrinter::setDoOutput( bench->doOutput() );
        bench->setInputSetId( inputSetId );
        bench->setMinTime( minTime );
        bench->setNumRepitions( numRep );
    }
    catch( std::exception& e )
    {
        benchmark::BenchmarkPrinter::warning( e.what() );
        return 1;
    }

    try
    {
        bench->run( std::cout );
    }
    catch( benchmark::BFException& be )
    {
        benchmark::BenchmarkPrinter::warning( be.what() );
        return 1;
    }
    catch( std::exception& e )
    {
        benchmark::BenchmarkPrinter::error( e.what() );
        return 1;
    }

    std::stringstream message;

    message << bench->getName() << "%," << bench->getInputSetId() << "%," << bench->getGroup() << "%,"
            << bench->getNumThreads() << "%," << common::typeSize( bench->getValueType() ) << "%," << bench->getExecutionFlops()
            << "%," << bench->getExecutionBandwidth() << "%," << bench->getSetupTime() << "%,"
            << bench->getMinExecutionTime() << "%," << bench->getMaxExecutionTime() << "%,"
            << bench->getExecutionTime() << "%," << bench->getTearDownTime() << "%,"
            << bench->getTotalExecutionTime() << "%," << bench->getExecutedTime() << "%,"
            << bench->getNumActualRepititons();

    BenchmarkPrinter::print( message.str() );

    return 0;
}
