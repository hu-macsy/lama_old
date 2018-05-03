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
 * @brief Print all benchmark classes registered in the Benchmark factory.
 * @author Thomas Brandes
 * @date 15.09.2017
 */

#include <scai/benchmark.hpp>
#include <scai/logging.hpp>
#include <scai/benchmark/BenchmarkPrinter.hpp>
#include <scai/benchmark/Benchmark.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/LibModule.hpp>

#include <vector>

using namespace scai;
using namespace std;

SCAI_LOG_DEF_LOGGER( logger, "RegisteredBenchmarks" )

int main( void )
{
    using namespace benchmark;

    std::string benchLibPath;

    if ( !common::Settings::getEnvironment( benchLibPath, "BENCHMARK_LIBRARY_PATH" ) )
    {
        throw BFError( "Set BENCHMARK_LIBRARY_PATH to the directory holding the "
                       "shared libraries of your benchmarks." );
    }

    common::LibModule::loadLibsInDir( benchLibPath.c_str() );

    vector<string> values;  // string is create type for benchmarks

    Benchmark::getCreateValues( values );

    std::stringstream message;

    for ( size_t i = 0; i < values.size(); ++i )
    {
        std::unique_ptr<Benchmark> bench( Benchmark::create( values[i], "" ) );
  
        SCAI_LOG_INFO( logger, "Benchmark " << i << " of " << values.size() << ": key = " << values[i] 
                                << ", Id = " << bench->getCreateId() << ", Name = " << bench->getName() )

        if ( i != 0 )
        {
            message << "%,";  // add separator after previous output
        }

        message << bench->getCreateId() << "%_:_%" << bench->getArgument();
    }

    benchmark::BenchmarkPrinter::print( message.str() );
}
