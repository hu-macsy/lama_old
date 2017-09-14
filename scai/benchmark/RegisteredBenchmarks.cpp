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
 * @author Jiri Kraus, Robin Rehrmann
 * @date 03.02.2011
 */

#include <scai/benchmark.hpp>
#include <scai/benchmark/BenchmarkPrinter.hpp>
#include <scai/benchmark/Benchmark.hpp>
#include <scai/benchmark/frame_stdlib.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/LibModule.hpp>
#include <scai/common/unique_ptr.hpp>

#include <vector>

using namespace scai;
using namespace std;

int main( void )
{
    using namespace bf;

    std::string benchLibPath;

    if ( !common::Settings::getEnvironment( benchLibPath, "BENCHMARK_LIBRARY_PATH" ) )
    {
        throw BFError( "Set BENCHMARK_LIBRARY_PATH to the directory holding the "
                       "shared libraries of your benchmarks." );
    }

    std::map<std::string,std::string> benchmarks;   // keep all pairs of id, gid

    common::LibModule::loadLibsInDir( benchLibPath.c_str() );

    vector<string> values;  // string is create type for benchmarks

    Benchmark::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        cout << "  Registered values[" << i << "] = " << values[i] << endl;

        common::unique_ptr<Benchmark> bench( Benchmark::create( values[i] ) );

        cout << "  Id = " << bench->getId() << ", Name = " << bench->getName() << endl;

        benchmarks.insert( std::make_pair( bench->getId(), bench->getName() ) );
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
