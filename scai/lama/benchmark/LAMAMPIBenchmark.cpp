/**
 * @file LAMAMPIBenchmark.cpp
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
 * @brief LAMAMPIBenchmark.cpp
 * @author Jiri Kraus
 * @date 13.05.2011
 */

#include <scai/common/OpenMP.hpp>

#include <scai/lama/benchmark/LAMAMPIBenchmark.hpp>
#include <scai/benchmark/Parser.hpp>

#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/NoCommunicator.hpp>

using namespace std;

namespace scai
{

SCAI_LOG_DEF_LOGGER( LAMAMPIBenchmark::logger, "Benchmark.MPIBenchmark" );

LAMAMPIBenchmark::LAMAMPIBenchmark( const string& name, const string& gId ) : benchmark::Benchmark( name, gId )

{
    mComm = dmemo::Communicator::getCommunicatorPtr();
    SCAI_LOG_INFO( logger,
                   "Communicator created for Benchmark( name = " << name << ", group = " << gId << ") : " << *mComm );
    readConfig();
}

void LAMAMPIBenchmark::readConfig()
{
    const char* envConfig = getenv( "LAMA_CONFIG" );

    if ( envConfig != NULL )
    {
        mConfig = envConfig;
    }

    std::vector<std::string> tokens;

    tokenize( tokens, mConfig, "," );

    int ntoken = static_cast<int>( tokens.size() );

    int size = mComm->getSize();
    int rank = mComm->getRank();

    if ( ntoken == 0 )
    {
        mConfig = "";
    }
    else if ( ntoken == 1 )
    {
        mConfig = tokens[0];
    }
    else if ( ntoken == size )
    {
        mConfig = tokens[rank];
    }
    else
    {
        COMMON_THROWEXCEPTION( "Mismatch of #processors = " << size << " and LAMA_CONFIG=" << mConfig );
    }

    // for convenience: make upper case of mConfig

    for ( size_t j = 0; j < mConfig.length(); j++ )
    {
        mConfig[j] = static_cast<std::string::value_type>( toupper( mConfig[j] ) );
    }

    SCAI_LOG_INFO( logger, "Process " << rank << " of " << size << " has config = " << mConfig );
}

void LAMAMPIBenchmark::getConfig( std::map<std::string, std::string>& tokens ) const
{
    std::vector<std::string> tmpTokens;

    tokenize( tmpTokens, mConfig, ":" );

    typedef std::vector<std::string>::const_iterator VecIter;
    VecIter end = tmpTokens.end();

    for ( VecIter it = tmpTokens.begin(); it != end; ++it )
    {
        std::vector<std::string> keyValue;
        tokenize( keyValue, *it, "=" );
        tokens[keyValue[0]] = keyValue[1];
    }
}

void LAMAMPIBenchmark::synchronize() const
{
    mComm->synchronize();
}

LAMAMPIBenchmark::~LAMAMPIBenchmark()
{
    mComm.reset();
}

std::string LAMAMPIBenchmark::getNumThreads() const
{
    std::stringstream threads;
    int numThreads = 0;
    #pragma omp parallel
    {
        #pragma omp master
        {
            numThreads = omp_get_num_threads();
        }
    }

    // sum up number of threads of all processors

    numThreads = mComm->sum( numThreads );
    threads << numThreads;

    SCAI_LOG_DEBUG( logger, "getNumThreads = " << numThreads );

    return threads.str();
}

bool LAMAMPIBenchmark::doOutput() const
{
    bool doOutput = mComm->getRank() == 0;  // only proc 0

    return doOutput;
}

}
