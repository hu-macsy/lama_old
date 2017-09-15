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
#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/NoCommunicator.hpp>

using namespace std;

namespace scai
{

SCAI_LOG_DEF_LOGGER( LAMAMPIBenchmark::logger, "Benchmark.MPIBenchmark" );

LAMAMPIBenchmark::LAMAMPIBenchmark()
    : benchmark::Benchmark()
{
    mComm = dmemo::Communicator::getCommunicatorPtr();
    SCAI_LOG_INFO( logger, "Communicator created for Benchmark: " << *mComm );
    readConfig();
}

LAMAMPIBenchmark::LAMAMPIBenchmark( const string& id, const string& gid )
    : benchmark::Benchmark( id, gid )
{
    mComm = dmemo::Communicator::getCommunicatorPtr();
    SCAI_LOG_INFO( logger,
                   "Communicator created for Benchmark (id = " << id << ", gid = " << gid << ") : " << *mComm );
    readConfig();
}

LAMAMPIBenchmark::LAMAMPIBenchmark( const LAMAMPIBenchmark& other )
    : benchmark::Benchmark( other ), mComm( other.mComm )
{
    readConfig();
}

static void tokenize( const string& str, vector<string>& tokens, const string& delimiters = " " )
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of( delimiters, lastPos );

    while ( string::npos != pos || string::npos != lastPos )
    {
        // Found a token, add it to the vector.
        tokens.push_back( str.substr( lastPos, pos - lastPos ) );
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of( delimiters, pos );
        // Find next "non-delimiter"
        pos = str.find_first_of( delimiters, lastPos );
    }
}

void LAMAMPIBenchmark::readConfig()
{
    const char* envConfig = getenv( "LAMA_CONFIG" );

    if ( envConfig != NULL )
    {
        config = envConfig;
    }

    std::vector<std::string> tokens;

    tokenize( config, tokens, "," );

    int ntoken = static_cast<int>( tokens.size() );

    int size = mComm->getSize();
    int rank = mComm->getRank();

    if ( ntoken == 0 )
    {
        config = "";
    }
    else if ( ntoken == 1 )
    {
        config = tokens[0];
    }
    else if ( ntoken == size )
    {
        config = tokens[rank];
    }
    else
    {
        COMMON_THROWEXCEPTION( "Mismatch of #processors = " << size << " and LAMA_CONFIG=" << config );
    }

    // for convenience: make upper case of config

    for ( size_t j = 0; j < config.length(); j++ )
    {
        config[j] = static_cast<std::string::value_type>( toupper( config[j] ) );
    }

    SCAI_LOG_INFO( logger, "Process " << rank << " of " << size << " has config = " << config );
}

void LAMAMPIBenchmark::getConfig( std::map<std::string, std::string>& tokens ) const
{
    std::vector<std::string> tmpTokens;
    tokenize( config, tmpTokens, ":" );
    typedef std::vector<std::string>::const_iterator VecIter;
    VecIter end = tmpTokens.end();

    for ( VecIter it = tmpTokens.begin(); it != end; ++it )
    {
        std::vector<std::string> keyValue;
        tokenize( *it, keyValue, "=" );
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

LAMAMPIBenchmark& LAMAMPIBenchmark::operator=( const LAMAMPIBenchmark& other )
{
    benchmark::Benchmark::operator ==( other );
    mComm = other.mComm;
    return *this;
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
    bool doOutput = true;

    if ( mComm->getRank() != 0 )
    {
        doOutput = false;
    }

    return doOutput;
}

}
