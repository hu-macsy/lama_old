/**
 * @file LAMAMPIBenchmark.cpp
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
 * @brief LAMAMPIBenchmark.cpp
 * @author Jiri Kraus
 * @date 13.05.2011
 * $Id$
 */
#include <bench/LAMAMPIBenchmark.hpp>
#include <lama/Communicator.hpp>
#include <lama/NoCommunicator.hpp>

using namespace std;

LAMA_LOG_DEF_LOGGER( LAMAMPIBenchmark::logger, "Benchmark.MPIBenchmark" );

LAMAMPIBenchmark::LAMAMPIBenchmark()
    : bf::Benchmark()
{
    mComm = lama::CommunicatorFactory::get( "MPI" );
    LAMA_LOG_INFO( logger, "MPI communicator created for Benchmark: " << *mComm );
    readConfig();
}

LAMAMPIBenchmark::LAMAMPIBenchmark( const string& id, const string& gid )
    : bf::Benchmark( id, gid )
{
    mComm = lama::CommunicatorFactory::get( "MPI" );
    LAMA_LOG_INFO( logger,
                   "MPI communicator created for Benchmark (id = " << id << ", gid = " << gid << ") : " << *mComm );
    readConfig();
}

LAMAMPIBenchmark::LAMAMPIBenchmark( const string& id, const string& gid, const string& arguments )
    : bf::Benchmark( id, gid, arguments )
{
    mComm = lama::CommunicatorFactory::get( "MPI" );
    LAMA_LOG_INFO( logger,
                   "MPI communicator created for Benchmark (id = " << id << ", gid = " << gid << ", args = " << arguments << ") : " << *mComm );
    readConfig();
}

LAMAMPIBenchmark::LAMAMPIBenchmark( const LAMAMPIBenchmark& other )
    : bf::Benchmark( other ), mComm( other.mComm )
{
    readConfig();
}

static void tokenize( const string& str, vector<string>& tokens, const string& delimiters = " " )
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of( delimiters, 0 );
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of( delimiters, lastPos );

    while( string::npos != pos || string::npos != lastPos )
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

    if( envConfig != NULL )
    {
        config = envConfig;
    }

    std::vector<std::string> tokens;

    tokenize( config, tokens, "," );

    int ntoken = static_cast<int>( tokens.size() );

    int size = mComm->getSize();
    int rank = mComm->getRank();

    if( ntoken == 0 )
    {
        config = "";
    }
    else if( ntoken == 1 )
    {
        config = tokens[0];
    }
    else if( ntoken == size )
    {
        config = tokens[rank];
    }
    else
    {
        LAMA_THROWEXCEPTION( "Mismatch of #processors = " << size << " and LAMA_CONFIG=" << config );
    }

    // for convenience: make upper case of config

    for( size_t j = 0; j < config.length(); j++ )
    {
        config[j] = static_cast<std::string::value_type>( toupper( config[j] ) );
    }

    LAMA_LOG_INFO( logger, "Process " << rank << " of " << size << " has config = " << config );
}

void LAMAMPIBenchmark::getConfig( std::map<std::string,std::string>& tokens ) const
{
    std::vector<std::string> tmpTokens;
    tokenize( config, tmpTokens, ":" );
    typedef std::vector<std::string>::const_iterator VecIter;
    VecIter end = tmpTokens.end();
    for( VecIter it = tmpTokens.begin(); it != end; ++it )
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
    bf::Benchmark::operator ==( other );
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

    LAMA_LOG_DEBUG( logger, "getNumThreads = " << numThreads );

    return threads.str();
}

bool LAMAMPIBenchmark::doOutput() const
{
    bool doOutput = true;
    if( mComm->getRank() != 0 )
    {
        doOutput = false;
    }
    return doOutput;
}
