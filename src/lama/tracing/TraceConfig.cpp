/**
 * @file TraceConfig.cpp
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
 * @brief Implementation of class TraceConfig
 * @author: Thomas Brandes
 * @date 01.01.2012
 * $Id$
 */

// hpp
#include <lama/tracing/TraceConfig.hpp>
#include <lama/tracing/VTInterface.hpp>

// others
#include <lama/tracing/RegionTable.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/ContextFactory.hpp>
#include <lama/ContextAccess.hpp>

#include <cstdlib>
#include <iostream>

using namespace lama;

namespace tracing
{

/* -------------------------------------------------------------------------- */

boost::shared_ptr<TraceConfig> TraceConfig::config;

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( TraceConfig::logger, "TraceConfig" );

/* -------------------------------------------------------------------------- */

static std::vector<std::string> split( const std::string& params, const char seperator )
{
    std::vector<std::string> args;
    size_t found = std::string::npos;

    do
    {
        const size_t prevFound = found + 1;
        found = params.find( seperator, prevFound );
        args.push_back( params.substr( prevFound, found - prevFound ) );
    }

    while ( found != std::string::npos );

    return args;
}

/* -------------------------------------------------------------------------- */

boost::shared_ptr<TraceConfig> TraceConfig::getInstancePtr()
{
    if ( !config )
    {
        config.reset( new TraceConfig() );
    }

    return config;
}

/* -------------------------------------------------------------------------- */

#if defined( LAMA_TRACE_LEVEL_VT )
extern "C" void VT_User_trace_on__();
extern "C" void VT_User_trace_off__();
#endif

/* -------------------------------------------------------------------------- */

void TraceConfig::setKey( const std::string& key, const std::string& value )
{
    LAMA_LOG_INFO( logger, "Set trace key " << key << " = " << value );

    if ( key == "PREFIX" )
    {
        mTraceFilePrefix = value;
    }
    else
    {
        LAMA_LOG_WARN( logger, key << " is unknown key for TRACE configuration" );
    }
}

/* -------------------------------------------------------------------------- */

void TraceConfig::enableVampirTrace( bool flag )
{
    mVampirTraceEnabled = flag;

#if defined( LAMA_TRACE_LEVEL_VT )
    if ( mCUDAContext )
    {
        // enable the context when tracing is switched on / off

        LAMA_LOG_INFO( logger, "enableVampirTrace: flag =  " << flag
                       << ", context = " << *mCUDAContext );

        LAMA_CONTEXT_ACCESS( mCUDAContext );

        VTInterface::enable( flag );
    }
    else
    {
        LAMA_LOG_INFO( logger, "enableVampirTrace: flag =  " << flag << ", no context" );

        VTInterface::enable( flag );
    }
#else
    LAMA_LOG_INFO( logger, "enableVampirTrace: flag =  " << flag << ", level VT disabled" );

    if ( mVampirTraceEnabled )
    {
        LAMA_LOG_WARN( logger, "TRACE:vt ignored, use LAMA_TRACE_LEVEL=VT for compilation." );
    }
#endif
}

/* -------------------------------------------------------------------------- */

void TraceConfig::setParam( const std::string& param )
{
    LAMA_LOG_INFO( logger, "Set trace config value : " << param );

    if ( param == "VT" )
    {
        mVampirTraceEnabled = true;
    }
    else if ( param == "THREAD" )
    {
        mThreadEnabled = true;
    }
    else if ( param == "TIME" )
    {
        mTimeTraceEnabled = true;
    }
    else
    {
        LAMA_LOG_WARN( logger, param << " is unknown option for TRACE" );
    }
}

/* -------------------------------------------------------------------------- */

TraceConfig::TraceConfig()
{
    mComm = CommunicatorFactory::get( "MPI" );
    // mCUDAContext = ContextFactory::getContext( Context::CUDA );

    // should be done explicitly
    // mCUDAContext = ContextFactory::getContext( Context::CUDA );

    mEnabled = false;
    mThreadEnabled = false;
    mVampirTraceEnabled = false;
    mTimeTraceEnabled = false;
    mTraceFilePrefix = "_";

    if ( getenv( LAMA_ENV_TRACE_CONFIG ) )
    {
        std::string params = getenv( LAMA_ENV_TRACE_CONFIG );

        std::vector<std::string> values = split( params, ':' );

        if ( values.size() != 1 || values[0] != "OFF" )
        {
            mEnabled = true;

            for ( size_t i = 0; i < values.size(); ++i )
            {
                std::vector<std::string> keys = split( values[i], '=' );

                std::string& key = keys[0];

                for ( size_t j = 0; j < key.length(); j++ )
                {
                    key[j] = static_cast<std::string::value_type>( toupper( key[j] ) );
                }

                if ( keys.size() == 1 )
                {
                    setParam( key );
                }
                else
                {
                    setKey( key, keys[1] );
                }
            }
        }
    }
    else
    {
        LAMA_LOG_WARN( logger, "LAMA_TRACE not set, tracing is disabled. Enable by LAMA_TRACE=time[:vt][:thread]" );
    }

    // enable/disable VampirTrace, action needed now

    enableVampirTrace( mVampirTraceEnabled );

    if ( mTraceFilePrefix == "_" )
    {
        // no prefix specified, so take default

        if ( getenv( "_" ) )
        {
            mTraceFilePrefix = getenv( "_" );
        }
        else
        {
            mTraceFilePrefix = "LAMA";
        }
    }
    else
    {
        // prefix set, use it also for VampirTrace
#if defined( LAMA_TRACE_LEVEL_VT )
        if ( mVampirTraceEnabled )
        {
            // 0 for do not overwrite existing value
            setenv( "VT_FILE_PREFIX", mTraceFilePrefix.c_str(), 0 );
        }
#endif
    }

    // save id of this main thread

    mMaster = Thread::getSelf();

    LAMA_LOG_INFO( logger, "ThreadConfig: enabled = " << mEnabled );
}

/* -------------------------------------------------------------------------- */

TraceConfig::~TraceConfig()
{
    LAMA_LOG_DEBUG( logger, "Entering Destructor." );

    if ( mVampirTraceEnabled )
    {
        // set tracing explicitly to off now

        enableVampirTrace( false );
    }

    if ( !mTimeTraceEnabled )
    {
        LAMA_LOG_INFO( logger, "~TraceConfig, no output file" );
        LAMA_LOG_DEBUG( logger, "Leaving Destructor." );
        return;
    }

    // now print all info in a file

    std::ostringstream fileName;

    fileName << mTraceFilePrefix << ".trace";

    if ( mComm->getSize() > 1 )
    {
        fileName << "." << mComm->getRank();
    }

    LAMA_LOG_INFO( logger, "~TraceConfig, output file = " << fileName.str() );

    FILE* f = fopen( fileName.str().c_str(), "w" );

    if ( f == NULL )
    {
        LAMA_LOG_ERROR( logger, "Could not open " << fileName.str() << " for writing time information." );
        return;
    }

    // Now write each RegionTable in file

    std::map<lama::Thread::Id,boost::shared_ptr<RegionTable> >::iterator it;

    for ( it = mRegionTables.begin(); it != mRegionTables.end(); it++ )
    {
        RegionTable& table = *it->second;
        table.printTimer( f );
    }

    fclose( f );

    LAMA_LOG_DEBUG( logger, "Leaving Destructor." );
}

/* -------------------------------------------------------------------------- */

static boost::mutex mapMutex; // needed to avoid conflicts for accesses on mRegionTables

/* -------------------------------------------------------------------------- */

RegionTable* TraceConfig::getRegionTable()
{
    if ( !mEnabled )
    {
        return NULL;
    }

    Thread::Id self = Thread::getSelf();

    if ( mThreadEnabled || self == mMaster )
    {
        return getInstance().getRegionTable( self );
    }
    else
    {
        return NULL;
    }
}

/* -------------------------------------------------------------------------- */

RegionTable* TraceConfig::getRegionTable( Thread::Id threadId )
{
    boost::shared_ptr<RegionTable> regionTable = mRegionTables[threadId];

    if ( !regionTable )
    {
        boost::mutex::scoped_lock scoped_lock( mapMutex );

        // this thread calls the first time a region

        regionTable.reset( new RegionTable( threadId ) );

        mRegionTables[threadId] = regionTable;
    }

    return regionTable.get();
}

} // namespace tracing

