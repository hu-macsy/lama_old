/**
 * @file TraceConfig.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @date 11.06.2015
 */

// hpp
#include <scai/tracing/TraceConfig.hpp>

// others
#include <scai/tracing/VTInterface.hpp>
#include <scai/tracing/TraceData.hpp>

#include <iostream>
#include <cstdlib>

namespace tracing
{

using common::Thread;
using common::shared_ptr;

/* -------------------------------------------------------------------------- *
 *   Static class variables                                                   *
 * -------------------------------------------------------------------------- */

shared_ptr<TraceConfig> TraceConfig::config;

bool TraceConfig::globalTraceFlag = true;

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( TraceConfig::logger, "TraceConfig" )

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

shared_ptr<TraceConfig> TraceConfig::getInstancePtr()
{
    if ( !config )
    {
        config.reset( new TraceConfig() );
    }

    return config;
}

/* -------------------------------------------------------------------------- */

#if defined( USE_VAMPIRTRACE )
extern "C" void VT_User_trace_on__();
extern "C" void VT_User_trace_off__();
#endif

/* -------------------------------------------------------------------------- */

void TraceConfig::setKey( const std::string& key, const std::string& value )
{
    SCAI_LOG_INFO( logger, "Set trace key " << key << " = " << value )

    if ( key == "PREFIX" )
    {
        mTraceFilePrefix = value;
    }
    else
    {
        SCAI_LOG_WARN( logger, key << " is unknown key for TRACE configuration" )
    }
}

/* -------------------------------------------------------------------------- */

void TraceConfig::enableVampirTrace( bool flag )
{
    mVampirTraceEnabled = flag;
#if defined( USE_VAMPIRTRACE )
    SCAI_LOG_INFO( logger, "enableVampirTrace: flag =  " << flag << ", no context" )
    VTInterface::enable( flag );
#else
    SCAI_LOG_INFO( logger, "enableVampirTrace: flag =  " << flag << ", VAMPIRTRACE not used" )

    if ( mVampirTraceEnabled )
    {
        SCAI_LOG_WARN( logger, "TRACE:vt ignored, define USE_VAMPIRTRACE for compilation." )
    }

#endif
}

/* -------------------------------------------------------------------------- */

void TraceConfig::setParam( const std::string& param )
{
    SCAI_LOG_INFO( logger, "Set trace config value : " << param )

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
    else if ( param == "CT" )
    {
        mCallTreeEnabled = true;
    }
    else
    {
        SCAI_LOG_WARN( logger, param << " is unknown option for TRACE" )
    }
}

/* -------------------------------------------------------------------------- */

TraceConfig::TraceConfig()
{
    // mComm = CommunicatorFactory::get();
    // mCUDAContext = ContextFactory::getContext( Context::CUDA );
    // should be done explicitly
    // mCUDAContext = ContextFactory::getContext( Context::CUDA );
    mEnabled = false;
    mThreadEnabled = false;
    mVampirTraceEnabled = false;
    mTimeTraceEnabled = false;
    mCallTreeEnabled = false;
    mTraceFilePrefix = "_";
    // value of environmentvariable:  param1:param2=valx:param3:param4=valy
    std::string params;
    const char* env = getenv( LAMA_ENV_TRACE_CONFIG );

    if ( env )
    {
        params = env;
        std::vector<std::string> values = split( params, ':' );

        if ( values.size() != 1 || values[0] != "OFF" )
        {
            mEnabled = true;

            for ( size_t i = 0; i < values.size(); ++i )
            {
                std::vector<std::string> keys = split( values[i], '=' );
                std::string& key = keys[0];

                // make upper case of key

                for ( size_t j = 0; j < key.length(); j++ )
                {
                    key[j] = static_cast<std::string::value_type>( toupper( key[j] ) );
                }

                if ( keys.size() == 1 )
                {
                    // is just a param
                    setParam( key );
                }
                else
                {
                    // is param=val
                    setKey( key, keys[1] );
                }
            }
        }
    }
    else
    {
        SCAI_LOG_WARN( logger,
                       LAMA_ENV_TRACE_CONFIG << " not set, tracing is disabled."
                       << " Enable by " << LAMA_ENV_TRACE_CONFIG << "=time|ct[:vt][:thread]" )
    }

    // enable/disable VampirTrace, action needed now
    enableVampirTrace( mVampirTraceEnabled );

    if ( mTraceFilePrefix == "_" )
    {
        // no prefix specified, so take default
        mTraceFilePrefix = "LAMA";
        // environment variable "_" contains name of last command
        const char* env = getenv( "_" );

        if ( env )
        {
            mTraceFilePrefix = env;
        }
    }
    else
    {
        // prefix set, use it also for VampirTrace
#if defined( USE_VAMPIRTRACE )
        if ( mVampirTraceEnabled )
        {
            // 0 for do not overwrite existing value
            setenv( "VT_FILE_PREFIX", mTraceFilePrefix.c_str(), 0 );
        }

#endif
    }

    // save id of this main thread
    mMaster = Thread::getSelf();
    SCAI_LOG_INFO( logger, "ThreadConfig: enabled = " << mEnabled )
}

/* -------------------------------------------------------------------------- */

TraceConfig::~TraceConfig()
{
    SCAI_LOG_DEBUG( logger, "Entering Destructor." )

    if ( mVampirTraceEnabled )
    {
        // set tracing explicitly to off now
        enableVampirTrace( false );
    }

    if ( !mTimeTraceEnabled )
    {
        SCAI_LOG_INFO( logger, "~TraceConfig, no output file" )
        SCAI_LOG_DEBUG( logger, "Leaving Destructor." )
        return;
    }

    // now print all info in a file
    std::ostringstream fileName;
    fileName << mTraceFilePrefix << ".trace";
    /* For parallel processes we should add suffix for rank

    if( mComm->getSize() > 1 )
    {
        fileName << "." << mComm->getRank();
    }

    */
    SCAI_LOG_INFO( logger, "~TraceConfig, output file = " << fileName.str() )
    std::ofstream outfile;
    outfile.open( fileName.str().c_str(), std::ios::out );

    if ( outfile.fail() )
    {
        // do not throw exception as only tracing caused problem
        SCAI_LOG_ERROR( logger, "Could not open " << fileName.str() << " for writing time information." )
        return;
    }

    // Now write each TraceData in file
    std::map<ThreadId, shared_ptr<TraceData> >::iterator it;

    for ( it = mTraceDataMap.begin(); it != mTraceDataMap.end(); it++ )
    {
        TraceData& data = *it->second;
        data.printTimer( outfile );
    }

    outfile.close();
    SCAI_LOG_DEBUG( logger, "~TraceConfig finished" )
}

/* -------------------------------------------------------------------------- */

static Thread::Mutex mapMutex; // needed to avoid conflicts for accesses on mTraceDataMap

/* -------------------------------------------------------------------------- */

TraceData* TraceConfig::getTraceData()
{
    if ( !mEnabled )
    {
        return NULL;
    }

    ThreadId self = common::Thread::getSelf();

    if ( mThreadEnabled || self == mMaster )
    {
        return getInstance().getTraceData( self );
    }
    else
    {
        return NULL;
    }
}

/* -------------------------------------------------------------------------- */

TraceData* TraceConfig::getTraceData( ThreadId threadId )
{
    // make sure that not two different threads try to allocate a table
    // read / write access at same time might also result in a crash

    Thread::ScopedLock lock( mapMutex );

    shared_ptr<TraceData> traceData;

    // Old stuff: shared_ptr<TraceData> traceData = mTraceDataMap[threadId];

    std::map<ThreadId, common::shared_ptr<TraceData> >::const_iterator it = mTraceDataMap.find( threadId );

    if ( it == mTraceDataMap.end() )
    {
        // this thread calls the first time a region

        traceData.reset( new TraceData( threadId, mThreadEnabled ) );
        mTraceDataMap.insert( std::pair<ThreadId, shared_ptr<TraceData> >( threadId, traceData ) );
    }
    else
    {
        traceData = it->second;
    }

    return traceData.get();
}

/* -------------------------------------------------------------------------- */

TraceConfig::TraceScope::TraceScope( bool flag )
{
    saveFlag = TraceConfig::globalTraceFlag;
    TraceConfig::globalTraceFlag = flag;
}

/** Destructor resets saved state of global trace flag. */

TraceConfig::TraceScope::~TraceScope()
{
    TraceConfig::globalTraceFlag = saveFlag;
}

} // namespace tracing

