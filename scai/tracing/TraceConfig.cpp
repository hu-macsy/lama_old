/**
 * @file TraceConfig.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Implementation of class TraceConfig
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include <scai/tracing/TraceConfig.hpp>

// local library
#include <scai/tracing/VTInterface.hpp>
#include <scai/tracing/TraceData.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/Settings.hpp>

// std
#include <iostream>
#include <cstdlib>
#include <memory>

using std::shared_ptr;

namespace scai
{

namespace tracing
{

/* -------------------------------------------------------------------------- *
 *   Static class variables                                                   *
 * -------------------------------------------------------------------------- */

shared_ptr<TraceConfig> TraceConfig::config;

bool TraceConfig::globalTraceFlag = true;

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( TraceConfig::logger, "TraceConfig" )

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
    mEnabled = false;
    mThreadEnabled = false;
    mVampirTraceEnabled = false;
    mTimeTraceEnabled = false;
    mCallTreeEnabled = false;
    mTraceFilePrefix = "_";
    // value of environmentvariable:  param1:param2=valx:param3:param4=valy
    std::vector<std::string> values;

    // get all values separated by :

    if ( scai::common::Settings::getEnvironment( values, SCAI_ENV_TRACE_CONFIG, ":" ) )
    {
        // SCAI_TRACE=key1:key2:key3=val3:key4
        if ( values.size() != 1 || values[0] != "OFF" )
        {
            mEnabled = true;

            for ( size_t i = 0; i < values.size(); ++i )
            {
                std::vector<std::string> keys;
                scai::common::Settings::tokenize( keys, values[i], "=" );
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
                       SCAI_ENV_TRACE_CONFIG << " not set, tracing is disabled."
                       << " Enable by " << SCAI_ENV_TRACE_CONFIG << "=time|ct[:vt][:thread]" )
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
    mMaster = std::this_thread::get_id();
    SCAI_LOG_INFO( logger, "ThreadConfig: enabled = " << mEnabled )
    SCAI_LOG_INFO( logger, "ThreadConfig: call tree enabled = " << mCallTreeEnabled )
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
    fileName << mTraceFilePrefix << ".time";
    /* For parallel processes we should add suffix for rank */
    int rank;

    if ( scai::common::Settings::getEnvironment( rank, "SCAI_RANK" ) )
    {
        fileName << "." << rank;
    }

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

static std::mutex mapMutex; // needed to avoid conflicts for accesses on mTraceDataMap

/* -------------------------------------------------------------------------- */

TraceData* TraceConfig::getTraceData()
{
    if ( !mEnabled )
    {
        return NULL;
    }

    std::thread::id self = std::this_thread::get_id();

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

TraceData* TraceConfig::getTraceData( std::thread::id threadId )
{
    // make sure that not two different threads try to allocate a table
    // read / write access at same time might also result in a crash
    std::unique_lock<std::mutex> lock( mapMutex );
    shared_ptr<TraceData> traceData;
    // Old stuff: shared_ptr<TraceData> traceData = mTraceDataMap[threadId];
    std::map<ThreadId, std::shared_ptr<TraceData> >::const_iterator it = mTraceDataMap.find( threadId );

    if ( it == mTraceDataMap.end() )
    {
        // this thread calls the first time a region
        try
        {
            traceData.reset( new TraceData( mTraceFilePrefix.c_str(), threadId, mThreadEnabled, mCallTreeEnabled ) );
        }
        catch ( common::Exception& ex )
        {
            std::cerr << "getTraceData: caugh exception, will continue without tracing" << std::endl;
            std::cerr << ex.what() << std::endl;
        }

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

} /* end namespace tracing */

} /* end namespace scai */
