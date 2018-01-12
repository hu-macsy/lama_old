/**
 * @file TraceConfig.hpp
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
 * @brief Definition of class that specifies the runtime configuration for tracing.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>

// interal scai libraries
#include <scai/logging.hpp>

// std
#include <string>
#include <map>
#include <memory>

namespace scai
{

namespace tracing
{

class TraceData;

/** Name of environment variable used to specify trace configuration. */

#define SCAI_ENV_TRACE_CONFIG "SCAI_TRACE"

/**
 * @brief This class is used to define/set the runtime configuration for tracing.
 *
 * This class is a singleton. The only one object holds member variables that specify
 * the current settings for tracing, e.g. what kind of information is collected or
 * whether it is enabled at all.
 *
 * The singleton object will be allocated on demand after program start. A static
 * object cannot be used as MPI or VampirTrace might be uninitialized.
 *
 * The destructor of the trace configuration is used to write out the collected
 * information at runtime, e.g. the collected time information of the regions.
 */

class COMMON_DLL_IMPORTEXPORT TraceConfig: private common::NonCopyable
{
public:

    /** Type definition for thread identification used in tracing. */

    typedef std::thread::id ThreadId;

    ~TraceConfig();

    /**
     * @brief Get reference to the actual trace configuration.
     *
     * @return the only instance of this class.
     */

    static TraceConfig& getInstance()
    {
        return *getInstancePtr();
    }

    /**
     * Get the actual trace configuration as a shared pointer.
     */
    static std::shared_ptr<TraceConfig> getInstancePtr();

    /**
     * Query if tracing is enabled.
     */
    bool isEnabled()
    {
        return mEnabled;
    }

    bool isThreadEnabled( ThreadId threadId )
    {
        return mThreadEnabled || threadId == mMaster;
    }

    bool isVampirTraceEnabled()
    {
        return mVampirTraceEnabled;
    }

    bool isTimeTraceEnabled()
    {
        return mTimeTraceEnabled;
    }

    bool isCallTreeEnabled()
    {
        return mCallTreeEnabled;
    }

    const char* getFilePrefix() const
    {
        return getInstance().mTraceFilePrefix.c_str();
    }

    /** Get trace data for the current thread.
     *
     *  @return pointers to the TraceData for the calling thread (might be NULL)
     */

    TraceData* getTraceData();

    void traceOff();

    /** Helper class for setting global trace flag in a scope
     *  (Constructor sets global flag, destructor resets it)
     */
    class TraceScope
    {
    public:

        /** Constructor sets global trace flag and saves old state. */

        TraceScope( bool flag );

        /** Destructor resets saved state of global trace flag. */

        ~TraceScope();

    private:

        bool saveFlag; // used to save state of global trace flag
    };

    static bool globalTraceFlag;

private:

    TraceConfig();

    void setParam( const std::string& param );

    void setKey( const std::string& key, const std::string& value );

    ThreadId mMaster; //!< id of mather thread

    bool mEnabled;

    bool mTimeTraceEnabled;

    bool mCallTreeEnabled;

    bool mVampirTraceEnabled;

    void enableVampirTrace( bool flag ); // switch trace on / off

    bool mThreadEnabled; //!< true if trace should also be done for threads

    std::string mTraceFilePrefix;

    /** Each thread will have its own table for region timing.
     *  Use of shared pointer for entry in map
     */

    std::map<ThreadId, std::shared_ptr<TraceData> > mTraceDataMap;

    /** Get the trace data by the id of a thread. */

    TraceData* getTraceData( ThreadId threadId );

    /** The only one instance allocated at program start. */

    static std::shared_ptr<TraceConfig> config;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace tracing */

} /* end namespace scai */
