/**
 * @file macros/logging.hpp
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
 * @brief Macro definitions for logging
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

/*************************************************************************
 *                                                                       *
 *  Compile time guards for LOGGING                                      *
 *                                                                       *
 *    make sure that the desired logging levels are enabled              *
 *                                                                       *
 *    SCAI_LOG_DEBUG_ENABLED  :  SCAI_LOG_DEBUG is done                  *
 *    SCAI_LOG_INFO_ENABLED   :  SCAI_LOG_INFO is done                   *
 *                                                                       *
 *  The compile time guards itself can be set by these macros:           *
 *                                                                       *
 *    SCAI_LOG_LEVEL_TRACE  - compile all                                *
 *    SCAI_LOG_LEVEL_DEBUG  - compile debug and higher                   *
 *    SCAI_LOG_LEVEL_INFO   - compile info and higher                    *
 *    SCAI_LOG_LEVEL_WARN   - compile warn and higher                    *
 *    SCAI_LOG_LEVEL_ERROR  - compile error and higher                   *
 *    SCAI_LOG_LEVEL_FATAL  - compile fatal only                         *
 *    SCAI_LOG_LEVEL_OFF    - removes all LOG statements at compile time *
 *                                                                       *
 *  Please note: These guards are only for compile time, so logging      *
 *               can still be switched off at runtime.                   *
 *                                                                       *
 ************************************************************************/

#include <iostream>

/*******************************************************
 *   SCAI_LOG_LEVEL_TRACE                                    *
 *******************************************************/

#if defined(SCAI_LOG_LEVEL_TRACE)

#define SCAI_LOG_TRACE_ENABLED
#define SCAI_LOG_DEBUG_ENABLED
#define SCAI_LOG_INFO_ENABLED
#define SCAI_LOG_WARN_ENABLED
#define SCAI_LOG_ERROR_ENABLED
#define SCAI_LOG_FATAL_ENABLED

/*******************************************************
 *   SCAI_LOG_LEVEL_DEBUG                                    *
 *******************************************************/

#elif defined(SCAI_LOG_LEVEL_DEBUG)

#undef SCAI_LOG_TRACE_ENABLED

#define SCAI_LOG_DEBUG_ENABLED
#define SCAI_LOG_INFO_ENABLED
#define SCAI_LOG_WARN_ENABLED
#define SCAI_LOG_ERROR_ENABLED
#define SCAI_LOG_FATAL_ENABLED

/*******************************************************
 *   SCAI_LOG_LEVEL_INFO                                     *
 *******************************************************/

#elif defined(SCAI_LOG_LEVEL_INFO)

#undef SCAI_LOG_TRACE_ENABLED
#undef SCAI_LOG_DEBUG_ENABLED

#define SCAI_LOG_INFO_ENABLED
#define SCAI_LOG_WARN_ENABLED
#define SCAI_LOG_ERROR_ENABLED
#define SCAI_LOG_FATAL_ENABLED

/*******************************************************
 *   SCAI_LOG_LEVEL_WARN                                     *
 *******************************************************/

#elif defined(SCAI_LOG_LEVEL_WARN)

#undef SCAI_LOG_TRACE_ENABLED
#undef SCAI_LOG_DEBUG_ENABLED
#undef SCAI_LOG_INFO_ENABLED

#define SCAI_LOG_WARN_ENABLED
#define SCAI_LOG_ERROR_ENABLED
#define SCAI_LOG_FATAL_ENABLED

/*******************************************************
 *   SCAI_LOG_LEVEL_ERROR                                    *
 *******************************************************/

#elif defined(SCAI_LOG_LEVEL_ERROR)

#undef SCAI_LOG_TRACE_ENABLED
#undef SCAI_LOG_DEBUG_ENABLED
#undef SCAI_LOG_INFO_ENABLED
#undef SCAI_LOG_WARN_ENABLED

#define SCAI_LOG_ERROR_ENABLED
#define SCAI_LOG_FATAL_ENABLED

/*******************************************************
 *   SCAI_LOG_LEVEL_FATAL                                    *
 *******************************************************/

#elif defined(SCAI_LOG_LEVEL_FATAL)

#undef SCAI_LOG_DEBUG_ENABLED
#undef SCAI_LOG_TRACE_ENABLED
#undef SCAI_LOG_INFO_ENABLED
#undef SCAI_LOG_WARN_ENABLED
#undef SCAI_LOG_ERROR_ENABLED

#define SCAI_LOG_FATAL_ENABLED

/*******************************************************
 *   SCAI_LOG_LEVEL_OFF                                      *
 *******************************************************/

#elif defined(SCAI_LOG_LEVEL_OFF)

#undef SCAI_LOG_DEBUG_ENABLED
#undef SCAI_LOG_TRACE_ENABLED
#undef SCAI_LOG_INFO_ENABLED
#undef SCAI_LOG_WARN_ENABLED
#undef SCAI_LOG_ERROR_ENABLED
#undef SCAI_LOG_FATAL_ENABLED

#else

/*******************************************************
 *   DEFAULT: The Default SCAI_LOG_FATAL_ENABLED is enabled *
 *******************************************************/

// turned off for master branch
// #pragma message("Please define SCAI_LOG_LEVEL_xxx with xxx = TRACE, DEBUG, INFO, WARN, ERROR, FATAL, or OFF") 
// #pragma message("Will use default SCAI_LOG_LEVEL_FATAL")

#undef SCAI_LOG_DEBUG_ENABLED
#undef SCAI_LOG_TRACE_ENABLED
#undef SCAI_LOG_INFO_ENABLED
#undef SCAI_LOG_WARN_ENABLED
#undef SCAI_LOG_ERROR_ENABLED

#define SCAI_LOG_FATAL_ENABLED

#endif

#ifdef SCAI_LOG_FATAL_ENABLED

/*******************************************************
 *   logging enabled, at least one -DSCAI_LOG_LEVEL_xxx      *
 *******************************************************/

#include <scai/logging/Level.hpp>
#include <scai/logging/SourceLocation.hpp>

#include <scai/logging/Logger.hpp>
#include <scai/logging/LoggerProvider.hpp>

#include <scai/common/thread.hpp>

/*******************************************************
 *   Definitions for logging                           *
 *******************************************************/

#define SCAI_LOG_DECL_STATIC_LOGGER(aLogger) static class scai::logging::Logger& aLogger;
#define SCAI_LOG_DEF_LOGGER(aLogger,name) scai::logging::Logger& aLogger = \
        scai::logging::LoggerProvider::getProvider().getInstance(std::string(name));
#define SCAI_LOG_DEF_TEMPLATE_LOGGER(temp,aLogger,name) temp scai::logging::Logger& aLogger = \
        scai::logging::LoggerProvider::getProvider().getInstance(std::string(name));
#define SCAI_LOG_USING(alogger) using alogger;

#include <sstream>

/*******************************************************
 *   SCAI_LOG_XXXXX_ON : Predicates                    *
 *******************************************************/

#define SCAI_LOG_TRACE_ON(logger) (logger.isTraceEnabled())
#define SCAI_LOG_DEBUG_ON(logger) (logger.isDebugEnabled())
#define SCAI_LOG_INFO_ON(logger) (logger.isInfoEnabled())
#define SCAI_LOG_WARN_ON(logger) (logger.isWarnEnabled())
#define SCAI_LOG_ERROR_ON(logger) (logger.isErrorEnabled())
#define SCAI_LOG_FATAL_ON(logger) (logger.isFatalEnabled())

/*******************************************************
 *   SCAI_LOG_SET_XXXX                                 *
 *******************************************************/

#define SCAI_LOG_SET_TRACE(logger) { logger.setLevel(scai::logging::TRACE, true); }
#define SCAI_LOG_SET_DEBUG(logger) { logger.setLevel(scai::logging::DEBUG, true); }
#define SCAI_LOG_SET_INFO(logger) { logger.setLevel(scai::logging::INFO, true); }
#define SCAI_LOG_SET_WARN(logger) { logger.setLevel(scai::logging::WARN, true); }
#define SCAI_LOG_SET_ERROR(logger) { logger.setLevel(scai::logging::ERROR, true); }
#define SCAI_LOG_SET_FATAL(logger) { logger.setLevel(scai::logging::FATAL, true); }

#else

/*******************************************************
 *   logging completely disabled                       *
 *******************************************************/

#define SCAI_LOG_DECL_STATIC_LOGGER( aLogger )
#define SCAI_LOG_DEF_LOGGER( aLogger, name)
#define SCAI_LOG_USING( aLogger )
#define SCAI_LOG_DEF_TEMPLATE_LOGGER( temp, aLogger, name)

#define SCAI_LOG_TRACE_ON(logger) (false)
#define SCAI_LOG_DEBUG_ON(logger) (false)
#define SCAI_LOG_INFO_ON(logger)  (false)
#define SCAI_LOG_WARN_ON(logger)  (false)
#define SCAI_LOG_ERROR_ON(logger) (false)
#define SCAI_LOG_FATAL_ON(logger) (false)

#endif // SCAI_LOG_FATAL_ENABLED

/*******************************************************
 *   SCAI_LOG_TRACE                                    *
 *******************************************************/

#ifdef SCAI_LOG_TRACE_ENABLED
#define SCAI_LOG_TRACE( logger, msg )                       \
    {                                                       \
        scai::logging::Logger& cLogger = logger;            \
        \
        if ( cLogger.isTraceEnabled() )                     \
        {                                                   \
            std::ostringstream omsg;                        \
            omsg << msg;                                    \
            cLogger.trace( LOG4SCAI_LOCATION, omsg.str());  \
        }                                                   \
    }
#else
#define SCAI_LOG_TRACE( logger, msg )                       \
    {                                                       \
        if ( false )                                        \
        {                                                   \
            std::cout << msg;                               \
        }                                                   \
    }
#endif

/*******************************************************
 *   SCAI_LOG_DEBUG                                    *
 *******************************************************/

#ifdef SCAI_LOG_DEBUG_ENABLED
#define SCAI_LOG_DEBUG( logger, msg )                       \
    {                                                       \
        scai::logging::Logger& cLogger = logger;            \
        if ( cLogger.isDebugEnabled() )                     \
        {                                                   \
            std::ostringstream omsg;                        \
            omsg << msg;                                    \
            cLogger.debug( LOG4SCAI_LOCATION, omsg.str());  \
        }                                                   \
    }
#else
#define SCAI_LOG_DEBUG( logger, msg )                       \
    {                                                       \
        if ( false )                                        \
        {                                                   \
            std::cout << msg;                               \
        }                                                   \
    }
#endif

/*******************************************************
 *   SCAI_LOG_INFO                                     *
 *******************************************************/

#ifdef SCAI_LOG_INFO_ENABLED
#define SCAI_LOG_INFO( logger, msg )                        \
    {                                                       \
        scai::logging::Logger& cLogger = logger;            \
        if ( cLogger.isInfoEnabled() )                      \
        {                                                   \
            std::ostringstream omsg;                        \
            omsg << msg;                                    \
            cLogger.info( LOG4SCAI_LOCATION, omsg.str());   \
        }                                                   \
    }
#else
#define SCAI_LOG_INFO( logger, msg )                        \
    {                                                       \
        if ( false )                                        \
        {                                                   \
            std::cout << msg;                               \
        }                                                   \
    }
#endif

/*******************************************************
 *   SCAI_LOG_WARN                                     *
 *******************************************************/

#ifdef SCAI_LOG_WARN_ENABLED
#define SCAI_LOG_WARN( logger, msg )                        \
    {                                                       \
        scai::logging::Logger& cLogger = logger;            \
        if ( cLogger.isWarnEnabled() )                      \
        {                                                   \
            std::ostringstream omsg;                        \
            omsg << msg;                                    \
            cLogger.warn( LOG4SCAI_LOCATION, omsg.str());   \
        }                                                   \
    }
#else
#define SCAI_LOG_WARN( logger, msg )                        \
    {                                                       \
        if ( false )                                        \
        {                                                   \
            std::cout << msg;                               \
        }                                                   \
    }
#endif

/*******************************************************
 *   SCAI_LOG_ERROR                                    *
 *******************************************************/

#ifdef SCAI_LOG_ERROR_ENABLED
#define SCAI_LOG_ERROR( logger, msg )                       \
    {                                                       \
        scai::logging::Logger& cLogger = logger;            \
        if ( cLogger.isErrorEnabled() )                     \
        {                                                   \
            std::ostringstream omsg;                        \
            omsg << msg;                                    \
            cLogger.error( LOG4SCAI_LOCATION, omsg.str());  \
        }                                                   \
    }
#else
#define SCAI_LOG_ERROR( logger, msg )                       \
    {                                                       \
        if ( false )                                        \
        {                                                   \
            std::cout << msg;                               \
        }                                                   \
    }
#endif

/*******************************************************
 *   SCAI_LOG_FATAL                                    *
 *******************************************************/

#ifdef SCAI_LOG_FATAL_ENABLED
#define SCAI_LOG_FATAL( logger, msg )                       \
    {                                                       \
        scai::logging::Logger& cLogger = logger;            \
        if ( cLogger.isFatalEnabled() )                     \
        {                                                   \
            std::ostringstream omsg;                        \
            omsg << msg;                                    \
            cLogger.fatal( LOG4SCAI_LOCATION, omsg.str());  \
        }                                                   \
    }
#else
#define SCAI_LOG_FATAL( logger, msg )                       \
    {                                                       \
        if ( false )                                        \
        {                                                   \
            std::cout << msg;                               \
        }                                                   \
    }
#endif

/*******************************************************
 *   SCAI_LOG_THREAD                                   *
 *******************************************************/

#ifdef SCAI_LOG_LEVEL_OFF

#define SCAI_LOG_THREAD( name )                      \
    {                                                \
        if ( false )                                 \
        {                                            \
            std::cout << "";                         \
        }                                            \
    }

#else

#include <scai/common/thread.hpp>

// macro defines a name for the current thread

#define SCAI_LOG_THREAD( name )                                                \
    {                                                                          \
        std::ostringstream oname;                                              \
        oname << name;                                                         \
        scai::common::thread::defineCurrentThreadName( oname.str().c_str() );  \
    }

#endif
