/**
 * @file logging.hpp
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

#pragma message("Please define SCAI_LOG_LEVEL_xxx with xxx = TRACE, DEBUG, INFO, WARN, ERROR, FATAL, or OFF")
#pragma message("Will use default SCAI_LOG_LEVEL_FATAL")

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

#include <scai/common/Thread.hpp>

/*******************************************************
 *   Definitions for logging                           *
 *******************************************************/

#define SCAI_LOG_DECL_STATIC_LOGGER(aLogger) static class logging::Logger& aLogger;
#define SCAI_LOG_DEF_LOGGER(aLogger,name) logging::Logger& aLogger = \
        logging::LoggerProvider::getProvider().getInstance(std::string(name));
#define SCAI_LOG_DEF_TEMPLATE_LOGGER(temp,aLogger,name) temp logging::Logger& aLogger = \
        logging::LoggerProvider::getProvider().getInstance(std::string(name));
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

#define SCAI_LOG_SET_TRACE(logger) { logger.setLevel(logging::TRACE, true); }
#define SCAI_LOG_SET_DEBUG(logger) { logger.setLevel(logging::DEBUG, true); }
#define SCAI_LOG_SET_INFO(logger) { logger.setLevel(logging::INFO, true); }
#define SCAI_LOG_SET_WARN(logger) { logger.setLevel(logging::WARN, true); }
#define SCAI_LOG_SET_ERROR(logger) { logger.setLevel(logging::ERROR, true); }
#define SCAI_LOG_SET_FATAL(logger) { logger.setLevel(logging::FATAL, true); }

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
#define SCAI_LOG_TRACE(logger,msg) { if (&logger && logger.isTraceEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.trace(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define SCAI_LOG_TRACE(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   SCAI_LOG_DEBUG                                    *
 *******************************************************/

#ifdef SCAI_LOG_DEBUG_ENABLED
#define SCAI_LOG_DEBUG(logger,msg) { if (&logger && logger.isDebugEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.debug(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define SCAI_LOG_DEBUG(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   SCAI_LOG_INFO                                     *
 *******************************************************/

#ifdef SCAI_LOG_INFO_ENABLED
#define SCAI_LOG_INFO(logger,msg) { if (&logger && logger.isInfoEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.info(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define SCAI_LOG_INFO(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   SCAI_LOG_WARN                                     *
 *******************************************************/

#ifdef SCAI_LOG_WARN_ENABLED
#define SCAI_LOG_WARN(logger,msg) { if (&logger && logger.isWarnEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.warn(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define SCAI_LOG_WARN(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   SCAI_LOG_ERROR                                    *
 *******************************************************/

#ifdef SCAI_LOG_ERROR_ENABLED
#define SCAI_LOG_ERROR(logger,msg) { if (&logger && logger.isErrorEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.error(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define SCAI_LOG_ERROR(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   SCAI_LOG_FATAL                                    *
 *******************************************************/

#ifdef SCAI_LOG_FATAL_ENABLED
#define SCAI_LOG_FATAL(logger,msg) { if (&logger && logger.isFatalEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.fatal(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define SCAI_LOG_FATAL(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   SCAI_LOG_THREAD                                   *
 *******************************************************/

#ifdef SCAI_LOG_LEVEL_OFF

#define SCAI_LOG_THREAD                          \
{                                                \
    if ( false )                                 \
    {                                            \
        std::cout << "";                         \
    }                                            \
}

#else

// macro defines a name for the current thread 

#define SCAI_LOG_THREAD( name )                                      \
{                                                                    \
    std::ostringstream oname;                                        \
    oname << name;                                                   \
    common::Thread::defineCurrentThreadName( oname.str().c_str() );  \
}

#endif

// just temporarily include all files
#include <scai/logging/AbstractLoggerCreator.hpp>
#include <scai/logging/GenLogger.hpp>
#include <scai/logging/GenLoggerCreator.hpp>
#include <scai/logging/Level.hpp>
#include <scai/logging/Logger.hpp>
#include <scai/logging/LoggerProvider.hpp>
#include <scai/logging/SourceLocation.hpp>
