/**
 * @file logging.hpp
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
 * @brief Macro definitions for logging within LAMA
 * @author brandes
 * @date 01.03.2011
 * $Id$
 */

#ifndef LAMA_LOGGING_HPP_
#define LAMA_LOGGING_HPP_

/************************************************************************
 *                                                                       *
 *  Compile time guards for LOGGING                                      *
 *                                                                       *
 *    make sure that the desired logging levels are enabled              *
 *                                                                       *
 *    LAMA_LOG_DEBUG_ENABLED  :  LAMA_LOG_DEBUG is done                            *
 *    LAMA_LOG_INFO_ENABLED   :  LAMA_LOG_INFO is done                             *
 *                                                                       *
 *  The compile time guards itself can be set by these macros:           *
 *                                                                       *
 *    LAMA_LOG_LEVEL_TRACE  - compile all                                     *
 *    LAMA_LOG_LEVEL_DEBUG  - compile debug and higher                        *
 *    LAMA_LOG_LEVEL_INFO   - compile info and higher                         *
 *    LAMA_LOG_LEVEL_WARN   - compile warn and higher                         *
 *    LAMA_LOG_LEVEL_ERROR  - compile error and higher                        *
 *    LAMA_LOG_LEVEL_FATAL  - compile fatal only                              *
 *    LAMA_LOG_LEVEL_OFF    - removes all LOG statements at compile time      *
 *                                                                       *
 *  Please note: These guards are only for compile time, so logging      *
 *               can still be switched off at runtime.                   *
 *                                                                       *
 ************************************************************************/

#include <iostream>

/*******************************************************
 *   LAMA_LOG_LEVEL_TRACE                                    *
 *******************************************************/

#if defined(LAMA_LOG_LEVEL_TRACE)

#define LAMA_LOG_TRACE_ENABLED
#define LAMA_LOG_DEBUG_ENABLED
#define LAMA_LOG_INFO_ENABLED
#define LAMA_LOG_WARN_ENABLED
#define LAMA_LOG_ERROR_ENABLED
#define LAMA_LOG_FATAL_ENABLED

/*******************************************************
 *   LAMA_LOG_LEVEL_DEBUG                                    *
 *******************************************************/

#elif defined(LAMA_LOG_LEVEL_DEBUG)

#undef LAMA_LOG_TRACE_ENABLED

#define LAMA_LOG_DEBUG_ENABLED
#define LAMA_LOG_INFO_ENABLED
#define LAMA_LOG_WARN_ENABLED
#define LAMA_LOG_ERROR_ENABLED
#define LAMA_LOG_FATAL_ENABLED

/*******************************************************
 *   LAMA_LOG_LEVEL_INFO                                     *
 *******************************************************/

#elif defined(LAMA_LOG_LEVEL_INFO)

#undef LAMA_LOG_TRACE_ENABLED
#undef LAMA_LOG_DEBUG_ENABLED

#define LAMA_LOG_INFO_ENABLED
#define LAMA_LOG_WARN_ENABLED
#define LAMA_LOG_ERROR_ENABLED
#define LAMA_LOG_FATAL_ENABLED

/*******************************************************
 *   LAMA_LOG_LEVEL_WARN                                     *
 *******************************************************/

#elif defined(LAMA_LOG_LEVEL_WARN)

#undef LAMA_LOG_TRACE_ENABLED
#undef LAMA_LOG_DEBUG_ENABLED
#undef LAMA_LOG_INFO_ENABLED

#define LAMA_LOG_WARN_ENABLED
#define LAMA_LOG_ERROR_ENABLED
#define LAMA_LOG_FATAL_ENABLED

/*******************************************************
 *   LAMA_LOG_LEVEL_ERROR                                    *
 *******************************************************/

#elif defined(LAMA_LOG_LEVEL_ERROR)

#undef LAMA_LOG_TRACE_ENABLED
#undef LAMA_LOG_DEBUG_ENABLED
#undef LAMA_LOG_INFO_ENABLED
#undef LAMA_LOG_WARN_ENABLED

#define LAMA_LOG_ERROR_ENABLED
#define LAMA_LOG_FATAL_ENABLED

/*******************************************************
 *   LAMA_LOG_LEVEL_FATAL                                    *
 *******************************************************/

#elif defined(LAMA_LOG_LEVEL_FATAL)

#undef LAMA_LOG_DEBUG_ENABLED
#undef LAMA_LOG_TRACE_ENABLED
#undef LAMA_LOG_INFO_ENABLED
#undef LAMA_LOG_WARN_ENABLED
#undef LAMA_LOG_ERROR_ENABLED

#define LAMA_LOG_FATAL_ENABLED

/*******************************************************
 *   LAMA_LOG_LEVEL_OFF                                      *
 *******************************************************/

#elif defined(LAMA_LOG_LEVEL_OFF)

#undef LAMA_LOG_DEBUG_ENABLED
#undef LAMA_LOG_TRACE_ENABLED
#undef LAMA_LOG_INFO_ENABLED
#undef LAMA_LOG_WARN_ENABLED
#undef LAMA_LOG_ERROR_ENABLED
#undef LAMA_LOG_FATAL_ENABLED

#else

/*******************************************************
 *   DEFAULT: The Default LAMA_LOG_FATAL_ENABLED is enabled *
 *******************************************************/

#pragma message("Please define LAMA_LOG_LEVEL_xxx with xxx = TRACE, DEBUG, INFO, WARN, ERROR, FATAL, or OFF")
#pragma message("Will use default LAMA_LOG_LEVEL_FATAL")

#undef LAMA_LOG_DEBUG_ENABLED
#undef LAMA_LOG_TRACE_ENABLED
#undef LAMA_LOG_INFO_ENABLED
#undef LAMA_LOG_WARN_ENABLED
#undef LAMA_LOG_ERROR_ENABLED

#define LAMA_LOG_FATAL_ENABLED

#endif

#ifdef LAMA_LOG_FATAL_ENABLED

/*******************************************************
 *   logging enabled, at least one -DLAMA_LOG_LEVEL_xxx      *
 *******************************************************/

#include "logging/Level.hpp"
#include "logging/SourceLocation.hpp"

#include "logging/Logger.hpp"
#include "logging/LoggerProvider.hpp"

/*******************************************************
 *   Defintions for logging                             *
 *******************************************************/

#define LAMA_LOG_DECL_STATIC_LOGGER(aLogger) static class log4lama::Logger& aLogger
#define LAMA_LOG_DEF_LOGGER(aLogger,name) log4lama::Logger& aLogger = \
        log4lama::LoggerProvider::getProvider().getInstance(std::string(name))
#define LAMA_LOG_DEF_TEMPLATE_LOGGER(temp,aLogger,name) temp log4lama::Logger& aLogger = \
        log4lama::LoggerProvider::getProvider().getInstance(std::string(name))
#define LAMA_LOG_USING(alogger) using alogger

#include <sstream>

/*******************************************************
 *   LAMA_LOG_XXXXX_ON : Predicates                          *
 *******************************************************/

#define LAMA_LOG_TRACE_ON(logger) (logger.isTraceEnabled())
#define LAMA_LOG_DEBUG_ON(logger) (logger.isDebugEnabled())
#define LAMA_LOG_INFO_ON(logger) (logger.isInfoEnabled())
#define LAMA_LOG_WARN_ON(logger) (logger.isWarnEnabled())
#define LAMA_LOG_ERROR_ON(logger) (logger.isErrorEnabled())
#define LAMA_LOG_FATAL_ON(logger) (logger.isFatalEnabled())

/*******************************************************
 *   LAMA_LOG_SET_XXXX                                       *
 *******************************************************/

#define LAMA_LOG_SET_TRACE(logger) { logger.setLevel(log4lama::TRACE, true); }
#define LAMA_LOG_SET_DEBUG(logger) { logger.setLevel(log4lama::DEBUG, true); }
#define LAMA_LOG_SET_INFO(logger) { logger.setLevel(log4lama::INFO, true); }
#define LAMA_LOG_SET_WARN(logger) { logger.setLevel(log4lama::WARN, true); }
#define LAMA_LOG_SET_ERROR(logger) { logger.setLevel(log4lama::ERROR, true); }
#define LAMA_LOG_SET_FATAL(logger) { logger.setLevel(log4lama::FATAL, true); }

#else

/*******************************************************
 *   logging completely disabled                        *
 *******************************************************/

#define LAMA_LOG_DECL_STATIC_LOGGER( aLogger )
#define LAMA_LOG_DEF_LOGGER( aLogger, name)
#define LAMA_LOG_USING( aLogger )
#define LAMA_LOG_DEF_TEMPLATE_LOGGER( temp, aLogger, name)

#define LAMA_LOG_TRACE_ON(logger) (false)
#define LAMA_LOG_DEBUG_ON(logger) (false)
#define LAMA_LOG_INFO_ON(logger)  (false)
#define LAMA_LOG_WARN_ON(logger)  (false)
#define LAMA_LOG_ERROR_ON(logger) (false)
#define LAMA_LOG_FATAL_ON(logger) (false)

#endif // LAMA_LOG_FATAL_ENABLED
/*******************************************************
 *   LAMA_LOG_TRACE                                          *
 *******************************************************/

#ifdef LAMA_LOG_TRACE_ENABLED
#define LAMA_LOG_TRACE(logger,msg) { if (&logger && logger.isTraceEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.trace(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define LAMA_LOG_TRACE(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   LAMA_LOG_DEBUG                                          *
 *******************************************************/

#ifdef LAMA_LOG_DEBUG_ENABLED
#define LAMA_LOG_DEBUG(logger,msg) { if (&logger && logger.isDebugEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.debug(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define LAMA_LOG_DEBUG(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   LAMA_LOG_INFO                                           *
 *******************************************************/

#ifdef LAMA_LOG_INFO_ENABLED
#define LAMA_LOG_INFO(logger,msg) { if (&logger && logger.isInfoEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.info(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define LAMA_LOG_INFO(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   LAMA_LOG_WARN                                           *
 *******************************************************/

#ifdef LAMA_LOG_WARN_ENABLED
#define LAMA_LOG_WARN(logger,msg) { if (&logger && logger.isWarnEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.warn(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define LAMA_LOG_WARN(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   LAMA_LOG_ERROR                                          *
 *******************************************************/

#ifdef LAMA_LOG_ERROR_ENABLED
#define LAMA_LOG_ERROR(logger,msg) { if (&logger && logger.isErrorEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.error(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define LAMA_LOG_ERROR(logger,msg) { if (false){ std::cout<<msg; } }
#endif

/*******************************************************
 *   LAMA_LOG_FATAL                                          *
 *******************************************************/

#ifdef LAMA_LOG_FATAL_ENABLED
#define LAMA_LOG_FATAL(logger,msg) { if (&logger && logger.isFatalEnabled()) \
        { std::ostringstream omsg; omsg << msg; logger.fatal(LOG4LAMA_LOCATION, omsg.str()); } }
#else
#define LAMA_LOG_FATAL(logger,msg) { if (false){ std::cout<<msg; } }
#endif

#endif // LAMA_LOGGING_HPP_
