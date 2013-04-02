/**
 * @file LoggerCreator.hpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief LoggerCreator.hpp
 * @author Kai Buschulte
 * @date 01.08.2012
 * $Id$
 */

#ifndef LAMA_LOGGERCREATOR_HPP_
#define LAMA_LOGGERCREATOR_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

#include <logging/logging.hpp>

//Timer
#include <lama/solver/logger/OpenMPTimer.hpp>

//Logger
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/solver/logger/FileLogger.hpp>

// boost
#include <boost/spirit/include/qi.hpp>
#include <boost/variant.hpp>

#include <string>
#include <ostream>

namespace lama
{

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

class LAMA_DLL_IMPORTEXPORT LoggerCreator
{
public:
    typedef qi::rule<std::string::const_iterator,LoggerPtr(),ascii::space_type> RuleType;
    typedef qi::symbols<char,LoggerPtr> LoggerInstanceMap;

    static RuleType& getSolverBoundRule();

    static qi::rule<std::string::const_iterator,void(),ascii::space_type>& getIndependentRule();

protected:
    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private:
    static LoggerCreator& getInstance();

    /**
     * @brief Add/register a logger with a given name/id to the Factory
     *
     * @param[in] name      The name to identify the instance
     * @param[in] logger  The instance that will be registered
     */
    void addLogger( const std::string& name, LoggerPtr logger );

    qi::symbols<char, LogLevel::LogLevel> mLogLevel;
    qi::symbols<char, LoggerWriteBehaviour::LoggerWriteBehaviour> mWriteBehaviour;

    qi::rule<std::string::const_iterator, std::string(), ascii::space_type> mRLoggerName;
    qi::rule<std::string::const_iterator, std::string(), ascii::space_type> mRFileName;
    qi::rule<std::string::const_iterator, Timer*(), ascii::space_type> mRTimer;

    RuleType mRSolverBoundLogger;
    qi::rule<std::string::const_iterator, std::string(), ascii::space_type> mRId;
    qi::rule<std::string::const_iterator, void(), ascii::space_type> mRIndependentLogger;
    qi::rule<std::string::const_iterator, Logger*(), ascii::space_type> mRCommonLogger;
    qi::rule<std::string::const_iterator, Logger*(), ascii::space_type> mRFileLogger;

    LoggerInstanceMap mLoggerInstanceMap;

    LoggerCreator();
};

} // namespace lama

#endif // LAMA_LOGGERCREATOR_HPP_
