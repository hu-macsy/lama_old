/**
 * @file LoggerCreator.cpp
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
 * @brief LoggerCreator.cpp
 * @author Kai Buschulte
 * @date 29.06.2012
 * $Id$
 */

// hpp
#include <lama/solver/creator/LoggerCreator.hpp>

// spirit
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( LoggerCreator::logger, "Solver.SolverCreator.LoggerCreator" );

LoggerCreator::LoggerCreator()
{
    namespace phoenix = boost::phoenix;
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    using qi::lit;
    using qi::_val;
    using qi::_1;
    using qi::_2;
    using qi::_3;
    using qi::_4;
    using qi::_5;

    using qi::on_error;
    using qi::fail;
    using qi::debug;

    using qi::double_;
    using qi::int_;

    using ascii::char_;

    using phoenix::val;
    using phoenix::construct;
    using phoenix::new_;
    using phoenix::bind;

    mRLoggerName = +char_( "0-9a-zA-Z<>" )[_val += _1];

    mRFileName = +char_( "0-9a-zA-Z." )[_val += _1];

    mLogLevel.add( "noLogging", LogLevel::noLogging )( "convergenceHistory", LogLevel::convergenceHistory )(
        "solverInformation", LogLevel::solverInformation )( "advancedInformation",
                LogLevel::advancedInformation )(
                    "completeInformation", LogLevel::completeInformation );

    mWriteBehaviour.add( "toConsoleOnly", LoggerWriteBehaviour::toConsoleOnly )(
        "toFileAndConsole", LoggerWriteBehaviour::toFileAndConsole )( "toFileOnly",
                LoggerWriteBehaviour::toFileOnly );

    mRTimer = lit( "OpenMPTimer" )[_val = new_<OpenMPTimer>()];

    mRCommonLogger = '('
                     > ( mRLoggerName > ',' > mLogLevel > ',' > mWriteBehaviour > ',' > mRTimer )[_val = new_<
                             CommonLogger>( _1, _2, _3, phoenix::construct<std::auto_ptr<Timer> >( _4 ) )] > ')';

    mRFileLogger =
        '('
        > ( mRLoggerName > ',' > mLogLevel > ',' > mWriteBehaviour > ',' > mRFileName > ','
            > mRTimer )[_val = new_<CommonLogger>(
                                   _1, _2, _3, _4, phoenix::construct<std::auto_ptr<Timer> >( _5 ) )]
        > ')';

    mRSolverBoundLogger = ( mLoggerInstanceMap[_val = _1]
                            | ( ( lit( "CommonLogger" ) > mRCommonLogger ) | ( lit( "FileLogger" ) > mRFileLogger ) )[_val =
                                        phoenix::construct<LoggerPtr>( _1 )] );

    mRId = char_( "a-zA-Z" )[_val = _1] >> *char_( "0-9a-zA-Z" )[_val += _1];

    mRIndependentLogger = ( ( lit( "CommonLogger" ) > mRId > mRCommonLogger )[phoenix::bind(
                                &LoggerCreator::addLogger, *this, _1, phoenix::construct<LoggerPtr>( _2 ) )]
                            | ( lit( "FileLogger" ) > mRId > mRFileLogger )[phoenix::bind( &LoggerCreator::addLogger, *this, _1,
                                    phoenix::construct<LoggerPtr>( _2 ) )

                                                                           ] ) > -lit( ';' );

    mRLoggerName.name( "LoggerName" );
    mRFileName.name( "FileName" );
    mRCommonLogger.name( "CommonLogger" );
    mRFileLogger.name( "FileLogger" );
    mRSolverBoundLogger.name( "Logger" );
    mRTimer.name( "Timer" );

    on_error<fail>( mRIndependentLogger, std::cout << val( "Error! Expecting " ) << _4 // what failed?
                    << val( " here: \"" ) << construct<std::string>( _3, _2 ) // iterators to error-pos, end
                    << val( "\" thrown behind: \"" ) << construct<std::string>( _1, _3 ) // iterators to start, error-pos
                    << "\"" << std::endl );

    on_error<fail>( mRSolverBoundLogger, std::cout << val( "Error! Expecting " ) << _4 // what failed?
                    << val( " here: \"" ) << construct<std::string>( _3, _2 ) // iterators to error-pos, end
                    << val( "\" thrown behind: \"" ) << construct<std::string>( _1, _3 ) // iterators to start, error-pos
                    << "\"" << std::endl );
}

LoggerCreator::RuleType& LoggerCreator::getSolverBoundRule()
{
    return getInstance().mRSolverBoundLogger;
}

qi::rule<std::string::const_iterator,void(),ascii::space_type>& LoggerCreator::getIndependentRule()
{
    return getInstance().mRIndependentLogger;
}

LoggerCreator& LoggerCreator::getInstance()
{
    static std::auto_ptr<LoggerCreator> instance;

    if ( !instance.get() )
    {
        instance = std::auto_ptr<LoggerCreator>( new LoggerCreator() );
    }
    return *instance;
}

void LoggerCreator::addLogger( const std::string& name, LoggerPtr loggerPtr )
{
    LAMA_ASSERT( !mLoggerInstanceMap.find( name ), "Logger with id " << name << " already registered." );

    mLoggerInstanceMap.add( name, loggerPtr );

    LAMA_LOG_INFO( logger, "Registered logger " << name << " describing " << loggerPtr->id() );
}

} // namespace lama
