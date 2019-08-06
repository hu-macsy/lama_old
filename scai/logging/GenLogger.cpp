/**
 * @file GenLogger.cpp
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
 * @brief Implementation of methods for generic logger.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

// local library
#include <scai/logging.hpp>

// internal scai libraries
#include <scai/common/exception/Exception.hpp>
#include <scai/common/Settings.hpp>

// std
#include <iostream>
#include <sstream>

#include <cstdio>          // FILE
#include <stdexcept>       // runtime_error
#include <cstring>

#undef DEBUGGING

using namespace std;

namespace scai
{

using common::Settings;

namespace logging
{

// default is not to flush

bool GenLogger::sFlush = false;

std::vector<std::string> GenLogger::formatTokens;

int ( *GenLogger::myPrintf ) ( const char* format, ... ) = ( int (* ) ( const char* format, ... ) )& printf ;

/********************************************************************
 *  Static variable: rootLogger for generic logging                  *
 ********************************************************************/

GenLogger* GenLogger::rootLogger = NULL;

/********************************************************************
 *  GenLogger:: getRoot()                                            *
 ********************************************************************/

Logger& GenLogger::getRoot()
{
    if ( rootLogger == NULL )
    {
        rootLogger = new GenLogger( "<root>", NULL );
#ifdef DEBUGGING
        printf( "root logger now available, do configure\n" );
#endif
        // Note: call configure after Root Logger is available
        GenLogger::configure();
    }

    return *rootLogger;
}

/********************************************************************
 *  GenLogger : Constructor                                          *
 ********************************************************************/

GenLogger::GenLogger( const std::string& name, Logger* parent )
    : Logger( name, parent )

{
#ifdef DEBUGGING
    printf( "GenLogger %s created\n", getFullName().c_str() );
#endif
}

/********************************************************************
 *  help routine: string2bool                                        *
 ********************************************************************/

static bool string2bool( const std::string& value )
{
    if ( value == "1" )
    {
        return true;
    }

    if ( value == "0" )
    {
        return false;
    }

    if ( value == "TRUE" )
    {
        return true;
    }

    if ( value == "FALSE" )
    {
        return false;
    }

    if ( value == "ON" )
    {
        return true;
    }

    if ( value == "OFF" )
    {
        return false;
    }

    throw std::runtime_error( "illegal boolean value" );
}

/********************************************************************
 *  help routine: evalEntry                                          *
 ********************************************************************/

static int evalEntry( char* line, int length, const char* /* filename */ )
{
    line[length] = '\0';
    string myLine = line;
    // check for an empty line
    string::size_type firstPos = myLine.find_first_not_of( " ", 0 );
#ifdef DEBUGGING
    printf( "pos of first relevant char = %lu\n", firstPos );
#endif

    if ( string::npos == firstPos )
    {
        return 0;
    }

    // check for a comment line
#ifdef DEBUGGING
    printf( "first relevant char = %c\n", myLine[firstPos] );
#endif

    if ( myLine[firstPos] == '#' )
    {
        return 0;
    }

    // check for an equal sign in the line
    string::size_type equalPos = myLine.find_first_of( "=", 0 );

    if ( string::npos == equalPos )
    {
        throw std::runtime_error( "no equal sign" );
    }

    // now find the name without blanks, e.g. "atom vec = " is only atom
    string::size_type lastPos = myLine.find_first_of( " =", firstPos );
    string name = myLine.substr( firstPos, lastPos - firstPos );
    firstPos = myLine.find_first_not_of( " ", equalPos + 1 );

    if ( string::npos == firstPos )
    {
        throw std::runtime_error( "no value after = " );
    }

    if ( myLine[firstPos] == '"' )
    {
        // look for matching " and do not
        firstPos++;
        lastPos = myLine.find_first_of( "\"", firstPos );

        if ( string::npos == lastPos )
        {
            lastPos = myLine.length();
        }
    }
    else
    {
        lastPos = myLine.find_first_of( " ", firstPos );

        if ( string::npos == lastPos )
        {
            lastPos = myLine.length();
        }

        for ( string::size_type i = firstPos; i <= lastPos; i++ )
        {
            myLine[i] = static_cast<string::value_type>( toupper( myLine[i] ) );
        }
    }

    string value = myLine.substr( firstPos, lastPos - firstPos );

    //Check for other options
    if ( name == "flush" )
    {
        GenLogger::setFlush( string2bool( value ) );
        return 1;
    }

    if ( name == "format" )
    {
        GenLogger::setFormat( value );
        return 1;
    }

    // take entries of SCAI_xxx as environment variables

    if ( strncmp( name.c_str(), "SCAI_", 5 ) == 0 )
    {
        // this is not a logging entry so take it as environment
        bool replace = false;   // do not override existing settings
        Settings::putEnvironment( name.c_str(), value.c_str(), replace );
        return 1;
    }

    // get the logger from the provider and set its level
    level::Level level = str2level( value );

    if ( level == level::MAXLEVEL )
    {
        throw std::runtime_error( "illegal log level" );
    }

    if ( name == "<root>" )
    {
        //Empty name references the rootLogger
        name = "";
    }

    LoggerProvider::getProvider().getInstance( name ).setLevel( level );
    return 1;
}

/********************************************************************
 *  Read logger configuration from a file
 ********************************************************************/

#define MAX_LINE_LENGTH 256

int GenLogger::readConfig( const char* fname )
{
    char buffer[MAX_LINE_LENGTH];
    FILE* configFile = fopen( fname, "r" ); // file with logger configuration

    if ( configFile == NULL )
    {
        SCAI_LOG_ERROR( ( *rootLogger ), "config: could not open config file " << fname );
        return 0;
    }

    int bufferLength = 0;
    char eof = EOF;
    int noEntries = 0; // number of relevant entries
    bool stop = false; // becomes true for termination of read loop

    while ( !stop )
    {
        char c = static_cast<char>( fgetc( configFile ) );

        if ( c == '\n' || c == eof )
        {
            /* new line, evaluate current line */
            try
            {
                noEntries += evalEntry( buffer, bufferLength, fname );
            }
            catch ( std::runtime_error& e )
            {
                SCAI_LOG_WARN( ( *rootLogger ),
                               "Config file '" << fname << "', ignored invalid line '" << buffer << "'" << ": " << e.what() );
            }

            bufferLength = 0;
            stop = ( c == eof );
        }
        else
        {
            buffer[bufferLength++] = c;

            if ( bufferLength == MAX_LINE_LENGTH )
            {
                SCAI_LOG_ERROR( ( *rootLogger ), "Config file '" << fname << "', too long line, stop reading" );
                stop = true;
            }
        }
    }

    fclose( configFile );
    return noEntries;
}

/********************************************************************
 *  Configuration of GenLogger done by reading configuration file    *
 ********************************************************************/

void GenLogger::configure()
{
    if ( !rootLogger )
    {
        throw std::runtime_error( "configure: rootLogger not available yet" );
    }

    rootLogger->setLevel( level::WARN ); // default setting for root

    // Set default format string, now as it might be used very soon

    if ( formatTokens.size() == 0 )
    {
        setFormat( "#date, #time #name @ #thread ( #func -> #file::#line ) #level #msg" );
    }

    std::string configFile;
    bool logDefined = Settings::getEnvironment( configFile, "SCAI_LOG" );

    if ( !logDefined )
    {
        // environment variable SCAI_LOG not set, so we try it at $HOME/.loggingrc
        if ( Settings::getEnvironment( configFile, "HOME" ) )
        {
            configFile += "/.loggingrc";
            FILE* fp = fopen ( configFile.c_str(), "r" );

            if ( fp != NULL )
            {
                fclose( fp );
                logDefined = true;   // file exists, so we take this as SCAI_LOG specification
            }
        }
    }

    if ( !logDefined )
    {
        SCAI_LOG_WARN( *rootLogger, "SCAI_LOG not set, no $HOME/.loggingrc, so use default configuration" );
        configFile.clear();
    }

    if ( configFile.length() == 0 )
    {
        rootLogger->setLevel( level::WARN );
    }
    else if ( configFile == level2str( level::OFF ) )
    {
        rootLogger->setLevel( level::OFF );
    }
    else if ( configFile == level2str( level::FATAL ) )
    {
        rootLogger->setLevel( level::FATAL );
    }
    else if ( configFile == level2str( level::SERROR ) )
    {
        rootLogger->setLevel( level::SERROR );
    }
    else if ( configFile == level2str( level::WARN ) )
    {
        rootLogger->setLevel( level::WARN );
    }
    else if ( configFile == level2str( level::INFO ) )
    {
        rootLogger->setLevel( level::INFO );
    }
    else if ( configFile == level2str( level::DEBUG ) )
    {
        rootLogger->setLevel( level::DEBUG );
    }
    else if ( configFile == level2str( level::TRACE ) )
    {
        rootLogger->setLevel( level::TRACE );
    }
    else
    {
        SCAI_LOG_INFO( ( *rootLogger ), "read configuration from file " << configFile );
        readConfig( configFile.c_str() );
    }

    rootLogger->traverse(); // traverse all loggers and might be print it
}

/********************************************************************
 *  Helper routine for logging via Python                            *
 ********************************************************************/

static void writeVal2( std::ostringstream& stream, int val )
{
    if ( val <= 9 )
    {
        stream << "0";
    }

    stream << val;
}

static void writeTime( std::ostringstream& stream )
{
    time_t timer;
    time ( &timer );
    struct tm* tp = localtime( &timer );
    writeVal2( stream, tp->tm_hour );
    stream << ":";
    writeVal2( stream, tp->tm_min );
    stream << ":";
    writeVal2( stream, tp->tm_sec );
}

static void writeDate( std::ostringstream& stream )
{
    time_t timer;
    time ( &timer );
    struct tm* tp = localtime( &timer );
    stream << ( 1900 + tp->tm_year ) << "-";
    writeVal2( stream, 1 + tp->tm_mon );
    stream << "-";
    writeVal2( stream, tp->tm_mday );
}

/********************************************************************
 *  general log routine                                              *
 ********************************************************************/

void GenLogger::log( const char* level, SourceLocation& loc, const string& msg )
{
    std::ostringstream output;

    for ( size_t i = 0; i < formatTokens.size(); ++i )
    {
        const std::string& token = formatTokens[i];

        if ( token[0] != '#' )
        {
            output << formatTokens[i];
        }
        else if ( token == "#NAME" )
        {
            output << getFullName();
        }
        else if ( token == "#TIME" )
        {
            writeTime( output );
        }
        else if ( token == "#DATE" )
        {
            writeDate( output );
        }
        else if ( formatTokens[i] == "#THREAD" )
        {
            output << *common::thread::getCurrentThreadName();
        }
        else if ( formatTokens[i] == "#FILE" )
        {
            output << loc.mFileName;
        }
        else if ( formatTokens[i] == "#LINE" )
        {
            output << loc.mLine;
        }
        else if ( formatTokens[i] == "#FUNC" )
        {
            output << loc.mFuncName;
        }
        else if ( formatTokens[i] == "#LEVEL" )
        {
            output << level;
        }
        else if ( formatTokens[i] == "#MSG" )
        {
            output << msg;
        }
        else if ( formatTokens[i] == "#STACK" )
        {
            // undocumented feature: print stack
            scai::common::Exception::addCallStack( output );
        }
        else
        {
            // ignore first character # and take it as environment variable
            const char* var = formatTokens[i].c_str() + 1;
            std::string val;

            if ( scai::common::Settings::getEnvironment( val, var ) )
            {
                output << val;
            }
            else
            {
                output << "${" << var << "}";
            }
        }
    }

    output << std::endl;

    myPrintf( "%s", output.str().c_str() );

    if ( sFlush )
    {
#ifdef DEBUGGING
        printf( "Flushed\n" );
#endif
        fflush( stdout );
    }
}

/********************************************************************
 *  Implementation of trace/debug/info/warn/error/fatal              *
 ********************************************************************/

void GenLogger::trace( SourceLocation loc, const string& msg )
{
    log( "TRACE", loc, msg );
}

void GenLogger::debug( SourceLocation loc, const string& msg )
{
    log( "DEBUG", loc, msg );
}

void GenLogger::info( SourceLocation loc, const string& msg )
{
    log( "INFO", loc, msg );
}

void GenLogger::warn( SourceLocation loc, const string& msg )
{
    log( "WARN", loc, msg );
}

void GenLogger::error( SourceLocation loc, const string& msg )
{
    log( "ERROR", loc, msg );
}

void GenLogger::fatal( SourceLocation loc, const string& msg )
{
    log( "FATAL", loc, msg );
}

/********************************************************************
 *  GenLogger::traverse()                                            *
 ********************************************************************/

void GenLogger::traverse()
{
    // SCAI_LOG_DEBUG(myLogger, "rootLogger " << getFullName() << ", level = "
    //            << getEffectiveLevel() << ", set = " << setFlag);
    for ( size_t i = 0; i < mSons.size(); i++ )
    {
        GenLogger* son = dynamic_cast<GenLogger*>( mSons[i] );

        if ( son )
        {
            son->traverse();
        }
    }
}

void GenLogger::setFlush( bool flush )
{
    GenLogger::sFlush = flush;
}

/********************************************************************
 *  GenLogger::setFormat( formatString )                             *
 ********************************************************************/

static void logTokenize( std::vector<std::string>& tokens, const std::string& input )
{
    tokens.clear();
    std::string::size_type lastPos = 0;
    std::string::size_type pos     = input.find_first_of( "#", lastPos );

    while ( std::string::npos != pos )
    {
        // found
        if ( lastPos < pos )
        {
            tokens.push_back( input.substr( lastPos, pos - lastPos ) );
        }

        lastPos = input.find_first_not_of( "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVXYZ_", pos + 1 );
        tokens.push_back( input.substr( pos, lastPos - pos ) );
        pos = input.find_first_of( "#", lastPos );
    }

    if ( lastPos < pos )
    {
        tokens.push_back( input.substr( lastPos, pos - lastPos ) );
    }
}

void GenLogger::setFormat( const std::string& format )
{
    logTokenize( formatTokens, format );

    // convert all tokens to upper case

    for ( size_t i = 0; i < formatTokens.size(); ++i )
    {
        std::string& val = formatTokens[i];

        if ( val.length() > 0  && val[0] == '#' )
        {
            for ( std::string::iterator p = val.begin(); val.end() != p; ++p )
            {
                *p = static_cast<string::value_type>( toupper( *p ) );
            }
        }
    }
}

} /* end namespace logging */

} /* end namespace scai */
