/**
 * @file GenLogger.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Implementation of methods for generic logger.
 * @author Thomas Brandes
 * @date 01.03.2011
 * @since 1.0.0
 */
#include <iostream>

#include <logging.hpp>
#include <cstdlib>         // import getenv
#include <cstdio>          // FILE
#include <stdexcept>       // runtime_error
#include <cstring>

// hpp
#include <logging/GenLogger.hpp>

#undef DEBUGGING

using namespace std;

namespace log4lama
{

// default is not to flush

bool GenLogger::sFlush = false;

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
    printf("GenLogger %s created\n", getFullName().c_str());
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

static int evalEntry( char* line, int length, const char* /* filename */)
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
    lastPos = myLine.find_first_of( " ", firstPos );
#ifdef DEBUGGING
    printf( "value at %lu - %lu\n", firstPos, lastPos );
#endif

    if ( string::npos == lastPos )
    {
        lastPos = myLine.length();
    }

    for ( string::size_type i = firstPos; i <= lastPos; i++ )
    {
        myLine[i] = static_cast<string::value_type>( toupper( myLine[i] ) );
    }

    string value = myLine.substr( firstPos, lastPos - firstPos );

    //Check for other options
    if ( name == "flush" )
    {
        GenLogger::setFlush( string2bool( value ) );
        return 1;
    }

    // get the logger from the provider and set its level
    Level level = str2level( value );

    if ( level == MAXLEVEL )
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
        LAMA_LOG_ERROR( ( *rootLogger ), "config: could not open config file " << fname );
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
                LAMA_LOG_WARN( ( *rootLogger ),
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
                LAMA_LOG_ERROR( ( *rootLogger ), "Config file '" << fname << "', too long line, stop reading" );

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

    rootLogger->setLevel( WARN ); // default setting for root

#ifdef LAMA_LOG
    const char* configFile = LAMA_LOG;
#else
    const char* configFile = getenv( "LAMA_LOG" );
#endif

    if ( configFile == NULL )
    {
        LAMA_LOG_WARN( ( *rootLogger ), "LAMA_LOG not set, use default configuration" );
    }
    else if ( strlen( configFile ) == 0 )
    {
        rootLogger->setLevel( WARN );
    }
    else if ( strcmp( configFile, level2str( OFF ) ) == 0 )
    {
        rootLogger->setLevel( OFF );
    }
    else if ( strcmp( configFile, level2str( FATAL ) ) == 0 )
    {
        rootLogger->setLevel( FATAL );
    }
    else if ( strcmp( configFile, level2str( SERROR ) ) == 0 )
    {
        rootLogger->setLevel( SERROR );
    }
    else if ( strcmp( configFile, level2str( WARN ) ) == 0 )
    {
        rootLogger->setLevel( WARN );
    }
    else if ( strcmp( configFile, level2str( INFO ) ) == 0 )
    {
        rootLogger->setLevel( INFO );
    }
    else if ( strcmp( configFile, level2str( DEBUG ) ) == 0 )
    {
        rootLogger->setLevel( DEBUG );
    }
    else if ( strcmp( configFile, level2str( TRACE ) ) == 0 )
    {
        rootLogger->setLevel( TRACE );
    }
    else
    {
        LAMA_LOG_INFO( ( *rootLogger ), "read configuration from file " << configFile );
        readConfig( configFile );
    }

    rootLogger->traverse(); // traverse all loggers and might be print it
}

/********************************************************************
 *  Helper routine for logging via Python                            *
 ********************************************************************/

void GenLogger::log( const char* level, SourceLocation& loc, const string& msg )
{
    printf( "%s (%s::%d,func=%s) %s: %s\n", getFullName().c_str(), loc.mFileName, loc.mLine, loc.mFuncName, level,
            msg.c_str() );

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
    // LAMA_LOG_DEBUG(myLogger, "rootLogger " << getFullName() << ", level = "
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

} //namespace log4lama

