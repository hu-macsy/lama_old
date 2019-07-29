/**
 * @file Settings.hpp
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
 * @brief Utilities for LAMA settings
 * @author Thomas Brandes
 * @date 19.01.2016
 */
#pragma once

#include <scai/common/config.hpp>

// std

#include <string>
#include <vector>
#include <sstream>

/** Namespace used for all projects of the LAMA software and other provided packages like logging, tracing, etc. */

namespace scai
{

namespace common
{

/**
 *  This singleton class provides methods to query settings of environment variables and
 *  command line arguments.
 *
 *  Note: This should be the only module to access environment variables directly
 *        as this operation is OS specific.
 */

class COMMON_DLL_IMPORTEXPORT Settings
{
public:

    /** Parse the command line arguments.
     *
     *  @param[in,out] argc number of command line arguments
     *  @param[in,out] argv arguments
     *
     *  The number of arguments and the arguments are reduced by those starting with --SCAI_...
     */
    static void parseArgs( int& argc, const char* argv[] );

    /** Set a value of its environment variable
     *
     *  @param[out]  val is variable that will be set
     *  @param[in]   envVarName is name of the environment variable
     *  @return      true if environment variable has been used to set flag
     */
    template<typename ValueType>
    static bool getEnvironment( ValueType& val, const char* envVarName );

    /** Set a string by value of its environment variable
     *
     *  @param[out]  val is string that will be set
     *  @param[in]   envVarName is name of the environment variable
     *  @return      true if environment variable was set and provided an integer vlaue
     */
    static bool getRankedEnvironment( std::string& val, const char* envVarName );

    /** Define an environment variable */

    static void putEnvironment( const char* envVarName, const char* val, bool replace = true );

    /** Define an environment variable */

    static void putEnvironment( const char* envVarName, const int val, bool replace = true );

    /** Get tokenized string from an environment variable
     *
     *  @param[out] vals is a vector of separated strings from the environment varialbe
     *  @param[in]  envVarName is name of the environment variable
     *  @param[in]  delimiters contains characters used for splitting
     *  @return     true if environment variable was set
     */
    static bool getEnvironment( std::vector<std::string>& vals, const char* envVarName, const char* delimiters );

    /** Help routine to tokenize a string by a given separator
     *
     *  @param[out] tokens is a vector of separated strings from the input string
     *  @param[in]  input is a string that will be tokenized
     *  @param[in]  delimiters contains all characters used for separation
     *
     */
    static void tokenize( std::vector<std::string>& tokens, const std::string& input, const std::string& delimiters = " " );

    /** This method sets globally the number of argument that is taken by a comma separated value list. */

    static void setRank( int rank );

    /** Print all environment variables starting with SCAI_ (only for debug, demo purpose) in an output stream */

    static void printEnvironment( std::ostream& out );

    /**
     *  @brief Read environment variables from an input file.
     * 
     *  @param[in] fileName name of the input file
     *  @param[in] name     string argument to find line in input file with specific settings
     *  @param[in] rank     additional int argument to find line in input file
     *  @returns            number of environment variables that have been set, -1 if no entry was found
     *
     *  \code
     *  # Example file with environment variables
     *  ...
     *  <name> <rank> DOMAIN=3 SCAI_CONTEXT=Host SCAI_NP=2x2 WEIGHT=4
     *  ...
     */
    static int readSettingsFile( const char* fileName, const char* name, int rank );

private:

    Settings();

    /** convert the string value to a boolean value, name only used for messages.
     *
     *  @param[out]  flag is boolean variable that will be set
     *  @param[in]   value is string to be converted
     *  @return      true if string could be converted, false if no legal value has been found
     */

    static bool convertYesNoString( bool& flag, const char* value );

    static int sRank;  //<!  specifies pos to take from comma separated values

    static const char* RANK_DELIMITER();

    /**
     *  @brief Function to set environment variable via string like "<var>=<value>"
     *
     *  The method throws an exception if the input string does not contain an assignment operator.
     */
    static void setEntry ( const std::string& setting, const char* fileName, int line );
};

/* ----------------------------------------------------------------------------- */

template<typename ValueType>
inline bool Settings::getEnvironment( ValueType& val, const char* envVarName )
{
    std::string envVal;

    if ( !Settings::getRankedEnvironment( envVal, envVarName ) )
    {   
        return false; // no initialization by environment
    }

    std::istringstream input( envVal );

    input >> val;

    return input.fail() == 0;
}

template<>
inline bool Settings::getEnvironment( std::string& val, const char* envVarName )
{
    // workaround with input stream not needed here

    return Settings::getRankedEnvironment( val, envVarName );
}

template<>
inline bool Settings::getEnvironment( bool& flag, const char* envVarName )
{
    std::string val;

    if ( !getEnvironment( val, envVarName ) )
    {
        return false;  // environment variable not set
    }

    return convertYesNoString( flag, val.c_str() );
}

} /* end namespace common */

} /* end namespace scai */
