/**
 * @file Settings.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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

    /** Set a flag by value of its environment variable
     *
     *  @param[out]  flag is boolean variable that will be set
     *  @param[in]   envVarName is name of the environment variable
     *  @return      true if environment variable has been used to set flag
     */
    static bool getEnvironment( bool& flag, const char* envVarName );

    /** Set a integer by value of its environment variable
     *
     *  @param[out]  val is integer variable that will be set
     *  @param[in]   envVarName is name of the environment variable
     *  @return      true if environment variable has been used to set flag
     */
    static bool getEnvironment( int& val, const char* envVarName );

    /** Set a string by value of its environment variable
     *
     *  @param[out]  val is string that will be set
     *  @param[in]   envVarName is name of the environment variable
     *  @return      true if environment variable was set and provided an integer vlaue
     */
    static bool getEnvironment( std::string& val, const char* envVarName );

    /** Define an environment variable */

    static void putEnvironment( const char* envVarName, const char* val, bool replace = true );

    /** Define an environment variable */

    static void putEnvironment( const char* envVarName, const int val, bool replace = true );

    /** Get tokenized string from an environment variable 
     *
     *  @param[out] vals is a vector of separated strings from the environment varialbe
     *  @param[in]  envVarName is name of the environment variable
     *  @param[in]  separator is the character used to separate
     *  @return     true if environment variable was set
     */
    static bool getEnvironment( std::vector<std::string>& vals, const char* envVarName, const char separator );

    /** Help routine to tokenize a string by a given separator
     *
     *  @param[out] values is a vector of separated strings from the input string
     *  @param[in]  input is a string that will be tokenized
     *  @param[in]  seperator is the character used to separate
     *
     */
    static void tokenize( std::vector<std::string>& values, const std::string& input, const char seperator );

    /** This method sets globally the number of argument that is taken by a comma separated value list. */

    static void setRank( int rank );

    /** Print all environment variables starting with SCAI_ (only for debug, demo purpose) */

    static void printEnvironment();

private:

    Settings();

    /** convert the string value to a boolean value, name only used for messages.
     *
     *  @param[out]  flag is boolean variable that will be set
     *  @param[in]   value is string to be converted
     *  @return      true if string could be converted, false if no legal value has been found
     */

    static bool convertYesNoString( bool& flag, const char* value );

    /** convert the string value to an int value
     *
     *  @param[out]  number is variable that will be set
     *  @param[in]   value is string to be converted
     *  @return      true if string could be converted, false if no legal value has been found
     */

    static bool convertValue( int& number, const char* value );

    static int sRank;  //<!  specifies pos to take from comma separated values

    static const char RANK_SEPARATOR_CHAR = ',';
};

} /* end namespace common */

} /* end namespace scai */
