/**
 * @file Settings.hpp
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
 * @brief Utilities for LAMA settings
 * @author Thomas Brandes
 * @date 19.06.2013
 * @since 1.0.1
 */
#pragma once

#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

class Communicator;

/**
 *  This singleton class provides methods to query environment variables.
 *
 *  Note: This should be the only module to access environment variables
 *        directly.
 */

class Settings
{
public:

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
     *  @return      true if environment variable has been used to set flag
     */
    static bool getEnvironment( std::string& val, const char* envVarName );

    /** Set a string by value of its environment variable
     *
     *  @param[out]  val is string that will be set
     *  @param[in]   envVarName is name of the environment variable
     *  @param[in]   comm is communicator used to bcast setting of environment variable
     *  @return      true if environment variable has been used to set flag
     */
    static bool getEnvironment( std::string& val, const char* envVarName, const Communicator& comm );

private:

    Settings();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** convert the string value to a boolean value, name only used for messages.
     *
     *  @param[out]  flag is boolean variable that will be set
     *  @param[in]   value is string to be converted
     *  @return      true if string could be converted, false if no legal value has been found
     */

static    bool convertYesNoString( bool& flag, const char* value );

    /** convert the string value to an int value
     *
     *  @param[out]  number is variable that will be set
     *  @param[in]   value is string to be converted
     *  @return      true if string could be converted, false if no legal value has been found
     */

    static bool convertValue( int& number, const char* value );

    static bool init();
};

} /* end namespace lama */

} /* end namespace scai */
