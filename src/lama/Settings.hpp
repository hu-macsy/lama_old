/**
 * @file LAMASettings.hpp
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
 * @brief Managing some settings for CUDA specified by environment variables
 * @author Thomas Brandes
 * @date 04.05.2013
 * @since 1.0.0
 */
#ifndef LAMA_LAMA_SETTINGS_HPP_
#define LAMA_LAMA_SETTINGS_HPP_

#include <logging/logging.hpp>

namespace lama
{

/** 
 *  This class provides only static methods.
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
    static bool getEnvironmentSetting( bool& flag, const char* envVarName );

    /** Set a integer by value of its environment variable
     *
     *  @param[out]  val is integer variable that will be set
     *  @param[in]   envVarName is name of the environment variable
     *  @return      true if environment variable has been used to set flag
     */
    static bool getEnvironmentSetting( int& flag, const char* envVarName );

private:

    Settings();

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    /** convert the string value to a boolean value, name only used for messages. 
     *
     *  @param[out]  flag is boolean variable that will be set
     *  @param[in]   value is string to be converted
     *  @return      true if string could be converted, false if no legal value has been found
     */

    static bool convertYesNoString( bool& flag, const char* value );

    /** convert the string value to an int value
     *
     *  @param[out]  int is variable that will be set
     *  @param[in]   value is string to be converted
     *  @return      true if string could be converted, false if no legal value has been found
     */

    static bool convertValue( int& flag, const char* value );

    static bool init();
};

} // namespace

#endif //  LAMA_LAMA_SETTINGS_HPP_
