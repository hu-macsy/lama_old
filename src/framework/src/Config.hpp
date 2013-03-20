/**
 * @file Config.h
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
 * @brief Config.h
 * @author: robin
 * @date 06.04.2011
 * $Id$
 */
/*
 * @file Config.h
 * @author: robin
 * Created on: 03.08.2010
 */

#ifndef LAMA_CONFIG_HPP_
#define LAMA_CONFIG_HPP_

#include <string>
#include <map>

#include <framework/src/config_framework.hpp>

namespace bf
{

/**
 * @brief This class holds static parameters for the benchmarks, which can be
 *        set by a main program, and asked by each benchmark, in case it is not
 *        possible to give this values as parameter values to the benchmarks.
 *        For example, this class could hold the path to input data.
 *        This class is a singleton, so there does exist only one instance of
 *        this class.
 */
class LAMABENCHFRAME_DLL_IMPORTEXPORT Config
{
public:
    virtual ~Config();
    /**
     * @brief returns the onliest instance of this class.
     *
     * @return the onliest instance of this class.
     */
    static Config& getInstance();

    /**
     * @brief returns the value of the given parameter.
     *
     * @param[in] param the parameter, to get the value of.
     *
     * @return A constant reference to the value of the parameter.
     */
    const std::string& getValueOf( const std::string& param );

    /**
     * @brief sets the given value for the given parameter.
     *
     * Sets the given value for the given parameter, in case the given
     * parameter does not exist, yet.
     *
     * @param[in]   param   the parameter to be set.
     * @param[in]   value   the value of the parameter.
     */
    void setValueFor( std::string param, std::string value );

private:
    Config();
    Config( const Config& cc );
    std::map<std::string,std::string> m_params;
};

} // namespace bf

#endif // LAMA_CONFIG_HPP_
