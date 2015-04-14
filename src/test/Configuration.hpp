/**
 * @file Configuration.hpp
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
 * @brief Configuration.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * @since 1.0.0
 */
#ifndef LAMA_CONFIGURATION_HPP_
#define LAMA_CONFIGURATION_HPP_

#include <string>
#include <lama/CommunicatorFactory.hpp>

class Configuration
{

public:
    virtual ~Configuration();
    static Configuration& getInstance();
    const std::string& getPath() const;
    const std::string& getCommunicatorType() const;
    void setCommunicatorType( const std::string& commType );

private:
    void setPath( const std::string& path );

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    Configuration    ();
    Configuration( const Configuration& cc );
    std::string mPath;
    std::string mCommType;

};

#endif // LAMA_CONFIGURATION_HPP_
