/**
 * @file Configuration.cpp
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
 * @brief Configuration.cpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * @since 1.0.0
 */
#include <test/Configuration.hpp>

#include <logging/logging.hpp>

#ifndef LAMA_TESTFILE_PATH
#define LAMA_TESTFILE_PATH "res/testfiles"
#endif

LAMA_LOG_DEF_LOGGER( Configuration::logger, "Configuration" );

Configuration::Configuration()
    : mPath( LAMA_TESTFILE_PATH ), mCommType( "none" )
{
}

Configuration::~Configuration()
{
}

Configuration& Configuration::getInstance()
{
    static Configuration cfg;
    return cfg;
}

const std::string& Configuration::getPath() const
{
    return mPath;
}

void Configuration::setPath( const std::string& path )
{
    LAMA_LOG_DEBUG( logger, "path = " << path );
    mPath = path;
}

const std::string& Configuration::getCommunicatorType() const
{
    return mCommType;
}

void Configuration::setCommunicatorType( const std::string& commType )
{
    mCommType = commType;
}
