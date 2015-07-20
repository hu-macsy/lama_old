/**
 * @file GenLoggerCreator.cpp
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
 * @brief Implemenation of methods for class GenLoggerCreator.
 * @author Thomas Brandes
 * @date 01.03.2011
 */

#include <iostream>

#include <cstdlib>         // import getenv
#include <cstdio>          // FILE
#include <stdexcept>       // runtime_error

#include <logging/logging.hpp>
#include <logging/GenLoggerCreator.hpp>
#include <logging/GenLogger.hpp>

#undef DEBUGGING

using namespace std;

namespace logging
{

// GenLoggerCreator becomes the 'static' global logger creator.

AbstractLoggerCreator& theLoggerCreator()
{
    return GenLoggerCreator::getTheCreator();
}

GenLoggerCreator* GenLoggerCreator::theCreator = NULL;

GenLoggerCreator& GenLoggerCreator::getTheCreator()
{
    if ( !theCreator )
    {
        theCreator = new GenLoggerCreator();
    }

    return *theCreator;
}

GenLoggerCreator::~GenLoggerCreator()
{
}

/********************************************************************
 *  GenLoggerCreator:: getRoot()                                     *
 ********************************************************************/

Logger& GenLoggerCreator::getRoot() const
{
    return GenLogger::getRoot();
}

Logger* GenLoggerCreator::create( const std::string& name, Logger* parent ) const
{
    return new GenLogger( name, parent );
}

} //namespace logging

