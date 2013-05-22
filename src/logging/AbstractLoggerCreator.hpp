/**
 * @file AbstractLoggerCreator.hpp
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
 * @brief Abstract class for logger creator.
 * @author Thomas Brandes
 * @date 02.03.2011
 * $Id$
 */
#ifndef LAMA_ABSTRACTLOGGERCREATER_HPP_
#define LAMA_ABSTRACTLOGGERCREATER_HPP_

#include <logging/Logger.hpp>

namespace log4lama
{

/** Abstract class from which classes must derive that can create loggers. */

class AbstractLoggerCreator
{

public:
    /**
     * @brief Destructor needed due to virtual functions
     */
    virtual ~AbstractLoggerCreator()
    {
    }

    /**
     * @brief get the root logger
     */

    virtual class Logger& getRoot() const = 0;

    /**
     * @brief Function to create a new instance at a given level.
     * @param[in] name is identification of logger, must not contain any dots
     * @param[in,out] parent is the parent logger, new logger is set as son
     */

    virtual class Logger* create( const std::string& name, class Logger* parent ) const =0;
};

// The implementation of this static method decides which logger
// creator will be used for the static loggers

AbstractLoggerCreator& theLoggerCreator();

}

#endif // LAMA_ABSTRACTLOGGERCREATER_HPP_
