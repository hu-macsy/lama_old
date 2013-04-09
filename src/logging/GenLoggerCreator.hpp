/**
 * @file GenLoggerCreator.hpp
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
 * @brief GenLoggerCreator.hpp
 * @author brandes
 * @date 01.03.2011
 * $Id$
 */

#ifndef LAMA_GEN_LOGGER_CREATOR_HPP_
#define LAMA_GEN_LOGGER_CREATOR_HPP_

// base classes
#include <logging/AbstractLoggerCreator.hpp>

namespace log4lama
{

class GenLoggerCreator: public AbstractLoggerCreator
{

public:

    virtual ~GenLoggerCreator();

    /**
     * @brief Function to create a new instance at a given level.
     * @param[in] name identifies the instance must not contain any dots
     * @param parent is existent parent logger
     * @return Pointer to the new created instance
     *
     * Note: the new created instance will be added as a new son to parent.
     */

    virtual Logger* create( const std::string& name, Logger* parent ) const;

    /** Getter for the root instance of all loggers. */

    virtual Logger& getRoot() const;

    /** Get the single instantiation of this class */

    static GenLoggerCreator& getTheCreator();

private:

    /** Only one instance of this object is generated.
     *
     */
    GenLoggerCreator()
    {
    }

    static GenLoggerCreator* theCreator; //!< singleton instantiation
};

}

#endif // LAMA_GEN_LOGGER_CREATOR_HPP_
