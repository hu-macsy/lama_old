/**
 * @file BFError.h
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
 * @brief BFError.h
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
 */
/*
 * BFError.h
 *
 *  Created on: 01.02.2011
 *      Author: rrehrman
 */

#pragma once

#include <string>
#include <exception>

#include <scai/common/config.hpp>

namespace bf
{

class COMMON_DLL_IMPORTEXPORT BFError: public std::exception
{
public:
    /**
     * @brief The default constructor to create an Exception without message.
     */
    BFError();
    /**
     * @brief The constructor to create an Exception with the given message.
     * @param[in] message   The message of the Exception.
     */
    BFError( const std::string& message );
    /**
     * @brief The destructor of the Exception.
     */
    virtual ~BFError() throw ();

    /**
     * @brief Returns the message of the Exception.
     * @return The message of the Exception.
     */
    virtual const char* what() const throw ();
private:
    /** The message of the Exception. */
    std::string m_message;
};

} // namespace bf
