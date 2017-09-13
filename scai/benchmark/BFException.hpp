/**
 * @file BFException.h
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
 * @brief BFException.h
 * @author jiri
 * @date 06.04.2011
 * $Id$
 */
/**
 * @file BFException.h
 * @author jiri
 * Created on: 11.05.2010
 */
#pragma once

#include <scai/common/config.hpp>
#include <scai/benchmark/BFError.hpp>

/**
 * @brief The namespace bf holds everything of the benchmark framework
 */
namespace bf
{

class COMMON_DLL_IMPORTEXPORT BFException: public BFError
{
public:
    /**
     * @brief Default constructor.
     */
    BFException();
    /**
     * @brief Constructor initializes BFException with given message.
     * @param[in] message The message of the Error.
     */
    BFException( const std::string& message );
    /**
     * @brief Destructor.
     */
    virtual ~BFException() throw ();
};

}

#ifndef LAMA_CHECK_BENCHMARK
/**
 * Checks an expression and prints out a warning at compile time, if
 * expression was not true.
 */
#define LAMA_CHECK_BENCHMARK(exp)                                              \
    if (!(exp))                                                                     \
    {                                                                               \
        std::ostringstream errorStr;                                                \
        errorStr<<"Warning: Incorrect Results.\n";                                  \
        throw bf::BFException(errorStr.str( ));                                \
    }

#endif
