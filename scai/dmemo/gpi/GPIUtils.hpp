/**
 * @file GPIUtils.hpp
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
 * @brief Utility macros for GPI calls
 * @author Lauretta Schubert
 * @date 25.02.2014
 * @since 1.1.0
 */

#pragma once

#include <scai/dmemo/gpi/GPIException.hpp>

#include <sstream>
#include <cstdio>
#include <iostream>
#include <GASPI.h>

/** This directive might be enabled for debugging  */

#define SCAI_GASPI_CALL( call )                                           \
    {                                                                     \
        gaspi_return_t status = call;                                     \
                                                                          \
        if ( status != GASPI_SUCCESS )                                    \
        {                                                                 \
            std::ostringstream errorStr;                                  \
            errorStr << "GPI error in line " << __LINE__ ;                \
            errorStr << " of file " << __FILE__ << ": ";                  \
            errorStr << #call << "\n";                                    \
            scai::common::Exception::addCallStack( errorStr );            \
            std::fprintf( stderr, "%s\n", errorStr.str().c_str() );       \
            gaspi_printf( "%s\n", errorStr.str().c_str() );               \
            throw GPIException( errorStr.str(), status );                 \
        }                                                                 \
    }

