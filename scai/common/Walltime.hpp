/**
 * @file common/Walltime.hpp
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
 * @brief Definition of static class that provides high precision walltime
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// std
#include <stdint.h>

namespace scai
{

namespace common
{

/** Make a common defintion for 8-byte integer as uint64_t can be implementation specific. */

typedef uint64_t INTEGER_8;

/**
 * @brief A simple static class that delivers walltime (used for logging and tracing)
 */
class COMMON_DLL_IMPORTEXPORT Walltime
{
public:

    /** Get the current walltime.
     *
     *  @return current walltime in ticks
     */
    static INTEGER_8 timestamp();

    /** Number of ticks per second, so timestamp() / timerate() gives time in seconds. */

    static INTEGER_8 timerate();

    /** Get the current walltime.
     *
     *  @return current walltime in seconds
     */
    static double get();

private:

    /** Private constructor for a static class. */

    Walltime();
};

} /* end namespace common */

} /* end namespace scai */
