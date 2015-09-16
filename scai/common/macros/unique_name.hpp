/**
 * @file unique_name.hpp
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
 * @brief Some macro utilities for generating symbol names.
 * @author Jiri Kraus
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once

/** Help macro to concatenate two symbols, can also deal with nested calls. */

#define LAMA_JOIN( symbol1, symbol2 ) _LAMA_DO_JOIN( symbol1, symbol2 )

/** Help macro to deal with nested calls of LAMA_JOIN */

#define _LAMA_DO_JOIN( symbol1, symbol2 ) _LAMA_DO_JOIN2( symbol1, symbol2 )

/** Furtherhelp macro to deal with nested calls of LAMA_JOIN */

#define _LAMA_DO_JOIN2( symbol1, symbol2 ) symbol1##symbol2

/** @brief Creates a unique symbol name by joining the prefix, the line and the postfix.
 *
 *  \code
 *  LAMA_UNIQUE_NAME( Interface, Registry ) -> Interface17Registry
 *  \endcode
 */

#define LAMA_UNIQUE_NAME( prefix, postfix )                                    \
    LAMA_JOIN( prefix , LAMA_JOIN( __LINE__ , postfix ) )
