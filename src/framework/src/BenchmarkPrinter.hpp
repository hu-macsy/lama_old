/**
 * @file print.h
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
 * @brief print.h
 * @author Jiri Kraus
 * @date 06.04.2011
 * $Id$
 */
/*
 * print.h
 *
 *  Created on: 31.01.2011
 *      Author: rrehrman
 */

#ifndef LAMA_BENCHMARKPRINTER_HPP_
#define LAMA_BENCHMARKPRINTER_HPP_

#include <iostream>

#include <framework/src/Benchmark.hpp>

namespace bf
{

class BenchmarkPrinter
{

public:

    /** Print message to std::cout */
    template<typename T>
    static inline void print( const T& message );

    /** Print message to std::cerr */
    template<typename T>
    static inline void error( const T& message );

    template<typename T>
    static inline void warning( const T& message );

    static void setDoOutput( const bool doOutput );

private:
    BenchmarkPrinter();
    BenchmarkPrinter( const BenchmarkPrinter& other );
    BenchmarkPrinter& operator=( const BenchmarkPrinter& other );

    static const char* const message_begin;
    static const char* const message_end;
    static const char* const error_begin;
    static const char* const error_end;
    static const char* const warning_begin;
    static const char* const warning_end;

    static bool doOutput;

};

/** Print message to std::cout */
template<typename T>
inline void BenchmarkPrinter::print( const T& message )
{
    if( doOutput )
    {
        std::cout << message_begin << std::endl;
        std::cout << message << std::endl;
        std::cout << message_end << std::endl;
    }
}

/** Print message to std::cerr */
template<typename T>
inline void BenchmarkPrinter::error( const T& message )
{
    if( doOutput )
    {
        std::cerr << error_begin << std::endl;
        std::cerr << message << std::endl;
        std::cerr << error_end << std::endl;
    }
}

template<typename T>
inline void BenchmarkPrinter::warning( const T& message )
{
    if( doOutput )
    {
        std::cerr << warning_begin << std::endl;
        std::cerr << message << std::endl;
        std::cerr << warning_end << std::endl;
    }
}

} // namespace bf

#endif // LAMA_BENCHMARKPRINTER_HPP_
