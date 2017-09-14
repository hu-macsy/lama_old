/**
 * @file BenchmarkPrinter.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief print.h
 * @author Jiri Kraus
 * @date 06.04.2011
 */
/*
 * print.h
 *
 *  Created on: 31.01.2011
 *      Author: rrehrman
 */

#pragma once

#include <iostream>

#include <scai/benchmark/Benchmark.hpp>

namespace scai
{

namespace bf
{

/** Static class for printing messages within benchmarks.
 *
 *  Output statements are encapsulated with special patterns so that these
 *  messages can be extracted later by the BenchmarkRunner.
 */

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

} // namespace scai
