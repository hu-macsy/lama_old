/*
 * @file scai/tracing/test/simpleTracing.cpp
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
 * @brief simple executable that is used by the tracing tests
 * @author Jan Ecker
 * @date 18.02.2016
 * @since 2.0.0
*/

#include <scai/tracing.hpp>
#include <scai/logging.hpp>
#include <cstdio>
#include <omp.h>

void subA( int& X )
{
    SCAI_REGION( "A" )
    ++X;
}

void subB( int& X )
{
    SCAI_REGION( "B" )
    X++;
}

int main()
{
    int X = 0;

    SCAI_REGION( "main" )

    #pragma omp parallel for
    for ( int i = 0; i < 10000; ++i )
    {
#ifndef UNNAMED_THREADS
        SCAI_LOG_THREAD( omp_get_thread_num() );
#endif
        SCAI_REGION_START( "main.loopA" )
        for ( int j = 0; j < 30; ++ j )
        {
            subA( X );
        }
        SCAI_REGION_END( "main.loopA" )

        SCAI_REGION_START( "main.loopB" )
        for ( int j = 0; j < 20; ++ j )
        {
            subB( X );
        }
        SCAI_REGION_END( "main.loopB" )
    }
}
