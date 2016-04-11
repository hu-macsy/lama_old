/**
 * @file OpenMP.cpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief Checking the OpenMP version
 * @author Eric Schricker
 * @date 07.04.2016
 * @since 2.0.0
 */

#include <stdio.h>
#include <omp.h>

#ifndef _OPENMP
    #define _OPENMP -1
#endif

int main(int argc, char *argv[])
{
    float version;

    switch( _OPENMP )
    {
        case 200505:
            version = 2.5; break;
        case 200805:
            version = 3.0; break;
        case 201107:
            version = 3.1; break;
        case 201307:
            version = 4.0; break;
        case 201511:
            version = 4.5; break;
        default:
            fprintf( stderr, "Unknown OpenMP-Version\n" );
            return -1;
    }
    
    printf("%3.1f", version );
}
