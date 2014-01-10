/**
 * @file cudaComputeCapability.c
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
 * @brief BiCGstabTest.cpp
 * @author Jan Ecker
 * @date 10.01.2014
 * @since 1.1.0
 *
 * Based on code by Florian Rathgeber<florian.rathgeber@gmail.com> on
 * http://www.cmake.org/Bug/print_bug_page.php?bug_id=11767
 */

#include <stdio.h>
#include <cuda_runtime.h>
#include <stdlib.h>

int main()
{
    int deviceCount;
    struct cudaDeviceProp properties;
    long int cuda_device;

    if ( cudaGetDeviceCount( &deviceCount ) != cudaSuccess )
        return 1;

    if ( getenv( "LAMA_DEVICE" ) ) {
        char *pEnd;
        cuda_device = strtol( getenv( "LAMA_DEVICE" ), &pEnd, 10 );
    } else {
        // TODO:  search for device with highest compute capability
        cuda_device = 0;
    }

    if ( deviceCount > cuda_device ) {
        cudaGetDeviceProperties( &properties, cuda_device );
        printf( "%d%d", properties.major, properties.minor );
        return 0;
    }

    // failure
    return 1;
}
