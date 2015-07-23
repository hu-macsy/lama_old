/**
 * @file Create.cpp
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
 * @brief Demo program for the Factory of LAMAArray.
 * @author: Thomas Brandes, Lauretta Schubert
 * @date 18.04.2012
 **/

#include <memory/LAMAArray.hpp>
#include <memory/WriteAccess.hpp>
#include <memory/ReadAccess.hpp>

#include <logging/logging.hpp>

#include <common/Exception.hpp>

#include <common/shared_ptr.hpp>

#include <iostream> 

using namespace std;
using namespace memory;
using namespace common;

LAMA_LOG_DEF_LOGGER( logger, "CreateTest" )

// Template instantiation of LAMArray

template class LAMAArray<double>;

int main()
{
    LAMA_LOG_THREAD( "Main" )

    static IndexType N =  100;

    LAMAArray<float> lamaArray ( N, 1.0 );

    shared_ptr<LAMAArray<float> > lamaArray1( lamaArray.clone() );

    *lamaArray1 = lamaArray;

    ReadAccess<float> read( lamaArray );
    ReadAccess<float> read1( *lamaArray1 );
   
    const float* data = read.get();
    const float* data1 = read1.get();

    for ( IndexType i = 0; i < N; ++i )
    {
        COMMON_ASSERT_EQUAL( data[i], data1[i], "" )
    }
  
    std::cout << "Create finished" << std::endl;

    std::cout << "LAMAArray<float>::initialized = " << LAMAArray<float>::initialized << std::endl;

    shared_ptr<ContextArray> lamaArray2( ContextArray::create( scalar::FLOAT ) );

    std::cout << "lamaArray2 = " << *lamaArray2 << std::endl;

    shared_ptr<ContextArray> lamaArray3( ContextArray::create( scalar::DOUBLE ) );

    std::cout << "lamaArray3 = " << *lamaArray3 << std::endl;
}
