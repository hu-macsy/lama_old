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
 * @brief Demo program for the Factory of HArray.
 * @author: Thomas Brandes, Lauretta Schubert
 * @date 18.04.2012
 **/

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <scai/logging.hpp>

#include <scai/common/macros/throw.hpp>

#include <scai/common/shared_ptr.hpp>

#include <iostream> 

using namespace std;
using namespace scai::hmemo;

SCAI_LOG_DEF_LOGGER( logger, "CreateTest" )

// Template instantiation of LAMArray

template class HArray<double>;

int main()
{
    SCAI_LOG_THREAD( "Main" )

    ContextPtr contextPtr = Context::getHostPtr();

    static IndexType N =  100;

    HArray<float> lamaArray ( N, 1.0 );

    scai::common::shared_ptr<HArray<float> > lamaArray1( HArray<float>::create( lamaArray.getValueType() ) );

    *lamaArray1 = lamaArray;

    ReadAccess<float> read( lamaArray, contextPtr );
    ReadAccess<float> read1( *lamaArray1, contextPtr );
   
    const float* data = read.get();
    const float* data1 = read1.get();

    for ( IndexType i = 0; i < N; ++i )
    {
        SCAI_ASSERT_EQUAL( data[i], data1[i], "" )
    }
  
    std::cout << "Create finished" << std::endl;

    scai::common::shared_ptr<_HArray> lamaArray2( _HArray::create( scai::common::scalar::FLOAT ) );

    std::cout << "lamaArray2 = " << *lamaArray2 << std::endl;

    scai::common::shared_ptr<_HArray> lamaArray3( _HArray::create( scai::common::scalar::DOUBLE ) );

    std::cout << "lamaArray3 = " << *lamaArray3 << std::endl;
}
