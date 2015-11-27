/**
 * @file CUSparseExample.cpp
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
 * @brief Using of CUSparseLibray with HArray and CUDAContext
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#include <scai/hmemo.hpp>

// include of CUDAContext implies include of cuBLAS, cuSPARSE

#include <scai/hmemo/cuda/CUDAContext.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <iostream>

SCAI_LOG_DEF_LOGGER( logger, "CudaExample" )

namespace scai
{
	extern cublasHandle_t CUDAContext_cublasHandle;

} /* end namespace scai */

using namespace scai;
using namespace hmemo;

template<typename ValueType>
void outArray( const HArray<ValueType>& array, const char* name )
{
    std::cout << name << "[ " << array.size() << " ] = {";

    ContextPtr contextPtr = Context::getContextPtr( common::context::Host );

    ReadAccess<ValueType> read( array, contextPtr );

    for ( IndexType i = 0; i < array.size(); ++i )
    {
        std::cout << " " << read.get()[i];
    }
    std::cout << " }" << std::endl;
}

int main()
{
    ContextPtr cuda = Context::getContextPtr( common::context::CUDA );

    /***********************************************************************
     *  Definition of input data via heterogeneous arrays                  *
     **********************************************************************/

    float a[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    float b[] = { 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 };

    int n  = sizeof( a ) / sizeof( float );
    int n1 = sizeof( b ) / sizeof( float );

    SCAI_ASSERT_EQUAL( n, n1, "mismatch of arrays for dot product" )

    // by using HArrayRef copy of elements on Host is not required

    HArrayRef<float> hA( n, a );
    HArrayRef<float> hB( n, b );

    outArray( hA, "hA" );
    outArray( hB, "hB" );

    /***********************************************************************
     *  Dotproduct                                                         *
     **********************************************************************/

    float dot;

    {
        ReadAccess<float> rA( hA, cuda );
        ReadAccess<float> rB( hB, cuda );

        SCAI_CONTEXT_ACCESS( cuda )

        SCAI_CUBLAS_CALL( cublasSdot( CUDAContext_cublasHandle, n,
                                      rA.get(), 1, rB.get(), 1, &dot),
                                      "cublasSDot<float>" );
  
    }

    std::cout << "dot result: " << dot << std::endl;
}
