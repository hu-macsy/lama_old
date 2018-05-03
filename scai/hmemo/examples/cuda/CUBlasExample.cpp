/**
 * @file hmemo/examples/cuda/CUBlasExample.cpp
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
 * @brief Using of CUSparseLibray with HArray and CUDAContext
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#include <scai/hmemo.hpp>

// include of CUDAContext implies include of cuBLAS, cuSPARSE

#include <scai/hmemo/cuda/CUDAContext.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>

#include <iostream>

SCAI_LOG_DEF_LOGGER( logger, "CudaExample" )

using namespace scai;
using namespace hmemo;

template<typename ValueType>
void outArray( const HArray<ValueType>& array, const char* name )
{
    std::cout << name << "[ " << array.size() << " ] = {";
    ContextPtr contextPtr = Context::getContextPtr( common::ContextType::Host );
    ReadAccess<ValueType> read( array, contextPtr );

    for ( IndexType i = 0; i < array.size(); ++i )
    {
        std::cout << " " << read.get()[i];
    }

    std::cout << " }" << std::endl;
}

int main()
{
    ContextPtr cuda = Context::getContextPtr( common::ContextType::CUDA );
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
        const common::CUDACtx& dev = common::CUDAAccess::getCurrentCUDACtx();
        SCAI_CUBLAS_CALL( cublasSdot( dev.getcuBLASHandle(), n,
                                      rA.get(), 1, rB.get(), 1, &dot ),
                          "cublasSDot<float>" );
    }
    std::cout << "dot result: " << dot << std::endl;
}
