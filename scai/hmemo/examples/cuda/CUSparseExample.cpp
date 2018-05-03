/**
 * @file CUSparseExample.cpp
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
     *  Definition of input data via LAMA arrays                           *
     **********************************************************************/
    int numRows = 4;
    int numColumns = 4;
    int numValues  = 9;
    int ia[] = { 0, 0, 1, 1, 1, 2, 2, 3, 3 };
    int ja[] = { 0, 1, 0, 1, 2, 2, 3, 1, 3 };
    float values[] = { 1.0, 2.0, 5.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
    // by using HArrayRef copy of elements on Host is not required
    HArrayRef<int> cooIA( numValues, ia );
    HArrayRef<int> csrJA( numValues, ja );
    HArrayRef<float> csrValues( numValues, values );
    outArray( cooIA, "cooIA" );
    outArray( csrJA, "cooJA" );
    outArray( csrValues, "csrValues" );
    /***********************************************************************
     *  Conversion COO -> CSR                                              *
     **********************************************************************/
    HArray<int> csrIA;
    {
        SCAI_CONTEXT_ACCESS( cuda )
        ReadAccess<int> readcooIA( cooIA, cuda );
        WriteOnlyAccess<int> writecsrIA( csrIA, cuda, numRows + 1 );
        const common::CUDACtx& dev = common::CUDAAccess::getCurrentCUDACtx();
        SCAI_CUSPARSE_CALL(
            cusparseXcoo2csr( dev.getcuSparseHandle(),
                              readcooIA.get(), numValues, numRows,
                              writecsrIA.get(),
                              CUSPARSE_INDEX_BASE_ZERO ),
            "coo2csr" )
    }
    std::cout << "Conversion COO2CSR done." << std::endl;
    outArray( csrIA, "csrIA" );
    /***********************************************************************
     *  Conversion CSR -> CSC                                              *
     **********************************************************************/
    HArray<int> cscIA;
    HArray<int> cscJA;
    HArray<float> cscValues;
    {
        SCAI_CONTEXT_ACCESS( cuda )
        ReadAccess<int> readcsrIA( csrIA, cuda );
        ReadAccess<int> readcsrJA( csrJA, cuda );
        ReadAccess<float> readcsrValues( csrValues, cuda );
        WriteOnlyAccess<int> writecscJA( cscJA, cuda, numColumns + 1 );
        WriteOnlyAccess<int> writecscIA( cscIA, cuda, numValues );
        WriteOnlyAccess<float> writecscValues( cscValues, cuda, numValues );
        const common::CUDACtx& dev = common::CUDAAccess::getCurrentCUDACtx();
        SCAI_CUSPARSE_CALL(
            cusparseScsr2csc( dev.getcuSparseHandle(),
                              numRows, numColumns, numValues,
                              readcsrValues.get(), readcsrIA.get(), readcsrJA.get(),
                              writecscValues.get(), writecscIA.get(), writecscJA.get(),
                              CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO ),
            "convertCSR2SCC<float>" )
    }
    std::cout << "Conversion CSR2CSC done." << std::endl;
    outArray( cscIA, "cscIA" );
    outArray( cscJA, "cscJA" );
    outArray( cscValues, "cscValues" );
}
