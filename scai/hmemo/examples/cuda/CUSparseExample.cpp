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
	extern cusparseHandle_t CUDAContext_cusparseHandle;
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
        ReadAccess<int> readcooIA( cooIA, cuda );
        WriteOnlyAccess<int> writecsrIA( csrIA, cuda, numRows + 1 );

        SCAI_CUSPARSE_CALL(
            cusparseXcoo2csr( CUDAContext_cusparseHandle,
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
        ReadAccess<int> readcsrIA( csrIA, cuda );
        ReadAccess<int> readcsrJA( csrJA, cuda );
        ReadAccess<float> readcsrValues( csrValues, cuda );
  
        WriteOnlyAccess<int> writecscJA( cscJA, cuda, numColumns + 1 );
        WriteOnlyAccess<int> writecscIA( cscIA, cuda, numValues );
        WriteOnlyAccess<float> writecscValues( cscValues, cuda, numValues );
     
        SCAI_CUSPARSE_CALL(
            cusparseScsr2csc( CUDAContext_cusparseHandle,
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
