/**
 * @file examples/cuda/DemoCUSparse.cpp
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
 * @brief ToDo: Missing description in ./examples/cuda/DemoCUSparse.cpp
 * @author Thomas Brandes
 * @date 14.03.2016
 */

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>

#include <iostream>
#include <memory>

/* --------------------------------------------------------------------- */

using namespace scai;
using namespace common;

/* --------------------------------------------------------------------- */

template<typename T>
T* myAllocate( int N )
{
    // allocate memory on the accessed device and copy host data to it
    CUdeviceptr d_pointer = 0;
    size_t size = sizeof( T ) * N;
    // allocate memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &d_pointer, sizeof( T ) * N ), "cuMemAlloc( size = " << size << " ) failed." )
    // initialize data with 0
    SCAI_CUDA_DRV_CALL( cuMemsetD8( d_pointer, 0, size ), "init" )
    return reinterpret_cast<T*>( d_pointer );
}

/* --------------------------------------------------------------------- */

template<typename T>
T* myAllocate( const T h_data[], int N )
{
    // allocate memory on the accessed device and copy host data to it
    CUdeviceptr pointer = 0;
    size_t size = sizeof( T ) * N;
    // allocate memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, sizeof( T ) * N ), "cuMemAlloc( size = " << size << " ) failed." )
    // transfer host data
    SCAI_CUDA_DRV_CALL( cuMemcpyHtoD( pointer, h_data, size ), "tranfer host->device" )
    return reinterpret_cast<T*>( pointer );
}

/* --------------------------------------------------------------------- */

template<typename T>
void myFree( const T* d_data )
{
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( d_data );
    SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << d_data << " ) failed" )
}


template<typename T>
void myFree( T h_data[], const T* d_data , int N )
{
    CUdeviceptr d_pointer = reinterpret_cast<CUdeviceptr>( d_data );
    size_t size = sizeof( T ) * N;
    // transfer data from device to host
    SCAI_CUDA_DRV_CALL( cuMemcpyDtoH( h_data, d_pointer, size  ), "tranfer device->host" )
    SCAI_CUDA_DRV_CALL( cuMemFree( d_pointer ), "cuMemFree( " << d_data << " ) failed" )
}

/* --------------------------------------------------------------------- */

void myCOO2CSR( int* d_csr_ia, const int* d_coo_ia, int numValues, int numRows )
{
    const CUDACtx& device = CUDAAccess::getCurrentCUDACtx();
    std::cout << "cusparseXcoo2csr, #values = " << numValues << ", #rows = " << numRows << std::endl;
    SCAI_CUSPARSE_CALL( cusparseXcoo2csr( device.getcuSparseHandle(),
                                          d_coo_ia, numValues, numRows, d_csr_ia,
                                          CUSPARSE_INDEX_BASE_ZERO ),
                        "coo2csr" )
}

/* --------------------------------------------------------------------- */

void myCSR2CSC( int* d_csc_ia, int* d_csc_ja, float* d_csc_values,
                const int* d_csr_ia, const int* d_csr_ja, const float* d_csr_values,
                int numRows, int numColumns, int numValues )
{
    const CUDACtx& device = CUDAAccess::getCurrentCUDACtx();
    std::cout << "cusparseScsr2csc, size = " << numRows << " x " << numColumns << ", #values = " << numValues << std::endl;
    SCAI_CUSPARSE_CALL( cusparseScsr2csc( device.getcuSparseHandle(),
                                          numRows, numColumns, numValues,
                                          d_csr_values, d_csr_ia, d_csr_ja,
                                          d_csc_values, d_csc_ia, d_csc_ja,
                                          CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO ),
                        "csr2csc" )
}

/* --------------------------------------------------------------------- */

int main( int argc, const char** argv )
{
    // at least --SCAI_DEVICE=id may be specified
    Settings::parseArgs( argc, argv );
    int nr = 0;   // take this as default
    Settings::getEnvironment( nr, "SCAI_DEVICE" );
    CUDACtx device( nr );
    /***********************************************************************
     *  Definition of input data                                           *
     *                                                                     *
     *     1.0  2.0   -    -                                               *
     *     5.0  3.0  4.0   -                                               *
     *     -     -   5.0  6.0                                              *
     *     -    7.0   -   8.0                                              *
     *                                                                     *
     **********************************************************************/
    int numRows = 4;
    int numColumns = 4;
    int ia[] = { 0, 0, 1, 1, 1, 2, 2, 3, 3 };
    int ja[] = { 0, 1, 0, 1, 2, 2, 3, 1, 3 };
    float values[] = { 1.0, 2.0, 5.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
    int numValues  = sizeof( ia ) / sizeof( int );
    SCAI_ASSERT_EQUAL( numValues, sizeof( ja ) / sizeof( int ), "COO size mismatch" )
    SCAI_ASSERT_EQUAL( numValues, sizeof( values ) / sizeof( float ), "COO size mismatch" )

    for ( int i = 0; i < numRows; ++i )
    {
        std::cout << "row[" << i << "] =";

        for ( int k = 0; k < numValues; ++k )
        {
            if ( ia[k] != i )
            {
                continue;
            }

            std::cout << " " << ja[k] << ":" << values[k];
        }

        std::cout << std::endl;
    }

    // not verified here: max( ia ) < numRows, max( ja ) < numColumns
    CUDAAccess access( device );
    int* d_coo_ia = myAllocate( ia, numValues );
    int* d_coo_ja = myAllocate( ja, numValues );
    float* d_coo_values = myAllocate( values, numValues );
    int* d_csr_ia = myAllocate<int>( numRows + 1 );
    int* d_csr_ja = d_coo_ja;
    float* d_csr_values = d_coo_values;
    myCOO2CSR( d_csr_ia, d_coo_ia, numValues, numRows );
    int* d_csc_ia = myAllocate<int>( numValues );
    int* d_csc_ja = myAllocate<int>( numColumns + 1 );
    float* d_csc_values = myAllocate<float>( numValues );
    myCSR2CSC( d_csc_ia, d_csc_ja, d_csc_values,
               d_csr_ia, d_coo_ja, d_coo_values,
               numRows, numColumns, numValues );
    // free memory on device and get host data
    std::unique_ptr<int[]> csc_ia( new int[numValues] );
    std::unique_ptr<int[]> csc_ja( new int[numColumns + 1] );
    std::unique_ptr<float[]> csc_values( new float[numValues] );
    myFree( d_coo_ia );
    myFree( d_csr_values );
    myFree( d_csr_ia );
    myFree( d_csr_ja );
    myFree( csc_ia.get(), d_csc_ia, numValues );
    myFree( csc_ja.get(), d_csc_ja, numColumns + 1 );
    myFree( csc_values.get(), d_csc_values, numValues );
    /***********************************************************************
     *   CSC ja:  0    2    5     7    9                                   *
     *                                                                     *
     *            1.0  2.0   -    -                                        *
     *            5.0  3.0  4.0   -                                        *
     *             -     -  5.0  6.0                                       *
     *             -   7.0   -   8.0                                       *
     *                                                                     *
     **********************************************************************/
    int   expected_ja[]     = { 0,    2,       5,    7,    9 };
    int   expected_ia[]     = { 0, 1, 0, 1, 3, 1, 2, 2, 3 };
    float expected_values[] = { 1, 5, 2, 3, 7, 4, 5, 6, 8 };

    for ( int j = 0; j < numColumns; ++j )
    {
        std::cout << "column[" << j << "] =";

        for ( int k = csc_ja[j]; k < csc_ja[j + 1]; ++k )
        {
            std::cout << " " << csc_ia[k] << ":" << csc_values[k];
        }

        std::cout << std::endl;
    }

    for ( int i = 0; i < numColumns + 1; ++i )
    {
        SCAI_ASSERT_EQUAL( expected_ja[i], csc_ja[i], "ja mismatch for i = " << i )
    }

    for ( int i = 0; i < numValues; ++i )
    {
        SCAI_ASSERT_EQUAL( expected_ia[i], csc_ia[i], "ia mismatch at index " << i )
    }

    for ( int i = 0; i < numValues; ++i )
    {
        SCAI_ASSERT_EQUAL( expected_values[i], csc_values[i], "values mismatch at index " << i )
    }
}
