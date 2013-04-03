/**
 * @file DenseMatrixOps.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Implementation of matrix operations for dense matrices.
 * @author Thomas Brandes
 * @date 04.01.2012
 * $Id$
 */

// hpp
#include <lama/matrix/DenseMatrixOps.hpp>

// others
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/HostWriteAccess.hpp>

#include <lama/distribution/NoDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>

// tracing
#include <lama/tracing.hpp>

// boost
#include <boost/scoped_array.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( DenseMatrixOps::logger, "DenseMatrixOps" )

template<typename T>
void DenseMatrixOps::invertReplicated( DenseMatrix<T>& matrix )
{
    LAMA_REGION( "DenseMatrixOps::invertReplicated" )
    LAMA_ASSERT_EQUAL_ERROR( 1, matrix.getDistribution().getNumPartitions() )
    LAMA_ASSERT_EQUAL_ERROR( 1, matrix.getColDistribution().getNumPartitions() )

    DenseStorage<T>& denseStorage = matrix.getLocalStorage();

    denseStorage.invert( denseStorage );
}

/* ------------------------------------------------------------------------ */

template<typename T>
void DenseMatrixOps::invertCyclic( DenseMatrix<T>& matrix )
{
    LAMA_REGION( "DenseMatrixOps::invertCyclic" )

    const Communicator& comm = matrix.getDistribution().getCommunicator();

    const CyclicDistribution* cyclicDist = dynamic_cast<const CyclicDistribution*>( matrix.getDistributionPtr().get() );

    LAMA_ASSERT_ERROR( cyclicDist, "no cyclic distribution: " << matrix.getDistribution() )

    const int nb = cyclicDist->chunkSize(); // blocking factor

    ContextPtr context = matrix.getContextPtr();

    const LAMAInterface& lamaInterface = context->getInterface();

    if ( !lamaInterface.getSCALAPACKInterface<T>().inverse )
    {
        LAMA_THROWEXCEPTION( "No SCALAPACK routines available" )
    }

    const int n = matrix.getNumRows();

    // assert square matrix

    LAMA_ASSERT_EQUAL_ERROR( matrix.getNumColumns(), n )

    DenseStorage<T>& denseStorage = matrix.getLocalStorage();

    const IndexType localSize = denseStorage.getData().size();

    LAMA_ASSERT_EQUAL_ERROR( localSize, denseStorage.getNumRows() * n )

    LAMA_LOG_INFO( logger, "local dense data = " << denseStorage << ", localSize = " << localSize )

    HostWriteAccess<T> localValues( denseStorage.getData() );

    T* data = localValues.get();

    LAMA_LOG_INFO( logger, "now call inverse" )

    lamaInterface.getSCALAPACKInterface<T>().inverse( n, nb, data, comm );
}

template<typename T>
void DenseMatrixOps::invert( DenseMatrix<T>& matrix )
{
    LAMA_REGION( "DenseMatrixOps::invert" )

    if ( matrix.getDistribution().getNumPartitions() == 1 && matrix.getColDistribution().getNumPartitions() == 1 )
    {
        LAMA_LOG_INFO( logger, "invert called for replicated matrices" )
        invertReplicated( matrix );
    }
    else
    {
        LAMA_LOG_INFO( logger, "invert called for distributed matrices" )
        const CyclicDistribution* cyclicDist =
            dynamic_cast<const CyclicDistribution*>( matrix.getDistributionPtr().get() );
        if ( cyclicDist && matrix.getColDistribution().getNumPartitions() == 1 )
        {
            LAMA_LOG_INFO( logger, "matrix already with cyclic distribution" )
            invertCyclic( matrix );
        }
        else
        {
            const IndexType blockSize = 64;
            DistributionPtr oldRowDist = matrix.getDistributionPtr();
            DistributionPtr oldColDist = matrix.getColDistributionPtr();

            DistributionPtr noDist( new NoDistribution( matrix.getNumRows() ) );
            DistributionPtr cycDist(
                new CyclicDistribution( matrix.getNumRows(), blockSize,
                                        oldRowDist->getCommunicatorPtr() ) );

            LAMA_LOG_INFO( logger, "redistibuting matrix" )
            matrix.redistribute( cycDist, noDist );

            LAMA_LOG_INFO( logger, "calling invert cyclic" )
            invertCyclic( matrix );

            LAMA_LOG_INFO( logger, "back-redistribution" )
            matrix.redistribute( oldRowDist, oldColDist );
        }
    }
}

template LAMA_DLL_IMPORTEXPORT
void DenseMatrixOps::invertCyclic( DenseMatrix<float>& matrix );

template LAMA_DLL_IMPORTEXPORT
void DenseMatrixOps::invertCyclic( DenseMatrix<double>& matrix );

template LAMA_DLL_IMPORTEXPORT
void DenseMatrixOps::invert( DenseMatrix<float>& matrix );

template LAMA_DLL_IMPORTEXPORT
void DenseMatrixOps::invert( DenseMatrix<double>& matrix );

template LAMA_DLL_IMPORTEXPORT
void DenseMatrixOps::invertReplicated( DenseMatrix<float>& matrix );

template LAMA_DLL_IMPORTEXPORT
void DenseMatrixOps::invertReplicated( DenseMatrix<double>& matrix );

} // namespace lama
