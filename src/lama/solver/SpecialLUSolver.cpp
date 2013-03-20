/**
 * @file SpecialLUSolver.cpp
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
 * @brief SpecialLUSolver.cpp
 * @author Robin Rehrmann
 * @date 18.07.2011
 * $Id$
 */

// hpp
#include <lama/solver/SpecialLUSolver.hpp>

// others
#include <lama/solver/h_frag.h>

#include <lama/CommunicatorFactory.hpp>
#include <lama/DenseStorage.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>
#include <lama/LAMAInterface.hpp>

#include <lama/OpenMP/OpenMPBLAS1.hpp>
#include <lama/OpenMP/OpenMPBLAS2.hpp>
#include <lama/OpenMP/OpenMPLAPACK.hpp>

#include <lama/distribution/CyclicDistribution.hpp>

#ifndef LAMA_BUILD_CUDA
#else
#include <lama/CUDA/CUDAContext.hpp>
#include <lama/CUDA/CUDAHostContext.hpp>
#include <lama/CUDA/CUDAHostContextManager.hpp>
#endif

#include <math.h>
#include <omp.h>

namespace lama
{

void computeROOTforwards( lama::PartitionId& root, const lama::IndexType& numProcs )
{
    root = ( root + 1 ) % numProcs;
}

void computeROOTbackwards( lama::PartitionId& root, const lama::IndexType& numProcs )
{
    if( root == 0 )
    {
        root = numProcs - 1;
    }
    else
    {
        --root;
    }
}

bool lowerFin( lama::IndexType& i, lama::PartitionId& root, const lama::IndexType& nP, const lama::IndexType& bH )
{
    computeROOTforwards( root, nP );
    ++i;
    return i < bH;
}

bool upperFin( lama::IndexType& i, lama::PartitionId& root, const lama::IndexType& nP, const lama::IndexType& )
{
    computeROOTbackwards( root, nP );
    --i;
    return i >= 0;
}

void lowerInit(
    lama::IndexType& begin,
    lama::IndexType& end,
    const lama::IndexType& i,
    const lama::IndexType& nP,
    const lama::IndexType& bW,
    const lama::PartitionId& me,
    const lama::PartitionId& root )
{
    begin = i / nP + ( me <= root ? 1 : 0 );
    end = bW;
}

void upperInit(
    lama::IndexType& begin,
    lama::IndexType& end,
    const lama::IndexType& i,
    const lama::IndexType& nP,
    const lama::IndexType&,
    const lama::PartitionId& me,
    const lama::PartitionId& root )
{
    begin = 0;
    end = ( i + 1 ) / nP + ( i + 1 ) % nP - ( me >= root ? 1 : 0 );
}

} // namespace

namespace lama
{

LAMA_LOG_DEF_LOGGER( LUSolver::logger, "LUSolver" );

const double LUSolver::epsilon = 1 + 1E-06;

LUSolver::LUSolver( const std::string& id )
    : Solver( id ), mTilesize( 64 ), mDev( -1 ), mComm()
{
}

LUSolver::LUSolver( const std::string& id, LoggerPtrlogger )
    : Solver( id, logger ), mTilesize( 64 ), mDev( -1 ), mComm()
{
}

LUSolver::LUSolver( const LUSolver& other )
    : Solver( other ), mTilesize( other.mTilesize ), mDev( other.mDev ), mComm()
{
}

LUSolver::LUSolverRuntime::LUSolverRuntime()
    : SolverRuntime(), mLUfactorization( 0 ), mPermutation( 0 )
{
}

LUSolver::~LUSolver()
{
}

LUSolver::LUSolverRuntime::~LUSolverRuntime()
{
    if( mLUfactorization != 0 )
    {
        delete mLUfactorization;
    }
}

void LUSolver::factorMatrixToLU( Matrix& matrix, std::vector<IndexType>& permutation )
{
    if( typeid( matrix ) == typeid(DenseMatrix<double> ) )
    {
        computeLUFactorization( dynamic_cast<DenseMatrix<double>&>( matrix ), permutation );
    }
    else if( typeid( matrix ) == typeid(DenseMatrix<float> ) )
    {
        computeLUFactorization( dynamic_cast<DenseMatrix<float>&>( matrix ), permutation );
    }
    else
    {
        LAMA_THROWEXCEPTION( "LU factorization cannot be done for matrix "<<matrix );
    }
}

inline void LUSolver::initializeCommunicator()
{
    if( mComm )
    {
        return;
    }
    mComm = CommunicatorFactory::get( "MPI" );
}

void LUSolver::initialize( const Matrix& coefficients )
{
    LAMA_REGION("initialize");
    LUSolverRuntime& runtime = getRuntime();

    if( runtime.mLUfactorization != 0 )
    {
        delete runtime.mLUfactorization;
    }

    runtime.mPermutation.resize( coefficients.getNumRows(), 0 );

    if( matrixIs<double>( coefficients ) )
    {
        DenseMatrix<double>* dLUfactorization = new DenseMatrix<double>( coefficients );

        computeLUFactorization( *dLUfactorization, runtime.mPermutation );

        runtime.mLUfactorization = dLUfactorization;
    }
    else if( matrixIs<float>( coefficients ) )
    {
        DenseMatrix<float>* sLUfactorization = new DenseMatrix<float>( coefficients );

        computeLUFactorization( *sLUfactorization, runtime.mPermutation );

        runtime.mLUfactorization = sLUfactorization;
    }

    Solver::initialize( coefficients );
}

template<typename T>
void LUSolver::computeLUFactorization( DenseMatrix<T> & matrix, std::vector<IndexType> & permutation )
{
    LAMA_REGION("computeLUFactorization");

    typedef T ValueType;
    const ValueType minus = static_cast<ValueType>( -1.0 );
    const ValueType one = static_cast<ValueType>( 1.0 );
    const IndexType inc1 = 1;

    if( matrix.getDistribution().getNumPartitions() == 1 && matrix.getColDistribution().getNumPartitions() == 1 )
    {
        const IndexType n = matrix.getNumRows();
        const IndexType m = matrix.getNumColumns();

        if( permutation.size() != static_cast<std::vector<IndexType>::size_type>( n ) )
        {
            permutation.resize( n, 0 );
        }

        const IndexType tilesize = computeTilesize( n, m );
        LAMA_LOG_DEBUG( logger, "Size of tiles is " << tilesize );

        boost::scoped_array<IndexType> ipiv( new IndexType[n] );
        for( IndexType i = 0; i < n; ++i )
        {
            permutation[i] = ipiv[i] = i;
        }

        HostWriteAccess<T> decomposition( matrix.getLocalValues().getData() );
        const IndexType tileEnd = std::min( m, n );

        /* In this loop, the LU-factorization of  the matrix is computed,  using BLAS and LAPACK routines. The algorithm
         * is blocked. There do exist five parts within the matrix, as described by the draft.
         *
         *                                                    M
         *                                +---------------------------------------+
         *                                |                                       |
         *                                |    A L R E A D Y   C O M P U T E D    |
         *                                |                                       |
         *                                |                                       |
         *                                + - - - - - - - - + - - - - + - - - - - +
         *                                | LASWP on the    |         |  T R S M  |
         *                                | left side                             |
         *                                |                 |  GETRF  |    ^      |    N
         *                                |                              LASWP    |
         *                                |                 |         + - -|- - - +
         *                                |                              right    |
         *                                |                 |         |    v      |
         *                                |                                       |
         *                                |                 |         |  G E M M  |
         *                                |                                       |
         *                                +-----------------+---------+-----------+
         *
         * The main block is  the diagonal block. It  begins at  position (i,i) and  has the size (n-i)x(tilesize). This
         * block is being computed  by the LAPACK  routine GETRF. That routine computes the complete LU-factorization of
         * this block and generates a permutation vector, by which it has swapped the rows in its block.
         *  Besides the GETRF-block  there is a leften side  and a righten side. The left side has already been computed
         * by the GETRF routine in previous steps. The onyl thing to do here, is to swap the rows, according to the just
         * generated  permutation vector.  This is done by the  LAPACK routine LASWP. The implementation of this routine
         * tough may call the BLAS routine SWAP,  because LASWP assumes  the given matrix  being stored in  column major
         * order. Translating the matrix in place before and after calling the routine LASWP, to satisfy this assumtion,
         * is too expensiv, so if the matrix was stored in row  major order, the rows are swapped using the BLAS level 1
         * routine SWAP.
         *  The right side of the GETRF block needs to swap its rows, too, also using the LAPACK routine LASWP.
         * Afterwards, it split into an upper and a lower block.
         *  The upper block begins at (i,i+tilesize) and is of size  (tilesize)x(m-i-tilesize). It uses the BLAS level 3
         * routine TRSM to solve the equation AX=B  where A is the lower  tridiagonal matrix computed by the GETRF block
         * and B is the matrix being computed by the TRSM block..
         *  The lower block begins at (i+tilesize,i+tilesize)  and is of  size (n-i-tilesize)x(m-i-tilesize). By calling
         * the BLAS level 3 routine GEMM it computes the difference of lower block of the GETRF block and the AXPY block
         * from its own.
         *  The last block of the  matrix is  the block  above the GETRF block. This  block has already been computed by
         * previous steps and does not need to be touched during the next steps, again.
         */
        for( IndexType diagonal = 0; diagonal < tileEnd; diagonal += tilesize )
        {
            const IndexType tile_right_side = ( diagonal + tilesize < m ? tilesize : m - diagonal );
            const IndexType tile_bottom = ( diagonal + tilesize < n ? tilesize : n - diagonal );

            OpenMPLAPACK::getrf( CblasRowMajor, n - diagonal, tile_right_side, &decomposition[diagonal * m + diagonal],
                                 n, &ipiv.get()[diagonal] );

            for( IndexType i = 0; i < tile_bottom; ++i )
            {
                if( ipiv[diagonal + i] != i )
                {
                    std::swap( permutation[diagonal + i], permutation[ipiv[diagonal + i] + diagonal] );
                }
            }

            // swap on the right side, then on the left side of the getrf-area.
            // splits up the omp threads for load-balancing
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    OpenMPLAPACK::laswp( CblasRowMajor, m - diagonal - tile_right_side,
                    &decomposition[diagonal * m + diagonal + tile_right_side], m, 0, tile_bottom,
                    ipiv.get() + diagonal, inc1 );
                }
                #pragma omp section
                {
                    OpenMPLAPACK::laswp( CblasRowMajor, diagonal, &decomposition[diagonal * m], m, 0, tile_bottom,
                    ipiv.get() + diagonal, inc1 );
                }
            }

            const T* const A = &decomposition[diagonal * n + diagonal];

            // compute the upper right side.
            #pragma omp parallel for schedule (LAMA_OMP_SCHEDULE)
            for( IndexType offset = diagonal + tilesize; offset < m; offset += tilesize )
            {
                const IndexType& bottom = tile_bottom;
                const IndexType right_side = ( offset + tilesize < m ? tilesize : m - offset );

                OpenMPBLAS3::trsm( CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, bottom, right_side,
                                   one, A, n, &decomposition[diagonal * n + offset], n );
            }

            // compute the lower right side.
            #pragma omp parallel for schedule (LAMA_OMP_SCHEDULE)
            for( IndexType iTile = diagonal + tilesize; iTile < m; iTile += tilesize )
            {
                const IndexType bottom = ( iTile + tilesize < n ? tilesize : n - iTile );
                const T* const A = &decomposition[iTile * m + diagonal];

                for( IndexType jTile = diagonal + tilesize; jTile < n; jTile += tilesize )
                {
                    const IndexType right_side = ( jTile + tilesize < m ? tilesize : m - jTile );

                    const T* const B = &decomposition[diagonal * m + jTile];
                    T* C = &decomposition[iTile * m + jTile];

                    OpenMPBLAS3::gemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, bottom, right_side, tile_bottom,
                                       minus, A, n, B, n, one, C, n, NULL );
                }
            }
        }

        decomposition.release();
    }
    /* distributed cyclic implementation */
    else if( matrix.getDistribution() == matrix.getColDistribution() && matrix.getNumRows() == matrix.getNumColumns()
             && typeid( matrix.getColDistribution() ) == typeid(CyclicDistribution) )
    {
        initializeCommunicator();
        const PartitionId myRank = mComm->getRank();
        const PartitionId numProcs = mComm->getSize();

        if( matrix.getNumTotalChunks() < numProcs )
        {
            LAMA_THROWEXCEPTION(
                "No LU decomposition has been computed, because there were too many processes ("<<numProcs<<") from which only "<<matrix.getNumTotalChunks( )<<" would be active." );
        }

        lama_swap<T> swap_cpu;
        swap_cpu.ctxt = ContextFactory::getContext( Context::Host );

        swap_cpu.func = &( OpenMPBLAS1::swap );
//        swap_cpu.func = lama_SWAP_cpu<T>;

        lama_swap<T> swap_cuda;
        lama_gemm<T> gemm;

        ContextPtr ctxtPtr;
        {
            const char* const cArr = getenv( "PLU_CONTEXT" );
            if( cArr == NULL || strcmp( cArr, "HOST" ) == 0 )
            {
                ctxtPtr = ContextFactory::getContext( Context::Host );
                swap_cuda = swap_cpu;
                gemm.stream = 0;
                gemm.func = lama_gemm<T>::__gemm;
                gemm.record = lama_gemm<T>::__recordDef;
                gemm.query = lama_gemm<T>::__queryDef;
                gemm.synchronize = lama_gemm<T>::__synchronizeDef;
            }
            else if( strcmp( cArr, "CUDA" ) == 0 )
            {
#ifdef LAMA_BUILD_CUDA
                //XXX CUDA
                LAMA_LOG_INFO( logger, "Using CUDA" );
                lama_init0_cuda();
                ctxtPtr = ContextFactory::getContext( Context::CUDA, mDev );
                const CUDAContext* cuda = dynamic_cast<const CUDAContext*>( ctxtPtr.get() );

                LAMA_ASSERT_ERROR( cuda, "dynamic cast for CUDAContext failed" );
                lama_gemm<T>::__syncTok = cuda->getComputeSyncToken();

                swap_cuda.ctxt = ctxtPtr;
                swap_cuda.func = lama_SWAP_cuda<T>;
                gemm.stream = lama_gemm<T>::__syncTok->getCUDAStream();
                gemm.func = lama_GEMMAsync_cuda<T>;
                gemm.record = gemm.__recordCuda;
                gemm.query = gemm.__queryCuda;
                gemm.synchronize = gemm.__synchronizeCuda;

                // support fast memory transfer Host->CUDA
                lama::CUDAHostContextManager::setAsCurrent( ctxtPtr );
#else
                LAMA_LOG_WARN( logger,"Specified PLU_CONTEXT as CUDA, but CUDA not found. Assuming host." );
                goto _default;
#endif
            }
            /* Add further Contexts, here. */
            else
            {
#ifndef LAMA_BUILD_CUDA
_default:
#endif
                ctxtPtr = ContextFactory::getContext( Context::Host );

                swap_cuda = swap_cpu;
                gemm.stream = 0;
                gemm.func = lama_gemm<T>::__gemm;
                gemm.record = lama_gemm<T>::__recordDef;
                gemm.query = lama_gemm<T>::__queryDef;
                gemm.synchronize = lama_gemm<T>::__synchronizeDef;
                LAMA_LOG_WARN( logger,
                               "Invalid value for environment variable PLU_CONTEXT: "<<cArr <<". Should be 'HOST' or not set, to be host, or 'CUDA', to be on cuda. Assuming host, now." );
            }
        }

        PartitionId ROOT = 0;
        boost::scoped_array<DenseStorage<T>*> blocks( new DenseStorage<T>*[matrix.getCyclicLocalValues().size()] );
        {
            typedef std::vector<typename DenseMatrix<T>::DenseStoragePtr> vec_type;

            const vec_type& vec = matrix.getCyclicLocalValues();
            for( typename vec_type::size_type i = 0; i < vec.size(); ++i )
            {
                blocks[i] = vec[i].get();
            }
        }

        const IndexType blocks_width = matrix.getNumLocalChunks(); // This might be a bit confusing, but keep in mind,
        const IndexType blocks_height = matrix.getNumTotalChunks(); // that the blocks are stored in column major order.

        const IndexType tilesize = ( *blocks.get() )->mNumColumns;
        const IndexType globalTilesize = mTilesize;
        mTilesize = mComm->max( tilesize );
        boost::scoped_array<IndexType> ipiv( new IndexType[globalTilesize] );

        if( permutation.size() != static_cast<std::vector<IndexType>::size_type>( blocks_width * globalTilesize ) )
        {
            permutation.resize( blocks_width * globalTilesize, 0 );
        }

        const T one = static_cast<T>( 1.0 );
        boost::scoped_array<CUevent> colEvents( new CUevent[blocks_height] );

#ifdef LAMA_BUILD_CUDA
        if( ctxtPtr->getType() == Context::CUDA )
            // initialize the CUDA events
            for( IndexType i = 0; i < blocks_height; ++i )
            {
                LAMA_CUDA_DRV_CALL( cuEventCreate( &colEvents[i], CU_EVENT_DISABLE_TIMING ),
                                    "Could not create event " << i << " for CUDA Queue." );
            }
#endif
        gemm.record( *colEvents.get() );

        LAMAArray<T> bArr( tilesize * tilesize );

        // loop over each block column
        for( IndexType blockCol = 0; blockCol < blocks_height; ++blockCol, ROOT = ( ROOT + 1 ) % numProcs )
        {
            gemm.synchronize( colEvents[blockCol] );
            const IndexType numBlockRows = blocks_width - ( blockCol / numProcs ) - ( myRank < ROOT ? 1 : 0 );

            {   // keep store local.
                DenseStorage<T>** store = NULL;
                if( numBlockRows > 0 )
                {
                    store = &blocks[( blockCol + 1 ) * blocks_width - numBlockRows];
                }
                else
                {
                    store = &blocks[blockCol * blocks_width];
                }

                pgetf2( numBlockRows, store, ipiv.get(), ROOT );
            }

            /* Loop over each block column from the left to the actual and swap. There have to be three distinctions:
             * There are more than one block rows to be updated per block column
             ->  swapping the rows.
             * There is exactly one block row per block column to be updated
             -> The height of the block row may differs from the height of the other block rows
             * There is no block row in any of the block columns
             -> do nothing.
             */
            if( numBlockRows > 1 )
            {
                for( IndexType blockCol_swapL = 0; blockCol_swapL < blockCol; ++blockCol_swapL )
                {
                    plaswp( &blocks[( blockCol_swapL + 1 ) * blocks_width - numBlockRows], ROOT, ipiv.get(), mTilesize,
                            swap_cpu );
                }
            }
            else if( numBlockRows == 1 )
            {
                const IndexType& n = blocks[( blockCol + 1 ) * blocks_width - numBlockRows]->mNumRows;
                for( IndexType blockCol_swapL = 0; blockCol_swapL < blockCol; ++blockCol_swapL )
                {
                    plaswp( &blocks[( blockCol_swapL + 1 ) * blocks_width - numBlockRows], ROOT, ipiv.get(), n,
                            swap_cpu );
                }
            }
            else // numBlockRows == 0
            {
                for( IndexType blockCol_swapL = 0; blockCol_swapL < blockCol; ++blockCol_swapL )
                {
                    plaswp( &blocks[( blockCol_swapL + 1 ) * blocks_width - numBlockRows], ROOT, ipiv.get(), 0,
                            swap_cpu );
                }
            }

            if( myRank == ROOT )
            {
                const IndexType aIndx = blockCol * blocks_width;
                const IndexType firstBlockRow = blockCol / numProcs;

                HostReadAccess<T> diag( blocks[( blockCol + 1 ) * blocks_width - numBlockRows]->getData() );

                // prefetch trsm-blocks to host
                for( IndexType blockColff = blocks_height - 1; blockColff >= blockCol + 1; --blockColff )
                {
                    blocks[( blockColff + 1 ) * blocks_width - numBlockRows]->prefetch( swap_cpu.ctxt );
                }

                for( IndexType blockColff = blockCol + 1; blockColff < blocks_height; ++blockColff )
                {
                    const IndexType blockIdx = ( blockColff + 1 ) * blocks_width - numBlockRows;
                    DenseStorage<T>* firstBlock = blocks[blockIdx];
                    const IndexType n = firstBlock->mNumRows, m = firstBlock->mNumColumns;

                    // 1. LASWP
                    if( numBlockRows > 1 )
                    {
                        plaswp( &blocks[blockIdx], ROOT, ipiv.get(), mTilesize, swap_cuda );
                    }
                    else // numBlockRows == 1
                    {
                        plaswp( &blocks[blockIdx], ROOT, ipiv.get(), n, swap_cuda );
                    }

                    // 2. TRSM
                    HostWriteAccess<T> trsm( firstBlock->getData() );

                    OpenMPBLAS3::trsm( CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, m, one,
                                       diag.get(), tilesize, trsm.get(), m );

                    mComm->bcast( trsm.get(), n * m, myRank );

                    // 3. GEMM
                    trsm.release();

                    if( numBlockRows - 1 > 0 )
                    {
                        ReadAccess<T> B( firstBlock->getData(), ctxtPtr );
                        const IndexType cIndx = blockColff * blocks_width;

                        // loop over each block row of block column
                        for( IndexType blockRow = 1; blockRow < numBlockRows; ++blockRow )
                        {
                            ReadAccess<T> A( blocks[aIndx + blockRow + firstBlockRow]->getData(), ctxtPtr );
                            WriteAccess<T> C( blocks[cIndx + blockRow + firstBlockRow]->getData(), ctxtPtr );

                            const IndexType nLoc =
                                blocks[blockColff * blocks_width + blockRow + firstBlockRow]->mNumRows;

                            gemm.func( CblasRowMajor, CblasNoTrans, CblasNoTrans, nLoc, m, tilesize, -one, A.get(),
                                       tilesize, B.get(), m, one, C.get(), m, gemm.stream );
                        } // for blockRow
                    } // if
                    gemm.record( colEvents[blockColff] );
                } // for blockColff

                // keep permutation global i.e. copy permutation array to the output vector.
                const IndexType permIndx = ( blockCol / numProcs ) * tilesize;
                for( IndexType i = 0; i < tilesize; ++i )
                {
                    permutation[permIndx + i] = ipiv[i];
                }
            }
            else // I'm not root.
            {
                const IndexType aIndx = blockCol * blocks_width;
                const IndexType firstBlockRow = blockCol / numProcs + ( myRank < ROOT ? 1 : 0 );

                // loop over each trailing block column
                for( IndexType blockColff = blockCol + 1; blockColff < blocks_height; ++blockColff )
                {
                    const IndexType blockIdx = numBlockRows > 0 ? ( blockColff + 1 ) * blocks_width - numBlockRows : 0;
                    DenseStorage<T>* firstBlock = blocks[blockIdx];
                    const IndexType m = firstBlock->mNumColumns;

                    // 1. LASWP
                    plaswp( &blocks[blockIdx], ROOT, ipiv.get(), mTilesize, swap_cuda );

                    // 2. Prefetch while root is doing TRSM
                    for( IndexType i = 0; i < numBlockRows; ++i )
                    {
                        blocks[blockColff * blocks_width + firstBlockRow + i]->prefetch( ctxtPtr );
                        blocks[blockCol * blocks_width + firstBlockRow + i]->prefetch( ctxtPtr );
                    }

                    // 3. GEMM
                    HostWriteAccess<T> bAccess( bArr );
                    mComm->bcast( bAccess.get(), bAccess.size(), ROOT ); // TODO why is this (1) OR that (2) expensive?
                    bAccess.release();

                    const IndexType cIndx = blockColff * blocks_width;

                    if( numBlockRows > 0 )
                    {
                        ReadAccess<T> B( bArr, ctxtPtr ); // TODO why is this (2) OR that (1) expensiv?

                        // loop over each block row of block column
                        for( IndexType blockRow = 0; blockRow < numBlockRows; ++blockRow )
                        {
                            ReadAccess<T> A( blocks[aIndx + blockRow + firstBlockRow]->getData(), ctxtPtr );
                            WriteAccess<T> C( blocks[cIndx + blockRow + firstBlockRow]->getData(), ctxtPtr );

                            const IndexType nLoc = blocks[cIndx + blockRow + firstBlockRow]->mNumRows;

                            gemm.func( CblasRowMajor, CblasNoTrans, CblasNoTrans, nLoc, m, tilesize, -one, A.get(),
                                       tilesize, B.get(), m, one, C.get(), m, gemm.stream );
                        }
                    }
                    gemm.record( colEvents[blockColff] );
                }
            }
        } // loop over each block column

        mTilesize = globalTilesize;

#ifdef LAMA_BUILD_CUDA
        if( ctxtPtr->getType() == Context::CUDA )
        {
            // free CUDA events.
            for( IndexType i = 0; i < blocks_height; ++i )
            {
                LAMA_CUDA_DRV_CALL( cuEventDestroy( colEvents[i] ),
                                    "Could not destroy event " << i << " from CUDA Queue." );
            }
            lama_gemm<T>::__syncTok.release();
        }
#endif
    }
    else
    {
        LAMA_THROWEXCEPTION( "LU-factorization not supported for " << matrix << "." );
    }
}

void LUSolver::solve( Vector& solution, const Vector& rhs )
{
    LAMA_REGION("solve");

    LUSolverRuntime& runtime = getRuntime();

    Matrix* luFactorization = runtime.mLUfactorization;

    if( typeid( solution ) == typeid( rhs )
            && ( ( matrixIs<double>( *luFactorization ) && typeid( solution) == typeid(DenseVector<double> ) )
                 || ( matrixIs<float>( *luFactorization )
                      && typeid( solution) == typeid(DenseVector<float> ) ) ) )
    {
        if( luFactorization->getColDistribution().getNumPartitions() != rhs.getDistribution().getNumPartitions()
                || luFactorization->getDistribution().getNumPartitions()
                != solution.getDistribution().getNumPartitions() )
        {
            LAMA_THROWEXCEPTION(
                "Solving the equation system failed due to not fitting partitioning of at least one " << "of the input vectors. The number of column partitions of the matrix was " << luFactorization->getColDistribution( ).getNumPartitions( ) << " whereas the number of partitions " << "of the right hand side vector was "<<rhs.getDistribution( ).getNumPartitions( ) << ". " << "The number of row partitions of the matrix was " << luFactorization->getDistribution( ).getNumPartitions( ) << " whereas the number of partitions of " << "the solution vector was "<<solution.getDistribution( ).getNumPartitions( ) << '.' );
        }

        if( luFactorization->getDistribution().getNumPartitions() == 1
                && luFactorization->getColDistribution().getNumPartitions() == 1 )
        {
            typedef std::vector<IndexType>::size_type size_type;

            if( matrixIs<double>( *luFactorization ) )
            {
                DenseMatrix<double>* denseLUFactorization = dynamic_cast<DenseMatrix<double>*>( luFactorization );
                const DenseVector<double>& denseRhs = dynamic_cast<const DenseVector<double>&>( rhs );
                DenseVector<double>& denseSol = dynamic_cast<DenseVector<double>&>( solution );

                HostReadAccess<double> matrix( denseLUFactorization->getLocalValues().getData() );

                HostWriteAccess<double> vector( denseSol.getLocalValues() );
                vector.resize( denseLUFactorization->getNumColumns() );

                HostReadAccess<double> x( denseRhs.getLocalValues() );

                // Do permutated copy
                for( IndexType i = 0; i < static_cast<IndexType>( runtime.mPermutation.size() ); ++i )
                {
                    vector[i] = x[runtime.mPermutation[i]];
                }

                // solve Ly = b
                OpenMPLAPACK::trtrs( CblasRowMajor, CblasLower, CblasNoTrans, CblasUnit,
                                     denseLUFactorization->getNumColumns(), 1, matrix.get(),
                                     denseLUFactorization->getNumRows(), vector.get(), 1 );

                // solve Rx = y
                OpenMPLAPACK::trtrs( CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                                     denseLUFactorization->getNumColumns(), 1, matrix.get(),
                                     denseLUFactorization->getNumRows(), vector.get(), 1 );
            }
            else if( matrixIs<float>( *luFactorization ) )
            {
                DenseMatrix<float>* denseLUFactorization = dynamic_cast<DenseMatrix<float>*>( luFactorization );
                const DenseVector<float>& denseRhs = dynamic_cast<const DenseVector<float>&>( rhs );
                DenseVector<float>& denseSol = dynamic_cast<DenseVector<float>&>( solution );

                HostReadAccess<float> matrix( denseLUFactorization->getLocalValues().getData() );

                HostWriteAccess<float> vector( denseSol.getLocalValues() );
                vector.resize( denseLUFactorization->getNumColumns() );

                HostReadAccess<float> x( denseRhs.getLocalValues() );

                // Do permutated copy
                for( IndexType i = 0; i < static_cast<IndexType>( runtime.mPermutation.size() ); ++i )
                {
                    vector[i] = x[runtime.mPermutation[i]];
                }

                // solve Ly = b
                OpenMPLAPACK::trtrs( CblasRowMajor, CblasLower, CblasNoTrans, CblasUnit,
                                     denseLUFactorization->getNumColumns(), 1, matrix.get(),
                                     denseLUFactorization->getNumRows(), vector.get(), 1 );

                // solve Rx = y
                OpenMPLAPACK::trtrs( CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                                     denseLUFactorization->getNumColumns(), 1, matrix.get(),
                                     denseLUFactorization->getNumRows(), vector.get(), 1 );
            }
        }
        /* distributed cyclic implementation */
        else if( luFactorization->getDistribution() == luFactorization->getColDistribution()
                 && luFactorization->getNumRows() == luFactorization->getNumColumns()
                 && typeid( luFactorization->getColDistribution() ) == typeid(CyclicDistribution) )
        {
            LAMA_THROWEXCEPTION( "Solving a distributed matrix is producing segmentation faults, yet." );
//            const IndexType myRank   = mComm->getRank( );
//            const IndexType numProcs = mComm->getSize( );
//
//            IndexType FIN = -1;
//            boost::scoped_array<IndexType> indxes( new IndexType[3] );
//
//
//            if( matrixIs<double>( *luFactorization ) )
//            {
//                DenseMatrix<double>* luFactorization = dynamic_cast<DenseMatrix<double>* >( luFactorization );
//                const DenseVector<double>& denseRhs = dynamic_cast<const DenseVector<double>& >( rhs );
//                DenseVector<double>& denseSol = dynamic_cast<DenseVector<double>& >( solution );
//                const IndexType tilesize = mComm->max( luFactorization->getCyclicLocalValues( )[0]->mNumRows );
//
//                HostReadAccess<double>  right( denseRhs.getLocalValues( ) );
//                HostWriteAccess<double> sol  ( denseSol.getLocalValues( ) );
//
//                // copy local values of rhs to solution
//                lama_COPY_cpu( denseRhs.size( ),right.get( ),1,sol.get( ),1 );
//
//                /*
//                 * Performing the swaps. For each swap within the lr-decomposed matrix during the lr-decomposition, a
//                 * swap needs to be done within the vector of the right side (which is static, therefore it will be done
//                 * on the solution vector, holding the values if the right hand side, by now). The swaps are not
//                 * performaned in the order of their appearance, but in the order of the processes. I.e. process by
//                 * process checks its local changes and performes them in communication with the other processes.
//                 */
//                // Loop over each process
//                for( PartitionId ROOT=0;ROOT<numProcs;++ROOT )
//                {
//                    if( myRank == ROOT )
//                    {
//                        const IndexType locSize = denseSol.getLocalValues( ).size( );
//                        const IndexType lastBlock = locSize/tilesize + ( locSize%tilesize != 0 ? 1 : 0 );
//
//                        // loop over each block of the permutation-vector. Skipp the last block. We need to do this for
//                        // each block seperated, because at each block the numbering of the rows starts with zero.
//                        for( IndexType block=0;block<lastBlock-1;++block )
//                        {
//                            for( IndexType i=0;i<tilesize;++i )
//                            {
//                                const IndexType indx = mPermutation[i];
//
//                                if( indx == i )
//                                    continue;
//
//                                IndexType proc = ( indx/tilesize )%numProcs;
//
//                                if( proc == 0 )
//                                {
//                                    std::swap( sol[block*tilesize + i],sol[block*tilesize + indx] );
//                                }
//                                else
//                                {
//                                    mComm->bcast( &proc,1,myRank );
//                                    const IndexType modI = indx%tilesize;
//
//                                    // indxes:
//                                    //  0 - the blockj holding the value
//                                    //  1 - the row within the block
//                                    //  2 - my local block, to start from.
//                                    indxes[0] = ((indx - modI)/tilesize - proc)/numProcs;
//                                    indxes[1] = modI;
//                                    indxes[2] = block;
//                                    mComm->swap( indxes.get( ),3,proc );
//                                    mComm->swap( &sol[block*tilesize + i],1,proc );
//                                }
//                            }
//                        }
//
//                        const IndexType lBsize = locSize%tilesize == 0 ? tilesize : locSize%tilesize;
//
//                        // The last block is seperated, because it has a different size.
//                        for( IndexType i=0;i<lBsize;++i )
//                        {
//                            const IndexType indx = mPermutation[i];
//
//                            if( indx == i )
//                                continue;
//
//                            IndexType proc = ( indx/tilesize )%numProcs;
//
//                            if( proc == 0 )
//                            {
//                                std::swap( sol[( lastBlock-1 )*tilesize + i],sol[( lastBlock-1 )*tilesize + indx] );
//                            }
//                            else
//                            {
//                                mComm->bcast( &proc,1,myRank );
//                                const IndexType modI = indx%tilesize;
//
//                                indxes[0] = ((indx - modI)/tilesize - proc)/numProcs;
//                                indxes[1] = modI;
//                                indxes[2] = lastBlock-1;
//                                mComm->swap( indxes.get( ),3,proc );
//                                mComm->swap( &sol[( lastBlock-1 )*tilesize + i],1,proc );
//                            }
//                        }
//                        mComm->bcast( &FIN,1,myRank );
//                    }
//                    else // I am not ROOT.
//                    {
//                        IndexType TOK = myRank;
//                        for( mComm->bcast( &TOK,1,ROOT );TOK!=FIN;mComm->bcast( &TOK,1,ROOT ) )
//                        {
//                            if( TOK == myRank )
//                            {
//                                mComm->swap( indxes.get( ),3,ROOT );
//                                const IndexType locBlock = indxes[0];
//                                const IndexType locRow   = indxes[1];
//                                const IndexType gloBlock = indxes[2] + ( myRank<ROOT ? 1 : 0 );
//
//                                mComm->swap( &sol[( gloBlock+locBlock )*tilesize + locRow],1,ROOT );
//                            }
//                        }
//                    }
//                }
//
//                right.release( );
//                sol.release( );
//
//                ptrsm( CblasLower,*luFactorization,denseSol );
//                ptrsm( CblasUpper,*luFactorization,denseSol );
//            }
//            else if( matrixIs<float>( *luFactorization ) )
//            {
//                // not implemented, yet, because it is not clear, whether the algorithm for double is bug-free (well,
//                // obviously, its not). If so, the algorithm just needs to be copied and pasted and changed to float
//                // instead of double.
//
//                // DenseMatrix<float>* luFactorization = dynamic_cast<DenseMatrix<float>* >( luFactorization );
//                // const DenseVector<float>& denseRhs = dynamic_cast<const DenseVector<float>& >( rhs );
//                // DenseVector<float>& denseSol = dynamic_cast<DenseVector<float>& >( solution );
//            }
        }
        else
        {
            LAMA_THROWEXCEPTION( "Solving equation system not supported for matrix " << (*luFactorization) << '.' );
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "No solution has been computed, because Matrix was in invalid format." )
    }
}

// XXX TRSM
template<typename T>
void LUSolver::ptrsm( const enum CBLAS_UPLO uplo, const DenseMatrix<T>& matrix, DenseVector<T>& solution )
{
    LAMA_REGION("ptrsm");
    const PartitionId myRank = mComm->getRank();
    const IndexType numProcs = mComm->getSize();

    PartitionId ROOT = 0;

    boost::scoped_array<DenseStorage<T>*> blocks( new DenseStorage<T>*[matrix.getCyclicLocalValues().size()] );
    {
        typedef std::vector<typename DenseMatrix<T>::DenseStoragePtr> vec_type;
        const vec_type& vec = matrix.getCyclicLocalValues();
        for( typename vec_type::size_type i = 0; i < vec.size(); ++i )
        {
            blocks[i] = vec[i].get();
        }
    }

    const IndexType blocks_width = matrix.getNumLocalChunks(); // This might be a bit confusing, but keep in mind,
    const IndexType blocks_height = matrix.getNumTotalChunks(); // that the blocks are stored in column major order.

    const IndexType tilesize = ( *blocks.get() )->mNumColumns;
    const T one = static_cast<T>( 1.0 );
    LAMAArray<T> r( tilesize, static_cast<T>( 0.0 ) );

    enum CBLAS_DIAG diag;
    typedef bool (*fin)( PartitionId&, IndexType&, const IndexType&, const IndexType& );
    fin is;
    typedef void (*initBeginEnd)(
        lama::IndexType&,
        lama::IndexType&,
        const lama::IndexType&,
        const lama::IndexType&,
        const lama::IndexType&,
        const lama::PartitionId&,
        const lama::PartitionId& );
    initBeginEnd init;
    IndexType blockCol;

    // Add the dynamic pointers for solving forwards or backwards.
    if( uplo == CblasLower )
    {
        diag = CblasUnit;
        is = lowerFin; // changes root and increments
        init = lowerInit;
        blockCol = 0;
    }
    else // uplo == CblasUpper
    {
        ROOT = ( numProcs - 1 + blocks_height % numProcs ) % numProcs;
        diag = CblasNonUnit;
        is = upperFin; // changes root and decrements
        init = upperInit;
        blockCol = blocks_height - 1;
    }

    do
    {
        if( myRank == ROOT )
        {
            DenseStorage<T>* trsm = blocks[blockCol * blocks_width + blockCol / numProcs];
            HostReadAccess<T> diagonal( trsm->getData() );
            HostWriteAccess<T> sol( solution.getLocalValues() );

            const IndexType n = trsm->mNumColumns;
            const IndexType idx = ( blockCol / numProcs ) * tilesize;

            // TRSM on diagonal block.
            OpenMPBLAS3::trsm( CblasRowMajor, CblasLeft, uplo, CblasNoTrans, diag, n, 1, one, diagonal.get(), n,
                               &sol[idx], 1 );
            {
                IndexType size = n;
                mComm->bcast( &size, 1, myRank );
            }

            mComm->bcast( &sol[idx], n, myRank );

            HostWriteAccess<T> rw( r );
            OpenMPBLAS1::copy( n, &sol[idx], 1, rw.get(), 1 );
        }
        else
        {
            IndexType size = tilesize;
            HostWriteAccess<T> rw( r );

            mComm->bcast( &size, 1, ROOT );
            mComm->bcast( rw.get(), size, ROOT );
        }

        IndexType blockRowBegin, blockRowEnd;
        init( blockRowBegin, blockRowEnd, blockCol, numProcs, blocks_width, myRank, ROOT );

        HostReadAccess<T> rw( r );

        // GEMV-udate
        for( IndexType blockRow = blockRowBegin; blockRow < blockRowEnd; ++blockRow )
        {
            LAMA_REGION("GEMV_CPU");

            DenseStorage<T>* actBlock = blocks[blockCol * blocks_width + blockRow];
            HostWriteAccess<T> sol( solution.getLocalValues() );
            HostReadAccess<T> mat( actBlock->getData() );

            const IndexType n = actBlock->mNumRows;
            const IndexType m = actBlock->mNumColumns;

            OpenMPBLAS2::gemv( CblasRowMajor, CblasNoTrans, n, m, -one, mat.get(), m, rw.get(), 1, one,
                               &sol[blockRow * tilesize], 1 );
        }
    } while( is( blockCol, ROOT, numProcs, blocks_height ) );
}

// XXX GETF2
template<typename T>
void LUSolver::pgetf2(
    const IndexType numBlockRows,
    DenseStorage<T>** const A,
    IndexType* const ipiv,
    const PartitionId ROOT )
{
    const LAMAInterface* lamaInterface = LAMAInterfaceRegistry::getRegistry().getInterface( Context::Host );

    LAMA_REGION("pgetf2");
    const T MINUS = static_cast<T>( -1.0 );
    const T sfmin = lamaInterface->getLAPACKInterface<T>().lamch( CblasSfmin );

    const PartitionId myRank = mComm->getRank();
    const PartitionId numProc = mComm->getSize();
    const IndexType& numCol = ( *A )->mNumColumns;
    const IndexType& numRows = ( *A )->mNumRows;

    if( myRank == ROOT ) // ROOT has the diagonal element
    {
        for( IndexType j = 0; j < numCol; ++j )
        {
            const int myLocMaxI = piamax_own( numBlockRows, A, j, j );
            const IndexType divI = myLocMaxI / numRows;
            const IndexType modI = myLocMaxI % A[divI]->mNumRows;

            T globalMax, localMax;
            T* m_row = NULL;

            if( divI != 0 )
            {
                HostWriteAccess<T> maxBlock( A[divI]->getData() );
                m_row = &maxBlock[modI * numCol];
            }

            HostWriteAccess<T> diagonal( ( *A )->getData() );

            globalMax = localMax = std::fabs( divI == 0 ? diagonal[modI * numCol + j] : m_row[j] );
            PartitionId maxOwner = myRank;

            mComm->maxloc( globalMax, maxOwner, myRank );

            // Do not swap via communication
            if( maxOwner == myRank || ( localMax > sfmin && std::fabs( globalMax / localMax ) < epsilon ) )
            {
                maxOwner = myRank;
                mComm->bcast( &maxOwner, 1, myRank );
                ipiv[j] = divI * numProc * mTilesize + modI;

                if( j != modI || divI != 0 || ( numCol == 1 && j != divI ) )
                {
                    if( divI == 0 )
                    {
                        OpenMPBLAS1::swap( numCol, &diagonal[j * numCol], 1, &diagonal[modI * numCol], 1 );
                    }
                    else
                    {
                        OpenMPBLAS1::swap( numCol, &diagonal[j * numCol], 1, m_row, 1 );
                    }
                }
            }
            else // do swap via communication
            {
                IndexType indx = myLocMaxI;

                mComm->bcast( &maxOwner, 1, myRank );
                mComm->swap( &indx, 1, maxOwner );
                mComm->swap( &diagonal[j * numCol], numCol, maxOwner );

                ipiv[j] = indx;
            }

            const T diag = diagonal[j * numCol + j];
            mComm->bcast( &diagonal[j * numCol + j], numCol - j, myRank );

            if( diag == 0 )
            {
                LAMA_THROWEXCEPTION( "Dividing through zero" );
            }

            if( j < numCol - 1 )
            {
                const T* const y = &diagonal[j * numCol + j + 1];

                T* const x = &diagonal[( j + 1 ) * numCol + j];
                T* const a = &diagonal[( j + 1 ) * numCol + j + 1];

                if( diag >= sfmin )
                {
                    OpenMPBLAS1::scal( numCol - 1 - j, 1 / diag, x, numCol );
                }
                else
                {
                    for( IndexType i = 0; i < numCol - 1 - j; ++i )
                    {
                        x[i * numCol] /= diag;
                    }
                }

                // like gemm but for vectors
                OpenMPBLAS2::ger( CblasRowMajor, numCol - 1 - j, numCol - 1 - j, MINUS, x, numCol, y, 1, a, numCol );

                for( IndexType block = 1; block < numBlockRows; ++block )
                {
                    HostWriteAccess<T> sub( A[block]->getData() );
                    const IndexType& numRow = A[block]->mNumRows;

                    T* const x = &sub[j];
                    T* const a = &sub[j + 1];

                    if( diag >= sfmin )
                    {
                        OpenMPBLAS1::scal( numRow, 1 / diag, x, numCol );
                    }
                    else
                    {
                        for( IndexType i = 0; i < numRow; ++i )
                        {
                            x[i * numCol] /= diag;
                        }
                    }
                    LAMA_REGION("GER_CPU");
                    OpenMPBLAS2::ger( CblasRowMajor, numRow, numCol - 1 - j, MINUS, x, numCol, y, 1, a, numCol );
                } // for
            }
            else
            {
                for( IndexType block = 1; block < numBlockRows; ++block )
                {
                    HostWriteAccess<T> sub( A[block]->getData() );
                    T* const x = &sub[j];
                    const IndexType& numRow = A[block]->mNumRows;

                    if( diag >= sfmin )
                    {
                        OpenMPBLAS1::scal( numRow, 1 / diag, x, numCol );
                    }
                    else
                    {
                        for( IndexType i = 0; i < numRow; ++i )
                        {
                            x[i * numCol] /= diag;
                        }
                    }
                } // for
            }
        }
    }
    else // myRank != ROOT
    {
        boost::scoped_array<T> y( new T[numCol] );

        if( numBlockRows > 0 )
        {
            for( IndexType j = 0; j < numCol; ++j )
            {
                const int myLocMaxI = piamax_own( numBlockRows, A, j );
                const IndexType divI = myLocMaxI / numCol;
                const IndexType modI = myLocMaxI % numCol;

                HostWriteAccess<T> maxBlock( A[divI]->getData() );
                {
                    PartitionId rank = myRank;
                    T myVal = std::fabs( maxBlock[modI * numCol + j] );
                    mComm->maxloc( myVal, rank, ROOT );
                }

                PartitionId maxOwner = myRank;
                mComm->bcast( &maxOwner, 1, ROOT );

                if( maxOwner == myRank )
                {
                    // divI*numProc                     : go to global block as if myRank was ROOT.
                    // (numProc+myRank-ROOT)%numProc    : add the missing blocks to the actual ROOT.
                    IndexType lIndx = ( divI * numProc + ( numProc + myRank - ROOT ) % numProc ) * mTilesize + modI;
                    mComm->swap( &lIndx, 1, ROOT );
                    mComm->swap( &maxBlock[modI * numCol], mTilesize, ROOT );
                }

                maxBlock.release();

                mComm->bcast( y.get() + j, numCol - j, ROOT );
                const T diag = y[j];

                if( diag == 0 )
                {
                    LAMA_THROWEXCEPTION( "Dividing through zero" );
                }

                if( j < numCol - 1 )
                {
                    for( IndexType block = 0; block < numBlockRows; ++block )
                    {
                        HostWriteAccess<T> sub( A[block]->getData() );
                        const IndexType& numRow = A[block]->mNumRows;

                        T* const x = &sub[j];
                        T* const a = &sub[j + 1];

                        if( diag >= sfmin )
                        {
                            OpenMPBLAS1::scal( numRow, 1 / diag, x, numCol );
                        }
                        else
                        {
                            for( IndexType i = 0; i < numRow; ++i )
                            {
                                x[i * numCol] /= diag;
                            }
                        }
                        LAMA_REGION("GER_CPU");
                        OpenMPBLAS2::ger( CblasRowMajor, numRow, numCol - 1 - j, MINUS, x, numCol, y.get() + j + 1, 1,
                                          a, numCol );
                    }
                }
                else
                {
                    for( IndexType block = 0; block < numBlockRows; ++block )
                    {
                        LAMA_REGION("SCAL_CPU");
                        HostWriteAccess<T> sub( A[block]->getData() );
                        const IndexType& numRow = A[block]->mNumRows;
                        T* const x = &sub[j];

                        if( diag >= sfmin )
                        {
                            OpenMPBLAS1::scal( numRow, 1 / diag, x, numCol );
                        }
                        else
                        {
                            for( IndexType i = 0; i < numRow; ++i )
                            {
                                x[i * numCol] /= diag;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            PartitionId rank = myRank;
            T zero = static_cast<T>( 0.0 );

            for( IndexType j = 0; j < numCol; ++j )
            {
                mComm->maxloc( zero, rank, ROOT );
                mComm->bcast( &rank, 1, ROOT );
                mComm->bcast( y.get(), numCol - j, ROOT );
            }
        }
    }
}

// XXX IAMAX
template<typename T>
int LUSolver::piamax_own(
    const IndexType numBlockCol,
    DenseStorage<T>** const local,
    const IndexType col,
    const IndexType locRow/*=0*/)
{
    LAMA_REGION("piamax_own");

    ContextPtr context = getCoefficients().getContextPtr();

    boost::scoped_array<IndexType> indices( new IndexType[numBlockCol] );
    boost::scoped_array<T> values( new T[numBlockCol] );

    const IndexType beginLoop = locRow > 0 ? 1 : 0;

    if( beginLoop == 1 )
    {
        DenseStorage<T>* storage = *local;
        HostReadAccess<T> read( storage->getData() );

        const IndexType lda = storage->mNumColumns;

        IndexType indx = OpenMPBLAS1::iamax( storage->mNumRows - locRow, &read[locRow * lda + col], lda ) + locRow;

        indices[0] = indx;
        values[0] = read[indx * lda + col];
    }

    // Receive from all blocks the indices of the maximum. Because the maximum of maxima is needed, the maxima and
    // their indices need to be stored in seperated arrays.
    for( IndexType block = beginLoop; block < numBlockCol; ++block )
    {
        DenseStorage<T>* storage = local[block];
        HostReadAccess<T> read( storage->getData() );
        IndexType indx = OpenMPBLAS1::iamax( storage->mNumRows, &read[col], storage->mNumColumns );

        indices[block] = indx + block * storage->mNumRows;
        values[block] = read[indx * storage->mNumColumns + col];
    }

    IndexType indx = OpenMPBLAS1::iamax( numBlockCol, values.get(), 1 );

    return indices[indx];
}

// XXX LASWP
template<typename T>
void LUSolver::plaswp(
    DenseStorage<T> * * const A,
    const PartitionId ROOT,
    const IndexType * const ipiv,
    const IndexType n,
    const lama_swap<T> swap )
{
    LAMA_REGION("plaswp");
    const PartitionId myRank = mComm->getRank();
    const PartitionId numProc = mComm->getSize();
    IndexType FIN = -1;
    boost::scoped_array<IndexType> indxes( new IndexType[2] );

    const IndexType& tilesize = ( *A )->mNumColumns;

    if( myRank == ROOT )
    {
        for( IndexType i = 0; i < n; ++i )
        {
            const IndexType indx = ipiv[i];

            if( indx == i )
            {
                continue;
            }

            IndexType proc = ( indx / mTilesize ) % numProc;

            if( proc == 0 )
            {
                WriteAccess<T> diagonal( ( *A )->getData(), swap.ctxt );

                if( ( indx / numProc ) / mTilesize == 0 )
                {
                    LAMA_REGION("LASWAP");
                    swap.func( tilesize, diagonal.get() + i * tilesize, 1,
                               diagonal.get() + ( indx % mTilesize ) * tilesize, 1 );
                }
                else // indx/tilesize != 0
                {
                    LAMA_REGION("LASWAP");
                    WriteAccess<T> other( A[( indx / numProc ) / mTilesize]->getData(), swap.ctxt );
                    swap.func( tilesize, diagonal.get() + i * tilesize, 1,
                               other.get() + ( indx % mTilesize ) * tilesize, 1 );
                }
            }
            else // swappable line is not local
            {
                PartitionId owner = ( proc + myRank ) % numProc;
                mComm->bcast( &owner, 1, myRank );
                const IndexType modI = indx % mTilesize;

                indxes[0] = ( ( indx - modI ) / mTilesize - proc ) / numProc;
                indxes[1] = modI;
                mComm->swap( indxes.get(), 2, owner );

                HostWriteAccess<T> diagonal( ( *A )->getData() );
                mComm->swap( &diagonal[i * tilesize], tilesize, owner );
            }
        }
        mComm->bcast( &FIN, 1, myRank );
    }
    else // I'm not ROOT
    {
        IndexType TOK = myRank;

        for( mComm->bcast( &TOK, 1, ROOT ); TOK != FIN; mComm->bcast( &TOK, 1, ROOT ) )
        {
            if( TOK == myRank )
            {
                mComm->swap( indxes.get(), 2, ROOT );

                HostWriteAccess<T> other( A[indxes[0]]->getData() );
                mComm->swap( &other[indxes[1] * tilesize], tilesize, ROOT );
            }
        }
    }
}

IndexType LUSolver::getTileSize()
{
    return mTilesize;
}

void LUSolver::setTileSize( const IndexType tilesize )
{
    if( tilesize <= 0 )
    {
        LAMA_LOG_WARN( logger, "Setting tilesize of " << tilesize << " left LUSolver unaffected." );
    }
    else
    {
        mTilesize = tilesize;
    }
}

IndexType LUSolver::computeTilesize( IndexType m, IndexType n )
{
    return ( mTilesize < std::min( m, n ) ? mTilesize : std::min( m, n ) );
}

void LUSolver::setDeviceNumber( const IndexType dev )
{
    if( dev < -1 )
    {
        LAMA_LOG_WARN( logger, "Setting tilesize of " << dev << " left LUSolver unaffected." );
    }
    else
    {
        mDev = dev;
    }
}

IndexType LUSolver::getDeviceNumber()
{
    return mDev;
}

LUSolver::LUSolverRuntime& LUSolver::getRuntime()
{
    return mLUSolverRuntime;
}

const LUSolver::LUSolverRuntime& LUSolver::getConstRuntime() const
{
    return mLUSolverRuntime;
}

SolverPtr LUSolver::copy()
{
    return SolverPtr( new LUSolver( *this ) );
}

} // namespace LAMA
