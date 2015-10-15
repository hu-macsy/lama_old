/**
 * @file MatrixCreator.cpp
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
 * @brief MatrixCreator.cpp
 * @author Thomas Brandes
 * @date 01.02.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/matutils/MatrixCreator.hpp>

// local library
#include <scai/lama/distribution/BlockDistribution.hpp>
#include <scai/lama/distribution/GeneralDistribution.hpp>

// internal scai libraries
#include <scai/hmemo/WriteAccess.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <cmath>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, MatrixCreator<ValueType>::logger, "MatrixCreator" )

static inline void getStencilProperties(
    IndexType& dimension,
    IndexType& length,
    IndexType& distance,
    IndexType stencilType )
{
    switch( stencilType )
    {
        case 27:
            length = 1;
            distance = 3;
            dimension = 3;
            break;

        case 19:
            length = 1;
            distance = 2;
            dimension = 3;
            break;

        case 7:
            length = 1;
            distance = 1;
            dimension = 3;
            break;

        case 9:
            length = 1;
            distance = 2;
            dimension = 2;
            break;

        case 5:
            length = 1;
            distance = 1;
            dimension = 2;
            break;

        case 3:
            length = 1;
            distance = 1;
            dimension = 1;
            break;

        default:
            length = 1;
            distance = 1;
            dimension = 1;
            COMMON_THROWEXCEPTION( "Invalid Stencil type = " << stencilType )
    }
}

template<typename ValueType>
bool MatrixCreator<ValueType>::supportedStencilType( const IndexType dimension, const IndexType stencilType )
{
    bool supported = false;

    try
    {
        IndexType dim;
        IndexType length;
        IndexType distance;
        getStencilProperties( dim, length, distance, stencilType );
        supported = ( dim == dimension );
    }
    catch( ... )
    {
        supported = false; // so it is not supported
    }

    return supported;
}

template<typename IndexType>
static inline IndexType getMatrixPosition(
    IndexType idX,
    IndexType idY,
    IndexType idZ,
    IndexType dimX,
    IndexType dimY,
    IndexType /* dimZ */)
{
    return idZ * dimX * dimY + idY * dimX + idX;
}

template<typename IndexType>
static inline int getNumNeighbors( IndexType id, IndexType dim, int shift )
{
    int neighbors = ::abs( shift );

    if( shift < 0 )
    {
        int leftBorder = id;

        if( leftBorder < neighbors )
        {
            neighbors = leftBorder;
        }
    }
    else
    {
        int rightBorder = dim - 1 - id;

        if( rightBorder < neighbors )
        {
            neighbors = rightBorder;
        }
    }

    // printf("%d of %d, shift = %d, has %d neighbors\n", id, dim, shift, neighbors );

    return neighbors;
}

template<typename IndexType>
static inline IndexType getNStencilValues(
    IndexType idX,
    IndexType idY,
    IndexType idZ,
    IndexType dimX,
    IndexType dimY,
    IndexType dimZ,
    IndexType length,
    IndexType maxDistance )
{
    SCAI_ASSERT_DEBUG( idX < dimX, "idX = " << idX << " out of range, dimX = " << dimX )
    SCAI_ASSERT_DEBUG( idY < dimY, "idY = " << idY << " out of range, dimY = " << dimY )
    SCAI_ASSERT_DEBUG( idZ < dimZ, "idZ = " << idZ << " out of range, dimZ = " << dimZ )

    IndexType nX = getNumNeighbors( idX, dimX, -length ) + getNumNeighbors( idX, dimX, length );
    IndexType nY = getNumNeighbors( idY, dimY, -length ) + getNumNeighbors( idY, dimY, length );
    IndexType nZ = getNumNeighbors( idZ, dimZ, -length ) + getNumNeighbors( idZ, dimZ, length );

    // SCAI_LOG_DEBUG( logger, idX << "," << idY << "," << idZ << ": neighbors = "
    //                   << nX << "," << nY << "," << nZ )

    // printf("%d,%d,%d has %d,%d,%d neighbors\n", idX, idY, idZ, nX, nY, nZ);

    // Note: nZ will be 0 for dim <= 2, nZ will be 0 for dim <= 1

    IndexType numValues = 1 + nX + nY + nZ; // maxDistance 1

    if( maxDistance >= 2 )
    {
        SCAI_ASSERT_ERROR( length == 1, "length > 1 for stencil not supported" )
        numValues += nX * nY + nX * nZ + nY * nZ;
    }

    if( maxDistance >= 3 )
    {
        SCAI_ASSERT_ERROR( length == 1, "length > 1 for stencil not supported" )
        // dimension cannot be 1, 2
        numValues += nX * nY * nZ;
    }

    return numValues;
}

template<typename IndexType,typename ValueType>
static inline void getStencil(
    std::vector<IndexType>& positions,
    std::vector<ValueType>& values,
    const IndexType idX,
    const IndexType idY,
    const IndexType idZ,
    const IndexType dimX,
    const IndexType dimY,
    const IndexType dimZ,
    const IndexType stencilType,
    const IndexType length,
    const IndexType maxDistance )
{
    positions.clear(); // reset number of entries to 0
    values.clear(); // reset

    positions.push_back( getMatrixPosition( idX, idY, idZ, dimX, dimY, dimZ ) );
    values.push_back( stencilType - static_cast<ValueType>(1.0) );

    IndexType leftX = getNumNeighbors( idX, dimX, -length );
    IndexType rightX = getNumNeighbors( idX, dimX, length );
    IndexType leftY = getNumNeighbors( idY, dimY, -length );
    IndexType rightY = getNumNeighbors( idY, dimY, length );
    IndexType leftZ = getNumNeighbors( idZ, dimZ, -length );
    IndexType rightZ = getNumNeighbors( idZ, dimZ, length );

    for( IndexType jz = idZ - leftZ; jz <= idZ + rightZ; ++jz )
    {
        IndexType distZ = std::abs( jz - idZ );

        for( IndexType jy = idY - leftY; jy <= idY + rightY; ++jy )
        {
            IndexType distYZ = distZ + std::abs( jy - idY );

            if( distYZ > maxDistance )
            {
                continue; // can already skip next loop
            }

            for( IndexType jx = idX - leftX; jx <= idX + rightX; ++jx )
            {
                IndexType distXYZ = distYZ + std::abs( jx - idX );

                if( distXYZ > maxDistance )
                {
                    continue;
                }

                // skip the diagonal element already added

                if( distXYZ == 0 )
                {
                    continue;
                }

                positions.push_back( getMatrixPosition( jx, jy, jz, dimX, dimY, dimZ ) );
                values.push_back( static_cast<ValueType>(-1.0) );
            }
        }
    }
}

template<typename ValueType>
void MatrixCreator<ValueType>::buildPoisson(
    CSRSparseMatrix<ValueType>& matrix,
    const IndexType dimension,
    const IndexType stencilType,
    const IndexType dimX,
    const IndexType dimY,
    const IndexType dimZ )
{
    // Calculate subdomains, subranges

    PartitionId gridSize[3] =
    { 1, 1, 1 };
    PartitionId gridRank[3] =
    { 0, 0, 0 };
    PartitionId dimLB[3] =
    { 0, 0, 0 };
    PartitionId dimUB[3] =
    { dimX - 1, dimY - 1, dimZ - 1 };

    // ToDo: take communicator from input set

    scai::lama::CommunicatorPtr comm = scai::lama::Communicator::get( "MPI" );

    // get rank of this processor

    if( dimension == 1 )
    {
        gridSize[0] = comm->getSize();
        gridRank[0] = comm->getRank();
        scai::lama::BlockDistribution::getRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
    }
    else if( dimension == 2 )
    {
        comm->factorize2( dimX, dimY, gridSize );
        comm->getGrid2Rank( gridRank, gridSize );
        scai::lama::BlockDistribution::getRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
        scai::lama::BlockDistribution::getRange( dimLB[1], dimUB[1], dimY, gridRank[1], gridSize[1] );
    }
    else if( dimension == 3 )
    {
        comm->factorize3( dimX, dimY, dimZ, gridSize );
        comm->getGrid3Rank( gridRank, gridSize );
        scai::lama::BlockDistribution::getRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
        scai::lama::BlockDistribution::getRange( dimLB[1], dimUB[1], dimY, gridRank[1], gridSize[1] );
        scai::lama::BlockDistribution::getRange( dimLB[2], dimUB[2], dimZ, gridRank[2], gridSize[2] );
    }

    SCAI_LOG_INFO( logger,
                   *comm << ": rank = (" << gridRank[0] << "," << gridRank[1] << "," << gridRank[2] << ") of (" << gridSize[0] << "," << gridSize[1] << "," << gridSize[2] << "), local range = [" << dimLB[0] << ":" << dimUB[0] << "," << dimLB[1] << ":" << dimUB[1] << "," << dimLB[2] << ":" << dimUB[2] << "] of " << dimX << " x " << dimY << " x " << dimZ )

    IndexType globalSize = dimX * dimY * dimZ; // number of rows, columns of full matrix

    std::vector<IndexType> myGlobalIndexes; // row indexes of this processor
    std::vector<IndexType> myIA; // number of entries in my rows
    IndexType myNNA = 0; // accumulated sum for number of my entries

    // compute local size on this processor

    IndexType localSize = 1;

    for( int i = 0; i < dimension; i++ )
    {
        // avoid negative sizes, can happen #procs >> #ndim
        if( dimLB[i] > dimUB[i] )
        {
            localSize = 0;
            break;
        }

        localSize *= ( dimUB[i] - dimLB[i] + 1 );
    }

    SCAI_LOG_DEBUG( logger, *comm << ": has local size = " << localSize )

    myGlobalIndexes.reserve( localSize );
    myIA.reserve( localSize );

    IndexType length;
    IndexType dimStencil;
    IndexType maxDistance;

    getStencilProperties( dimStencil, length, maxDistance, stencilType );

    SCAI_ASSERT_EQUAL_ERROR( dimStencil, dimension )

    SCAI_LOG_INFO( logger,
                   "stencil type = " << stencilType << " -> dim = " << dimStencil << ", direction length = " << length << ", max distance = " << maxDistance )

    // compute global indexes this processor is responsibile for and number of non-zero values
    // of the rows owned by this processor

    for( IndexType idZ = dimLB[2]; idZ <= dimUB[2]; ++idZ )
    {
        for( IndexType idY = dimLB[1]; idY <= dimUB[1]; ++idY )
        {
            for( IndexType idX = dimLB[0]; idX <= dimUB[0]; ++idX )
            {
                const IndexType iPos = getMatrixPosition( idX, idY, idZ, dimX, dimY, dimZ );

                myGlobalIndexes.push_back( iPos );

                const IndexType numNonZeros = getNStencilValues( idX, idY, idZ, dimX, dimY, dimZ, length, maxDistance );

                myIA.push_back( numNonZeros );

                // accumulate number of stencil values

                myNNA += numNonZeros;
            }
        }
    }

    SCAI_LOG_INFO( logger, *comm << ": has local " << localSize << " rows, nna = " << myNNA )
    // allocate and fill local part of the distributed matrix

    scai::lama::DistributionPtr distribution = scai::lama::DistributionPtr(
            new scai::lama::GeneralDistribution( globalSize, myGlobalIndexes, comm ) );

    SCAI_LOG_INFO( logger, "distribution = " << *distribution )

    // create new local CSR data ( # local rows x # columns )

    scai::lama::CSRStorage<ValueType> localMatrix;
    localMatrix.allocate( localSize, globalSize );

    // Allocate local matrix with correct sizes and correct first touch in case of OpenMP
    // ToDo: localMatrix( localSize, numColumns, numNonZeros, &myIA[0] );

    hmemo::LAMAArray<IndexType> csrIA;
    hmemo::LAMAArray<IndexType> csrJA;
    hmemo::LAMAArray<ValueType> csrValues;

    {
        hmemo::WriteOnlyAccess<IndexType> ia( csrIA, localSize + 1 );
        hmemo::WriteOnlyAccess<IndexType> ja( csrJA, myNNA );
        hmemo::WriteOnlyAccess<ValueType> values( csrValues, myNNA );

        ia[0] = 0;

        std::vector<IndexType> colIndexes;
        std::vector<double> colValues;
        ;
        colIndexes.reserve( stencilType );
        colValues.reserve( stencilType );

        // compute global indexes this processor is responsibile for and number of non-zero values
        // Important: same loop order as above

        IndexType rowCounter = 0; // count local rows
        IndexType nnzCounter = 0; // count local non-zero elements

        for( IndexType idZ = dimLB[2]; idZ <= dimUB[2]; ++idZ )
        {
            for( IndexType idY = dimLB[1]; idY <= dimUB[1]; ++idY )
            {
                for( IndexType idX = dimLB[0]; idX <= dimUB[0]; ++idX )
                {
                    // get column positions and values of matrix, diagonal element is first

                    getStencil( colIndexes, colValues, idX, idY, idZ, dimX, dimY, dimZ, stencilType, length,
                                maxDistance );

                    // check colIndexes.size() against number of values given by getNStencilValues

                    SCAI_ASSERT_EQUAL_DEBUG( (int ) colIndexes.size(),
                                             getNStencilValues( idX, idY, idZ, dimX, dimY, dimZ, length, maxDistance ) );

                    SCAI_ASSERT_EQUAL_DEBUG( colIndexes.size(), colValues.size() )

                    SCAI_LOG_TRACE( logger,
                                    *comm << ": at " << idX << " x " << idY << " x " << idZ << ", local row : " << rowCounter << ", global row : " << getMatrixPosition( idX, idY, idZ, dimX, dimY, dimZ ) << ": " << colIndexes.size() << " entries: " << colIndexes[0] << ": " << colValues[0] << ", ..." )

                    ia[rowCounter + 1] = ia[rowCounter] + static_cast<IndexType>( colIndexes.size() );

                    for( size_t k = 0; k < colIndexes.size(); ++k )
                    {
                        ja[nnzCounter] = colIndexes[k];
                        values[nnzCounter] = static_cast<ValueType>( colValues[k] );
                        ++nnzCounter;
                    }

                    ++rowCounter;
                }
            }
        }
    }

    localMatrix.swap( csrIA, csrJA, csrValues );

    SCAI_LOG_DEBUG( logger, "replace owned data with " << localMatrix )

    matrix.assign( localMatrix, distribution, distribution ); // builds also halo

    // but now the local part of matrixA should have the diagonal property as global column
    // indexes have been localized

    SCAI_ASSERT_DEBUG( matrix.getLocalStorage().hasDiagonalProperty(), "CSR data has not diagonal property" )

    SCAI_LOG_INFO( logger, "built matrix A = " << matrix )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixCreator<ValueType>::buildPoisson1D(
    CSRSparseMatrix<ValueType>& matrix,
    const IndexType stencilType,
    const IndexType dim )
{
    SCAI_LOG_INFO( logger, "build Poisson1D" << stencilType << "P( " << dim << ")" )

    buildPoisson( matrix, 1, stencilType, dim, 1, 1 );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixCreator<ValueType>::buildPoisson2D(
    CSRSparseMatrix<ValueType>& matrix,
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2 )
{
    SCAI_LOG_INFO( logger, "build Poisson2D" << stencilType << "P( " << dim1 << ", " << dim2 << ")" )

    buildPoisson( matrix, 2, stencilType, dim1, dim2, 1 );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixCreator<ValueType>::buildPoisson3D(
    CSRSparseMatrix<ValueType>& matrix,
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2,
    const IndexType dim3 )
{
    SCAI_LOG_INFO( logger, "build Poisson3D" << stencilType << "P( " << dim1 << ", " << dim2 << ", " << dim3 << ")" )

    buildPoisson( matrix, 3, stencilType, dim1, dim2, dim3 );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixCreator<ValueType>::fillRandom( Matrix& matrix, double density )
{
    int seed = 15191;

    // Shape and distribution of matrix is not changed

    const Distribution& dist = matrix.getDistribution();

    if( dist.getNumPartitions() == 1 )
    {
        // all processors must have the same random numbers

        std::srand( seed );
    }
    else
    {
        // processors must generate different numbers

        PartitionId rank = dist.getCommunicator().getRank();
        std::srand( seed + rank );
    }

    const IndexType localRowSize = dist.getLocalSize();
    const IndexType colSize = matrix.getNumColumns();

    const IndexType expectedEntries = static_cast<IndexType>( localRowSize * colSize * density + 30.0 );

    std::vector<IndexType> csrIA( localRowSize + 1 );
    std::vector<IndexType> csrJA;
    std::vector<ValueType> csrValues;

    csrJA.reserve( expectedEntries );
    csrValues.reserve( expectedEntries );

    IndexType numValues = 0;

    csrIA[0] = numValues;

    for( int i = 0; i < localRowSize; ++i )
    {
        for( int j = 0; j < colSize; ++j )
        {
            ValueType value = static_cast<ValueType>( rand() ) / static_cast<ValueType>( RAND_MAX );

            if( value < density )
            {
                value = static_cast<ValueType>( rand() ) / static_cast<ValueType>( RAND_MAX );
                csrJA.push_back( j );
                csrValues.push_back( value );
                ++numValues;
            }
        }

        csrIA[i + 1] = numValues;
    }

    CSRStorage<ValueType> localCSRStorage;

    localCSRStorage.setRawCSRData( localRowSize, colSize, numValues, &csrIA[0], &csrJA[0], &csrValues[0] );

    SCAI_LOG_DEBUG( logger, "replace owned data with " << localCSRStorage )

    // The new matrix data has the same row distribution as the input
    // matrix, also take over the original column distribution

    matrix.assign( localCSRStorage, matrix.getDistributionPtr(), matrix.getColDistributionPtr() ); // builds also halo
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixCreator<ValueType>::buildRandom(
    CSRSparseMatrix<ValueType>& matrix,
    const IndexType size,
    const double density )
{
    CommunicatorPtr comm = scai::lama::Communicator::get( "MPI" );

    DistributionPtr dist( new BlockDistribution( size, comm ) );
    matrix.allocate( dist, dist );
    fillRandom( matrix, density );
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

#define LAMA_MATRIX_CREATOR_INSTANTIATE(z, I, _)                                     \
    template class COMMON_DLL_IMPORTEXPORT MatrixCreator<ARITHMETIC_HOST_TYPE_##I> ;

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_MATRIX_CREATOR_INSTANTIATE, _ )

#undef LAMA_MATRIX_CREATOR_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
