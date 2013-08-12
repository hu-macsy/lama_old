/**
 * @file MatrixCreator.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
#include <lama/matutils/MatrixCreator.hpp>

// others
#include <lama/HostWriteAccess.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/GeneralDistribution.hpp>


#include <sys/time.h>
#include <cmath>

namespace lama

{

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, MatrixCreator<ValueType>::logger, "MatrixCreator" )

static inline void getStencilProperties( IndexType& dimension, IndexType& length, IndexType& distance, IndexType stencilType )
{
    switch ( stencilType )
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
        LAMA_THROWEXCEPTION( "Invalid Stencil type = " << stencilType )
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
    catch ( ... )
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

    if ( shift < 0 )
    {
        int leftBorder = id;

        if ( leftBorder < neighbors )
        {
            neighbors = leftBorder;
        }
    }
    else
    {
        int rightBorder = dim - 1 - id;
        if ( rightBorder < neighbors )
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
    LAMA_ASSERT_DEBUG( idX < dimX, "idX = " << idX << " out of range, dimX = " << dimX )
    LAMA_ASSERT_DEBUG( idY < dimY, "idY = " << idY << " out of range, dimY = " << dimY )
    LAMA_ASSERT_DEBUG( idZ < dimZ, "idZ = " << idZ << " out of range, dimZ = " << dimZ )

    IndexType nX = getNumNeighbors( idX, dimX, -length ) + getNumNeighbors( idX, dimX, length );
    IndexType nY = getNumNeighbors( idY, dimY, -length ) + getNumNeighbors( idY, dimY, length );
    IndexType nZ = getNumNeighbors( idZ, dimZ, -length ) + getNumNeighbors( idZ, dimZ, length );

    // LAMA_LOG_DEBUG( logger, idX << "," << idY << "," << idZ << ": neighbors = "
    //                   << nX << "," << nY << "," << nZ )

    // printf("%d,%d,%d has %d,%d,%d neighbors\n", idX, idY, idZ, nX, nY, nZ);

    // Note: nZ will be 0 for dim <= 2, nZ will be 0 for dim <= 1

    IndexType numValues = 1 + nX + nY + nZ; // maxDistance 1

    if ( maxDistance >= 2 )
    {
        LAMA_ASSERT_ERROR( length == 1, "length > 1 for stencil not supported" )
        numValues += nX * nY + nX * nZ + nY * nZ;
    }
    if ( maxDistance >= 3 )
    {
        LAMA_ASSERT_ERROR( length == 1, "length > 1 for stencil not supported" )
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
    values.push_back( stencilType - 1.0 );

    IndexType leftX = getNumNeighbors( idX, dimX, -length );
    IndexType rightX = getNumNeighbors( idX, dimX, length );
    IndexType leftY = getNumNeighbors( idY, dimY, -length );
    IndexType rightY = getNumNeighbors( idY, dimY, length );
    IndexType leftZ = getNumNeighbors( idZ, dimZ, -length );
    IndexType rightZ = getNumNeighbors( idZ, dimZ, length );

    for ( IndexType jz = idZ - leftZ; jz <= idZ + rightZ; ++jz )
    {
        IndexType distZ = std::abs( jz - idZ );

        for ( IndexType jy = idY - leftY; jy <= idY + rightY; ++jy )
        {
            IndexType distYZ = distZ + std::abs( jy - idY );

            if ( distYZ > maxDistance )
            {
                continue; // can already skip next loop
            }

            for ( IndexType jx = idX - leftX; jx <= idX + rightX; ++jx )
            {
                IndexType distXYZ = distYZ + std::abs( jx - idX );

                if ( distXYZ > maxDistance )
                {
                    continue;
                }

                // skip the diagonal element already added

                if ( distXYZ == 0 )
                {
                    continue;
                }

                positions.push_back( getMatrixPosition( jx, jy, jz, dimX, dimY, dimZ ) );
                values.push_back( -1.0 );
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

    lama::PartitionId gridSize[3] =
    { 1, 1, 1 };
    lama::PartitionId gridRank[3] =
    { 0, 0, 0 };
    lama::PartitionId dimLB[3] =
    { 0, 0, 0 };
    lama::PartitionId dimUB[3] =
    { dimX - 1, dimY - 1, dimZ - 1 };

    // ToDo: take communicator from input set

    lama::CommunicatorPtr comm = lama::CommunicatorFactory::get( "MPI" );

    // get rank of this processor

    if ( dimension == 1 )
    {
        gridSize[0] = comm->getSize();
        gridRank[0] = comm->getRank();
        lama::BlockDistribution::getRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
    }
    else if ( dimension == 2 )
    {
        comm->factorize2( dimX, dimY, gridSize );
        comm->getGrid2Rank( gridRank, gridSize );
        lama::BlockDistribution::getRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
        lama::BlockDistribution::getRange( dimLB[1], dimUB[1], dimY, gridRank[1], gridSize[1] );
    }
    else if ( dimension == 3 )
    {
        comm->factorize3( dimX, dimY, dimZ, gridSize );
        comm->getGrid3Rank( gridRank, gridSize );
        lama::BlockDistribution::getRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
        lama::BlockDistribution::getRange( dimLB[1], dimUB[1], dimY, gridRank[1], gridSize[1] );
        lama::BlockDistribution::getRange( dimLB[2], dimUB[2], dimZ, gridRank[2], gridSize[2] );
    }

    LAMA_LOG_INFO( logger,
                   *comm << ": rank = (" << gridRank[0] << "," << gridRank[1] << "," << gridRank[2] << ") of (" << gridSize[0] << "," << gridSize[1] << "," << gridSize[2] << "), local range = [" << dimLB[0] << ":" << dimUB[0] << "," << dimLB[1] << ":" << dimUB[1] << "," << dimLB[2] << ":" << dimUB[2] << "] of " << dimX << " x " << dimY << " x " << dimZ )

    IndexType globalSize = dimX * dimY * dimZ; // number of rows, columns of full matrix

    std::vector<IndexType> myGlobalIndexes; // row indexes of this processor
    std::vector<IndexType> myIA; // number of entries in my rows
    IndexType myNNA = 0; // accumulated sum for number of my entries

    // compute local size on this processor

    IndexType localSize = 1;

    for ( int i = 0; i < dimension; i++ )
    {
        // avoid negative sizes, can happen #procs >> #ndim
        if ( dimLB[i] > dimUB[i] )
        {
            localSize = 0;
            break;
        }
        localSize *= ( dimUB[i] - dimLB[i] + 1 );
    }

    LAMA_LOG_DEBUG( logger, *comm << ": has local size = " << localSize )

    myGlobalIndexes.reserve( localSize );
    myIA.reserve( localSize );

    IndexType length;
    IndexType dimStencil;
    IndexType maxDistance;

    getStencilProperties( dimStencil, length, maxDistance, stencilType );

    LAMA_ASSERT_EQUAL_ERROR( dimStencil, dimension )

    LAMA_LOG_INFO( logger,
                   "stencil type = " << stencilType << " -> dim = " << dimStencil << ", direction length = " << length << ", max distance = " << maxDistance )

    // compute global indexes this processor is responsibile for and number of non-zero values
    // of the rows owned by this processor

    for ( IndexType idZ = dimLB[2]; idZ <= dimUB[2]; ++idZ )
    {
        for ( IndexType idY = dimLB[1]; idY <= dimUB[1]; ++idY )
        {
            for ( IndexType idX = dimLB[0]; idX <= dimUB[0]; ++idX )
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

    LAMA_LOG_INFO( logger, *comm << ": has local " << localSize << " rows, nna = " << myNNA )
    // allocate and fill local part of the distributed matrix

    lama::DistributionPtr distribution = lama::DistributionPtr(
            new lama::GeneralDistribution( globalSize, myGlobalIndexes, comm ) );

    LAMA_LOG_INFO( logger, "distribution = " << *distribution )

    // create new local CSR data ( # local rows x # columns )

    lama::CSRStorage<ValueType> localMatrix;
    localMatrix.allocate( localSize, globalSize );

    // Allocate local matrix with correct sizes and correct first touch in case of OpenMP
    // ToDo: localMatrix( localSize, numColumns, numNonZeros, &myIA[0] );

    lama::LAMAArray<IndexType> csrIA;
    lama::LAMAArray<IndexType> csrJA;
    lama::LAMAArray<ValueType> csrValues;

    {
        lama::HostWriteAccess<IndexType> ia( csrIA );
        lama::HostWriteAccess<IndexType> ja( csrJA );
        lama::HostWriteAccess<ValueType> values( csrValues );

        ia.resize( localSize + 1 );
        ja.resize( myNNA );
        values.resize( myNNA );

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

        for ( IndexType idZ = dimLB[2]; idZ <= dimUB[2]; ++idZ )
        {
            for ( IndexType idY = dimLB[1]; idY <= dimUB[1]; ++idY )
            {
                for ( IndexType idX = dimLB[0]; idX <= dimUB[0]; ++idX )
                {
                    // get column positions and values of matrix, diagonal element is first

                    getStencil( colIndexes, colValues, idX, idY, idZ, dimX, dimY, dimZ, stencilType, length,
                                maxDistance );

                    // check colIndexes.size() against number of values given by getNStencilValues

                    LAMA_ASSERT_EQUAL_DEBUG(
                        (int) colIndexes.size(),
                        getNStencilValues( idX, idY, idZ, dimX, dimY, dimZ, length, maxDistance ) );

                    LAMA_ASSERT_EQUAL_DEBUG( colIndexes.size(), colValues.size() )

                    LAMA_LOG_TRACE( logger,
                                    *comm << ": at " << idX << " x " << idY << " x " << idZ << ", local row : " << rowCounter << ", global row : " << getMatrixPosition( idX, idY, idZ, dimX, dimY, dimZ ) << ": " << colIndexes.size() << " entries: " << colIndexes[0] << ": " << colValues[0] << ", ..." )

                    ia[rowCounter + 1] = ia[rowCounter] + static_cast<IndexType>( colIndexes.size() );

                    for ( size_t k = 0; k < colIndexes.size(); ++k )
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

    LAMA_LOG_DEBUG( logger, "replace owned data with " << localMatrix )

    matrix.assign( localMatrix, distribution, distribution ); // builds also halo

    // but now the local part of matrixA should have the diagonal property as global column
    // indexes have been localized

    LAMA_ASSERT_DEBUG( matrix.getLocalStorage().hasDiagonalProperty(), "CSR data has not diagonal property" )

    LAMA_LOG_INFO( logger, "built matrix A = " << matrix )
}

/* ------------------------------------------------------------------------- */

template<typename T>
void MatrixCreator<T>::buildPoisson1D( CSRSparseMatrix<T>& matrix, const IndexType stencilType, const IndexType dim )
{
    LAMA_LOG_INFO( logger, "build Poisson1D" << stencilType << "P( " << dim << ")" )

    buildPoisson( matrix, 1, stencilType, dim, 1, 1 );
}

/* ------------------------------------------------------------------------- */

template<typename T>
void MatrixCreator<T>::buildPoisson2D(
    CSRSparseMatrix<T>& matrix,
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2 )
{
    LAMA_LOG_INFO( logger, "build Poisson2D" << stencilType << "P( " << dim1 << ", " << dim2 << ")" )

    buildPoisson( matrix, 2, stencilType, dim1, dim2, 1 );
}

/* ------------------------------------------------------------------------- */

template<typename T>
void MatrixCreator<T>::buildPoisson3D(
    CSRSparseMatrix<T>& matrix,
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2,
    const IndexType dim3 )
{
    LAMA_LOG_INFO( logger, "build Poisson3D" << stencilType << "P( " << dim1 << ", " << dim2 << ", " << dim3 << ")" )

    buildPoisson( matrix, 3, stencilType, dim1, dim2, dim3 );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixCreator<ValueType>::fillRandom( Matrix& matrix, double density )
{
    int seed;

    timeval time;
    gettimeofday(&time, NULL);
    seed = ((unsigned long long)time.tv_sec * 1000000) + time.tv_usec;

    // Shape and distribution of matrix is not changed

    const Distribution& dist = matrix.getDistribution();

    if ( dist.getNumPartitions() == 1 )
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

    const IndexType expectedEntries = static_cast<IndexType>( localRowSize * colSize * density );

    std::vector<IndexType> csrIA( localRowSize + 1 );
    std::vector<IndexType> csrJA;
    std::vector<ValueType> csrValues;

    csrJA.reserve( expectedEntries );
    csrValues.reserve( expectedEntries );

    IndexType numValues = 0;

    csrIA[0] = numValues;


    IndexType rowDensity = (IndexType)(localRowSize * density);


//#pragma omp parallel for
    for ( int i = 0; i < localRowSize; ++i )
    {
        IndexType numValuesRow = 0;
        while ( numValuesRow < rowDensity )
        {
            IndexType randJ = rand() % colSize;

            for ( IndexType j = 0; j < numValuesRow; ++j )
            {
                if ( csrJA[i*rowDensity+j] == randJ )
                {
                    continue;
                }
            }
            ValueType value = static_cast<ValueType>( rand() ) / static_cast<ValueType>( RAND_MAX );
            csrJA.push_back( randJ );
            csrValues.push_back( value );
            ++numValuesRow;
        }
        csrIA[i + 1] = (i+1)*rowDensity;
    }



    CSRStorage<ValueType> localCSRStorage;

    localCSRStorage.setRawCSRData( localRowSize, colSize, csrIA[localRowSize], &csrIA[0], &csrJA[0], &csrValues[0] );

    LAMA_LOG_DEBUG( logger, "replace owned data with " << localCSRStorage )

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
    CommunicatorPtr comm = lama::CommunicatorFactory::get( "MPI" );

    DistributionPtr dist( new BlockDistribution( size, comm ) );
    matrix.allocate( dist, dist );
    fillRandom( matrix, density );
}

/* ------------------------------------------------------------------------- */

// Template instanstiations for float + double
template class LAMA_DLL_IMPORTEXPORT MatrixCreator<float> ;
template class LAMA_DLL_IMPORTEXPORT MatrixCreator<double> ;

} // namespace lama

