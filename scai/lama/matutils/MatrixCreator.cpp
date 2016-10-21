/**
 * @file MatrixCreator.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief MatrixCreator.cpp
 * @author Thomas Brandes
 * @date 01.02.2012
 */

// hpp
#include <scai/lama/matutils/MatrixCreator.hpp>

// local library
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

// internal scai libraries
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/common/macros/loop.hpp>

// std
#include <cmath>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( MatrixCreator::logger, "MatrixCreator" )

int MatrixCreator::mSeed = 15191;

static inline void getStencilProperties(
    IndexType& dimension,
    IndexType& length,
    IndexType& distance,
    IndexType stencilType )
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
            COMMON_THROWEXCEPTION( "Invalid Stencil type = " << stencilType )
    }
}

bool MatrixCreator::supportedStencilType( const IndexType dimension, const IndexType stencilType )
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
    IndexType /* dimZ */ )
{
    return idZ * dimX * dimY + idY * dimX + idX;
}

template<typename IndexType>
static inline int getNumNeighbors( IndexType id, IndexType dim, int shift )
{
    int neighbors = common::Math::abs( shift );

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

    if ( maxDistance >= 2 )
    {
        SCAI_ASSERT_ERROR( length == 1, "length > 1 for stencil not supported" )
        numValues += nX * nY + nX * nZ + nY * nZ;
    }

    if ( maxDistance >= 3 )
    {
        SCAI_ASSERT_ERROR( length == 1, "length > 1 for stencil not supported" )
        // dimension cannot be 1, 2
        numValues += nX * nY * nZ;
    }

    // std::cout << "getNStencilValues( idX = " << idX << ", idY = " << idY << ", idZ = " << idZ
    //                         << ", dimX = " << dimX << ", dimY = " << dimY << ", dimZ = " << dimZ 
    //                         << ", len = " << length << ", maxDist = " << maxDistance  << " ) -> " << numValues << std::endl;

    return numValues;
}

static inline void getStencil(
    std::vector<IndexType>& positions,
    std::vector<int>& values,
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
    values.push_back( stencilType - 1 );

    int leftX = getNumNeighbors( idX, dimX, -length );
    int rightX = getNumNeighbors( idX, dimX, length );
    int leftY = getNumNeighbors( idY, dimY, -length );
    int rightY = getNumNeighbors( idY, dimY, length );
    int leftZ = getNumNeighbors( idZ, dimZ, -length );
    int rightZ = getNumNeighbors( idZ, dimZ, length );

    // std::cout << "range [ " << leftX << ":" << rightX << ", " << leftY << ":" << rightY << ", " << leftZ << ":" << rightZ << "]"
    //           << " for " << idX << ":" << idY << ":" << idZ << std::endl;

    for ( int dz = -leftZ; dz <= rightZ; ++dz )
    {
        IndexType distZ = common::Math::abs( dz );

        for ( int dy = -leftY; dy <= rightY; ++dy )
        {
            IndexType distYZ = distZ + common::Math::abs( dy );

            if ( distYZ > maxDistance )
            {
                continue; // can already skip next loop
            }

            for ( int dx = -leftX; dx <= rightX; ++dx )
            {
                IndexType distXYZ = distYZ + common::Math::abs( dx);

                if ( distXYZ > maxDistance )
                {
                    continue;
                }

                // skip the diagonal element already added

                if ( distXYZ == 0 )
                {
                    continue;
                }

                positions.push_back( getMatrixPosition( idX + dx, idY + dy, idZ + dz, dimX, dimY, dimZ ) );
                values.push_back( -1 );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildPoisson(
    Matrix& matrix,
    const IndexType dimension,
    const IndexType stencilType,
    const IndexType dimX,
    const IndexType dimY,
    const IndexType dimZ )
{
    // Calculate subdomains, subranges
    PartitionId gridSize[3] = { 1, 1, 1 };
    PartitionId gridRank[3] = { 0, 0, 0 };
    PartitionId dimLB[3]    = { 0, 0, 0 };
    PartitionId dimUB[3]    = { dimX, dimY, dimZ };

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr( );

    // get rank of this processor

    if ( dimension == 1 )
    {
        gridSize[0] = comm->getSize();
        gridRank[0] = comm->getRank();
        dmemo::BlockDistribution::getLocalRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
    }
    else if ( dimension == 2 )
    {
        comm->factorize2( gridSize, static_cast<double>( dimX ), static_cast<double>( dimY ) );
        comm->getGrid2Rank( gridRank, gridSize );
        dmemo::BlockDistribution::getLocalRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
        dmemo::BlockDistribution::getLocalRange( dimLB[1], dimUB[1], dimY, gridRank[1], gridSize[1] );
    }
    else if ( dimension == 3 )
    {
        comm->factorize3( gridSize, static_cast<double>( dimX ), static_cast<double>( dimY ), static_cast<double>( dimZ ) );
        comm->getGrid3Rank( gridRank, gridSize );
        dmemo::BlockDistribution::getLocalRange( dimLB[0], dimUB[0], dimX, gridRank[0], gridSize[0] );
        dmemo::BlockDistribution::getLocalRange( dimLB[1], dimUB[1], dimY, gridRank[1], gridSize[1] );
        dmemo::BlockDistribution::getLocalRange( dimLB[2], dimUB[2], dimZ, gridRank[2], gridSize[2] );
    }

    SCAI_LOG_INFO( logger,
                   *comm << ": rank = (" << gridRank[0] << "," << gridRank[1] << "," << gridRank[2] 
                    << ") of (" << gridSize[0] << "," << gridSize[1] << "," << gridSize[2] 
                    << "), local range = [" << dimLB[0] << ":" << dimUB[0] << "," << dimLB[1] << ":" << dimUB[1] << "," << dimLB[2] << ":" << dimUB[2] 
                    << "] of " << dimX << " x " << dimY << " x " << dimZ )

    IndexType globalSize = dimX * dimY * dimZ; // number of rows, columns of full matrix
    std::vector<IndexType> myGlobalIndexes; // row indexes of this processor
    std::vector<IndexType> myIA; // number of entries in my rows
    IndexType myNNA = 0; // accumulated sum for number of my entries
    // compute local size on this processor
    IndexType localSize = 1;

    for ( IndexType i = 0; i < dimension; i++ )
    {
        localSize *= dimUB[i] - dimLB[i];
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

    for ( IndexType idZ = dimLB[2]; idZ < dimUB[2]; ++idZ )
    {
        for ( IndexType idY = dimLB[1]; idY < dimUB[1]; ++idY )
        {
            for ( IndexType idX = dimLB[0]; idX < dimUB[0]; ++idX )
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
    hmemo::HArrayRef<IndexType> indexes(  myGlobalIndexes );
    dmemo::DistributionPtr distribution( new dmemo::GeneralDistribution( globalSize, indexes, comm ) );
    SCAI_LOG_INFO( logger, "distribution = " << *distribution )
    // create new local CSR data ( # local rows x # columns )
    scai::lama::CSRStorage<double> localMatrix;
    localMatrix.allocate( localSize, globalSize );
    // Allocate local matrix with correct sizes and correct first touch in case of OpenMP
    // ToDo: localMatrix( localSize, numColumns, numNonZeros, &myIA[0] );
    hmemo::HArray<IndexType> csrIA;
    hmemo::HArray<IndexType> csrJA;
    hmemo::HArray<double> csrValues;
    {
        hmemo::WriteOnlyAccess<IndexType> ia( csrIA, localSize + 1 );
        hmemo::WriteOnlyAccess<IndexType> ja( csrJA, myNNA );
        hmemo::WriteOnlyAccess<double> values( csrValues, myNNA );
        ia[0] = 0;

        std::vector<IndexType> colIndexes;
        std::vector<int> colValues;     

        colIndexes.reserve( stencilType );
        colValues.reserve( stencilType );
        // compute global indexes this processor is responsibile for and number of non-zero values
        // Important: same loop order as above
        IndexType rowCounter = 0; // count local rows
        IndexType nnzCounter = 0; // count local non-zero elements

        for ( IndexType idZ = dimLB[2]; idZ < dimUB[2]; ++idZ )
        {
            for ( IndexType idY = dimLB[1]; idY < dimUB[1]; ++idY )
            {
                for ( IndexType idX = dimLB[0]; idX < dimUB[0]; ++idX )
                {
                    // get column positions and values of matrix, diagonal element is first
                    getStencil( colIndexes, colValues, idX, idY, idZ, dimX, dimY, dimZ, stencilType, length,
                                maxDistance );
                    // check colIndexes.size() against number of values given by getNStencilValues
                    SCAI_ASSERT_EQUAL_DEBUG( ( IndexType ) colIndexes.size(),
                                             getNStencilValues( idX, idY, idZ, dimX, dimY, dimZ, length, maxDistance ) );
                    SCAI_ASSERT_EQUAL_DEBUG( colIndexes.size(), colValues.size() )
                    SCAI_LOG_TRACE( logger,
                                    *comm << ": at " << idX << " x " << idY << " x " << idZ << ", local row : " << rowCounter << ", global row : " << getMatrixPosition( idX, idY, idZ, dimX, dimY, dimZ ) << ": " << colIndexes.size() << " entries: " << colIndexes[0] << ": " << colValues[0] << ", ..." )
                    ia[rowCounter + 1] = ia[rowCounter] + static_cast<IndexType>( colIndexes.size() );

                    for ( size_t k = 0; k < colIndexes.size(); ++k )
                    {
                        ja[nnzCounter] = colIndexes[k];
                        values[nnzCounter] = static_cast<double>( colValues[k] );
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

    // but now the local part of matrixA should have the diagonal property as global column // indexes have been localized
    // is not for each storage format the case 
    // SCAI_ASSERT_DEBUG( matrix.getLocalStorage().hasDiagonalProperty(), "local storage data has not diagonal property: " << matrix )

    SCAI_LOG_INFO( logger, "built matrix A = " << matrix )
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildPoisson1D(
    Matrix& matrix,
    const IndexType stencilType,
    const IndexType dim )
{
    SCAI_LOG_INFO( logger, "build Poisson1D" << stencilType << "P( " << dim << ")" )
    buildPoisson( matrix, 1, stencilType, dim, 1, 1 );
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildPoisson2D(
    Matrix& matrix,
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2 )
{
    SCAI_LOG_INFO( logger, "build Poisson2D" << stencilType << "P( " << dim1 << ", " << dim2 << ")" )
    buildPoisson( matrix, 2, stencilType, dim1, dim2, 1 );
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildPoisson3D(
    Matrix& matrix,
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2,
    const IndexType dim3 )
{
    SCAI_LOG_INFO( logger, "build Poisson3D" << stencilType << "P( " << dim1 << ", " << dim2 << ", " << dim3 << ")" )
    buildPoisson( matrix, 3, stencilType, dim1, dim2, dim3 );
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::fillRandom( Matrix& matrix, double density )
{
    // Shape and distribution of matrix is not changed
    const dmemo::Distribution& dist = matrix.getRowDistribution();

    if ( dist.getNumPartitions() == 1 )
    {
        // all processors must have the same random numbers
        std::srand( mSeed );
    }
    else
    {
        // processors must generate different numbers
        PartitionId rank = dist.getCommunicator().getRank();
        std::srand( mSeed + rank );
    }

    const IndexType localRowSize = dist.getLocalSize();

    const IndexType colSize = matrix.getNumColumns();
    const IndexType expectedEntries = static_cast<IndexType>( localRowSize * colSize * density + 30.0 );

    std::vector<IndexType> csrIA( localRowSize + 1 );
    std::vector<IndexType> csrJA;   // take here a vector, more convenient for push_back

    csrJA.reserve( expectedEntries );

    IndexType numValues = 0;

    csrIA[0] = numValues;

    for ( IndexType i = 0; i < localRowSize; ++i )
    {
        for ( IndexType j = 0; j < colSize; ++j )
        {
            double value = static_cast<double>( rand() ) / static_cast<double>( RAND_MAX );

            if ( value < density )
            {
                csrJA.push_back( j );
                ++numValues;
            }
        }

        csrIA[i + 1] = numValues;
    }

    // now we draw the non-zero values

    common::unique_ptr<hmemo::_HArray> csrValues( hmemo::_HArray::create( matrix.getValueType() ) );

    hmemo::_HArray& values = *csrValues;

    // draw the non-zero values, now with fill rate 1.0f

    utilskernel::HArrayUtils::setRandom( values, numValues, 1.0f );
    
    // some tricky stuff to avoid an additional copy

    _MatrixStorage& localMatrix = const_cast<_MatrixStorage&>( matrix.getLocalStorage() );

    hmemo::HArrayRef<IndexType> ia( csrIA );
    hmemo::HArrayRef<IndexType> ja( csrJA );

    localMatrix.setCSRData( localRowSize, colSize, numValues, ia, ja, *csrValues );

    SCAI_LOG_DEBUG( logger, "local random part: " << localMatrix )

    // The new matrix data has the same row distribution as the input
    // matrix, also take over the original column distribution to build halo

    matrix.assign( localMatrix, matrix.getRowDistributionPtr(), matrix.getColDistributionPtr() ); 
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildRandom(
    Matrix& matrix,
    const IndexType size,
    const double density )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr( );
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( size, comm ) );
    matrix.allocate( dist, dist );
    fillRandom( matrix, density );
}

/* ------------------------------------------------------------------------- */

/** Help routine to replicate matrix storage in diagonal blocks. */

template<typename ValueType>
static void replicateStorageDiag( 
    MatrixStorage<ValueType>& out,
    const MatrixStorage<ValueType>& in,
    const IndexType nRepeat )
{
    using namespace hmemo;

    // replication is done via CSR storage data, so get it form input storage

    HArray<IndexType> inIA;
    HArray<IndexType> inJA;
    HArray<ValueType> inValues;

    in.buildCSRData( inIA, inJA, inValues );

    IndexType nRows   = inIA.size() - 1;
    IndexType nCols   = in.getNumColumns();
    IndexType nValues = inJA.size();

    HArray<IndexType> outIA;
    HArray<IndexType> outJA;
    HArray<ValueType> outValues;

    // Replication is done at the host

    ContextPtr ctx = Context::getHostPtr();

    {
        IndexType offset = 0;   // current offset in new JA, Values array

        ReadAccess<IndexType> rIA( inIA, ctx );
        ReadAccess<IndexType> rJA( inJA, ctx );
        ReadAccess<ValueType> rValues( inValues, ctx );

        WriteOnlyAccess<IndexType> wIA( outIA, ctx, nRows * nRepeat + 1 );
        WriteOnlyAccess<IndexType> wJA( outJA, ctx, nValues * nRepeat );
        WriteOnlyAccess<ValueType> wValues( outValues, ctx, nValues * nRepeat );

        for ( IndexType iRepeat = 0; iRepeat < nRepeat; iRepeat++ )
        {
            for ( IndexType i = 0; i < nRows; i++ )
            {
                wIA[ iRepeat * nRows + i ] = offset + rIA[i];   // the current row offset
            }

            IndexType colOffset = iRepeat * nCols;

            for ( IndexType k = 0; k < nValues; ++k )
            {
                wJA[ offset ] = rJA[ k ] + colOffset;
                wValues[ offset ] = rValues[ k ];
                offset++;
            }
        }

        wIA[ nRepeat * nRows ] = offset;  // final offset

        SCAI_ASSERT_EQUAL( offset, nRepeat * nValues, "size mismatch" );
    }

    out.setCSRData( nRows * nRepeat, nCols * nRepeat, nValues * nRepeat,
                    outIA, outJA, outValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
static void replicateStorage( 
    MatrixStorage<ValueType>& out,
    const MatrixStorage<ValueType>& in,
    const IndexType nRepeatRow,
    const IndexType nRepeatCol )
{
    IndexType nRepeat = nRepeatRow * nRepeatCol;

    using namespace hmemo;
    HArray<IndexType> inIA;
    HArray<IndexType> inJA;
    HArray<ValueType> inValues;

    in.buildCSRData( inIA, inJA, inValues );

    IndexType nRows   = inIA.size() - 1;
    IndexType nCols   = in.getNumColumns();
    IndexType nValues = inJA.size();

    HArray<IndexType> outIA;
    HArray<IndexType> outJA;
    HArray<ValueType> outValues;

    // Replication is done at the host

    ContextPtr ctx = Context::getHostPtr();

    {
        IndexType offset = 0;   // current offset in new JA, Values array

        ReadAccess<IndexType> rIA( inIA, ctx );
        ReadAccess<IndexType> rJA( inJA, ctx );
        ReadAccess<ValueType> rValues( inValues, ctx );

        WriteOnlyAccess<IndexType> wIA( outIA, ctx, nRows * nRepeatRow + 1 );
        WriteOnlyAccess<IndexType> wJA( outJA, ctx, nValues * nRepeat );
        WriteOnlyAccess<ValueType> wValues( outValues, ctx, nValues * nRepeat );

        for ( IndexType iRepeat = 0; iRepeat < nRepeatRow; iRepeat++ )
        {
            for ( IndexType i = 0; i < nRows; i++ )
            {
                wIA[ iRepeat * nRows + i ] = offset;   // the current row offset

                IndexType rowOffset  = rIA[i];
                IndexType nRowValues = rIA[i+1] - rIA[i];

                IndexType colOffset = 0;

                for ( IndexType jRepeat = 0; jRepeat < nRepeatCol; jRepeat++ )
                {
                    for ( IndexType jj = 0; jj < nRowValues; ++jj )
                    {
                        wJA[ offset ] = rJA[ rowOffset + jj ] + colOffset;
                        wValues[ offset ] = rValues[ rowOffset + jj ];
                        offset++;
                    }

                    colOffset += nCols;
                }
            }
        }

        wIA[ nRepeatRow * nRows ] = offset;   // the final offset

        SCAI_ASSERT_EQUAL( offset, nRepeat * nValues, "serious offset/size mismatch" )
    }

    out.setCSRData( nRows * nRepeatRow, nCols * nRepeatCol, nValues * nRepeat,
                    outIA, outJA, outValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixCreator::buildReplicatedDiag( 
    SparseMatrix<ValueType>& matrix,
    const MatrixStorage<ValueType>& storage,
    const IndexType nRepeat )
{
    using namespace dmemo;

    // create the distribution

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr( );

    IndexType nRows     = storage.getNumRows() * nRepeat;
    IndexType nCols     = storage.getNumColumns() * nRepeat;

    IndexType nChunks;  // will be number of chunks for this processor

    // this replication will never split any of the storages 
    // bit tricky: use a cyclic( 1 ) distribution of nrepeat to get the local size of this processor

    {
        CyclicDistribution cdist( nRepeat, 1, comm );
        nChunks = cdist.getNumLocalChunks();
    }

    // we will take a general block distribution for the rows

    dmemo::DistributionPtr rowDist( new dmemo::GenBlockDistribution( nRows, nChunks * storage.getNumRows(), comm ) );

    // we will take also a general block distribution for the columns to avoid the translation into global indexes

    dmemo::DistributionPtr colDist( new dmemo::GenBlockDistribution( nCols, nChunks * storage.getNumColumns(), comm ) );

    SCAI_LOG_DEBUG( logger, *comm << ": row dist = " << *rowDist )
    SCAI_LOG_DEBUG( logger, *comm << ": col dist = " << *colDist )

    // Allocate the correct size

    matrix.allocate( rowDist, colDist );

    // now build the local part

    if ( nChunks >= 1 )
    {
        const MatrixStorage<ValueType>& local = matrix.getLocalStorage();

        // we only change the local part, this is safe!

        replicateStorageDiag( const_cast<MatrixStorage<ValueType>& >( local ), storage, nChunks );
    }

    // The halo part remains empty, so we are done
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixCreator::buildReplicated( SparseMatrix<ValueType>& matrix,
        const MatrixStorage<ValueType>& storage,
        const IndexType nRepeatRow,
        const IndexType nRepeatCol )
{
    using namespace dmemo;

    // create the distribution

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr( );

    IndexType chunkSize = storage.getNumRows();
    IndexType nGlobal   = storage.getNumRows() * nRepeatRow;

    IndexType nLocal;   // will be #local chunks * chunkSize
    IndexType nChunks;  // will be #local chunks * chunkSi

    // bit tricky: use a cyclic( nLocal)  distribution to get the local size of this processor

    {
        CyclicDistribution cdist( nGlobal, chunkSize, comm );
        nChunks = cdist.getNumLocalChunks();
        nLocal  = cdist.getLocalSize();

        SCAI_ASSERT_EQUAL( nLocal, nChunks * chunkSize, "serious mismatch" )
    }

    // we will take a general block distribution

    dmemo::DistributionPtr rowDist( new dmemo::GenBlockDistribution( nGlobal, nLocal, comm ) );
    dmemo::DistributionPtr colDist( new dmemo::NoDistribution( storage.getNumColumns() * nRepeatCol  ) );

    SCAI_LOG_DEBUG( logger, *comm << ": rowDist for replicated matrix = " << *rowDist )

    // Allocate the correct size

    matrix.allocate( rowDist, colDist );

    // now build the local part

    if ( nChunks >= 1 )
    {
        const MatrixStorage<ValueType>& local = matrix.getLocalStorage();

        // we only change the local part, this is safe!

        replicateStorage( const_cast<MatrixStorage<ValueType>& >( local ), storage, nChunks, nRepeatCol );
    }
}

/* ========================================================================= */
/*       Template Method Instantiations                                      */
/* ========================================================================= */

#define MATRIX_CREATOR_SPECIFIER( ValueType )                                                  \
    template void MatrixCreator::buildReplicated( SparseMatrix<ValueType>& matrix,             \
                                                  const MatrixStorage<ValueType>& storage,     \
                                                  const IndexType nRepeatRow,                  \
                                                  const IndexType nRepeatCol );                \
    template void MatrixCreator::buildReplicatedDiag( SparseMatrix<ValueType>& matrix,         \
                                                      const MatrixStorage<ValueType>& storage, \
                                                      const IndexType nRepeat ) ;

SCAI_COMMON_LOOP( MATRIX_CREATOR_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef MATRIX_CREATOR_SPECIFIER

} /* end namespace lama */

} /* end namespace scai */

