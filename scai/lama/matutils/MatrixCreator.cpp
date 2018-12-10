/**
 * @file MatrixCreator.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief MatrixCreator.cpp
 * @author Thomas Brandes
 * @date 01.02.2012
 */

// hpp
#include <scai/lama/matutils/MatrixCreator.hpp>

// local library
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
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

/* ------------------------------------------------------------------------- */ 

template<typename IndexType>
static inline IndexType getNStencilValues(
    const IndexType pos[],
    const common::Grid& grid,
    IndexType length,
    IndexType maxDistance )
{
    SCAI_ASSERT_DEBUG( grid.validPos( pos ), "Invalid pos = " << pos[0] << "-" << pos[1] << "-" << pos[2] << " in " << grid )

    IndexType nX = getNumNeighbors( pos[0], grid.size( 0 ), -length ) + getNumNeighbors( pos[0], grid.size( 0 ), length );
    IndexType nY = getNumNeighbors( pos[1], grid.size( 1 ), -length ) + getNumNeighbors( pos[1], grid.size( 1 ), length );
    IndexType nZ = getNumNeighbors( pos[2], grid.size( 2 ), -length ) + getNumNeighbors( pos[2], grid.size( 2 ), length );

    // Note: nZ will be 0 for dim <= 2, nY will be 0 for dim <= 1

    IndexType numValues = 1 + nX + nY + nZ; // maxDistance 1

    if ( maxDistance >= 2 )
    {
        SCAI_ASSERT_EQ_ERROR( 1, length, "length > 1 for stencil not supported" )
        numValues += nX * nY + nX * nZ + nY * nZ;
    }

    if ( maxDistance >= 3 )
    {
        SCAI_ASSERT_EQ_ERROR( 1, length, "length > 1 for stencil not supported" )

        // dimension cannot be 1, 2

        numValues += nX * nY * nZ;
    }

    /*
    std::cout << "getNStencilValues( pos = " << pos[0] << ", " << pos[1] << ", " << pos[2]
                            << ", grid = " << grid
                            << ", len = " << length << ", maxDist = " << maxDistance  << " ) -> " << numValues << std::endl;
    */

    return numValues;
}

static inline void getStencil(
    std::vector<IndexType>& positions,
    std::vector<int>& values,
    const IndexType pos[],
    const common::Grid& grid,
    const IndexType stencilType,
    const IndexType length,
    const IndexType maxDistance )
{
    positions.clear(); // reset number of entries to 0
    values.clear(); // reset

    positions.push_back( grid.linearPos( pos ) );
    values.push_back( stencilType - 1 );

    int leftX = getNumNeighbors( pos[0], grid.size( 0 ), -length );
    int rightX = getNumNeighbors( pos[0], grid.size( 0 ), length );
    int leftY = getNumNeighbors( pos[1], grid.size( 1 ), -length );
    int rightY = getNumNeighbors( pos[1], grid.size( 1 ), length );
    int leftZ = getNumNeighbors( pos[2], grid.size( 2 ), -length );
    int rightZ = getNumNeighbors( pos[2], grid.size( 2 ), length );

    // std::cout << "range [ " << leftX << ":" << rightX << ", " << leftY << ":" << rightY << ", " << leftZ << ":" << rightZ << "]"
    //           << " for " << pos[0] << ":" << pos[1] << ":" << pos[2] << std::endl;

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
                IndexType distXYZ = distYZ + common::Math::abs( dx );

                if ( distXYZ > maxDistance )
                {
                    continue;
                }

                // skip the diagonal element already added

                if ( distXYZ == 0 )
                {
                    continue;
                }

                // reinterpretation as 3D grid is safe as all 3 dimensions are set

                const common::Grid3D& grid3D = reinterpret_cast<const common::Grid3D&>( grid );

                positions.push_back( grid3D.linearPos( pos[0] + dx, pos[1] + dy, pos[2] + dz ) );

                values.push_back( -1 );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildPoisson(
    _Matrix& matrix,
    const IndexType dimension,
    const IndexType stencilType,
    const IndexType dimX,
    const IndexType dimY,
    const IndexType dimZ )
{
    using namespace dmemo;   // simplifies working with distributons, communicator

    IndexType sizes[] = { dimX, dimY, dimZ };

    common::Grid globalGrid( dimension, sizes );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr( );   // get the default communicator

    std::shared_ptr<GridDistribution> gridDistribution( new GridDistribution( globalGrid, comm ) );

    SCAI_LOG_INFO( logger, *comm << ": " << *gridDistribution )

    IndexType myNNA = 0; // accumulated sum for number of my entries

    IndexType length;
    IndexType dimStencil;
    IndexType maxDistance;

    getStencilProperties( dimStencil, length, maxDistance, stencilType );

    SCAI_ASSERT_EQUAL_ERROR( dimStencil, dimension )

    SCAI_LOG_INFO( logger, "stencil type = " << stencilType << " -> dim = " << dimStencil 
                            << ", direction length = " << length << ", max distance = " << maxDistance )

    // Loop over all local grid points and count number of matrix values needed
    // Actual we could make a guess here with getLocalSize * nPoints

    const common::Grid& localGrid = gridDistribution->getLocalGrid();
    const IndexType localSize = gridDistribution->getLocalSize();

    IndexType localGridPos[3]  = { 0, 0, 0 };
    IndexType globalGridPos[3] = { 0, 0, 0 };

    for ( IndexType i = 0; i < localSize; ++i )
    {
        localGrid.gridPos( localGridPos, i );
        gridDistribution->local2global( globalGridPos, localGridPos );

        const IndexType numNonZeros = getNStencilValues( globalGridPos, globalGrid, length, maxDistance );

        myNNA += numNonZeros;
    }

    // create new local CSR data ( # local rows x # columns )

    CSRStorage<DefaultReal> localMatrix;

    localMatrix.allocate( localSize, globalGrid.size() );

    hmemo::HArray<IndexType> csrIA;
    hmemo::HArray<IndexType> csrJA;
    hmemo::HArray<DefaultReal> csrValues;

    {
        hmemo::WriteOnlyAccess<IndexType> ia( csrIA, localSize + 1 );
        hmemo::WriteOnlyAccess<IndexType> ja( csrJA, myNNA );
        hmemo::WriteOnlyAccess<DefaultReal> values( csrValues, myNNA );

        ia[0] = 0;

        std::vector<IndexType> colIndexes;
        std::vector<int> colValues;

        colIndexes.reserve( stencilType );
        colValues.reserve( stencilType );

        // compute global indexes this processor is responsibile for and number of non-zero values
        // Important: same loop order as above
        IndexType rowCounter = 0; // count local rows
        IndexType nnzCounter = 0; // count local non-zero elements

        for ( IndexType i = 0; i < localSize; ++i )
        {
            localGrid.gridPos( localGridPos, i );
            gridDistribution->local2global( globalGridPos, localGridPos );

            // get column positions and values of matrix, diagonal element is first

            getStencil( colIndexes, colValues, globalGridPos, globalGrid, stencilType, length, maxDistance );

            // check colIndexes.size() against number of values given by getNStencilValues

            SCAI_ASSERT_EQ_DEBUG( ( IndexType ) colIndexes.size(), 
                                    getNStencilValues( globalGridPos, globalGrid, length, maxDistance ), 
                                    "serious mismatch" );

            SCAI_ASSERT_EQ_DEBUG( colIndexes.size(), colValues.size(), "serious mismatch" )

            SCAI_LOG_DEBUG( logger,
                            *comm << ": at " << globalGridPos[0] << "-" << globalGridPos[1] << "-" << globalGridPos[2] 
                             << ", local row : " << rowCounter << ", global row : " << globalGrid.linearPos( globalGridPos ) 
                             << ": " << colIndexes.size() << " entries: " << colIndexes[0] << ": " << colValues[0] << ", ..." )

            ia[rowCounter + 1] = ia[rowCounter] + static_cast<IndexType>( colIndexes.size() );

            for ( size_t k = 0; k < colIndexes.size(); ++k )
            {
                ja[nnzCounter] = colIndexes[k];
                values[nnzCounter] = static_cast<DefaultReal>( colValues[k] );
                ++nnzCounter;
            }

            ++rowCounter;
        }
    }

    localMatrix.swap( csrIA, csrJA, csrValues );
    SCAI_LOG_DEBUG( logger, "replace owned data with " << localMatrix )
    matrix.assignLocal( localMatrix, gridDistribution );   
    matrix.redistribute( gridDistribution, gridDistribution );     // builds also halo

    SCAI_LOG_INFO( logger, "built matrix A = " << matrix )
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildPoisson1D(
    _Matrix& matrix,
    const IndexType stencilType,
    const IndexType dim )
{
    SCAI_LOG_INFO( logger, "build Poisson1D" << stencilType << "P( " << dim << ")" )
    buildPoisson( matrix, 1, stencilType, dim, 1, 1 );
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildPoisson2D(
    _Matrix& matrix,
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2 )
{
    SCAI_LOG_INFO( logger, "build Poisson2D" << stencilType << "P( " << dim1 << ", " << dim2 << ")" )
    buildPoisson( matrix, 2, stencilType, dim1, dim2, 1 );
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildPoisson3D(
    _Matrix& matrix,
    const IndexType stencilType,
    const IndexType dim1,
    const IndexType dim2,
    const IndexType dim3 )
{
    SCAI_LOG_INFO( logger, "build Poisson3D" << stencilType << "P( " << dim1 << ", " << dim2 << ", " << dim3 << ")" )
    buildPoisson( matrix, 3, stencilType, dim1, dim2, dim3 );
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::randomCSRPattern( 
    std::vector<IndexType>& csrIA, 
    std::vector<IndexType>& csrJA, 
    const IndexType numRows, 
    const IndexType numCols, 
    float density )
{
    const IndexType expectedEntries = static_cast<IndexType>( numRows * numCols * density + 30.0 );

    csrIA.resize( numRows + 1 );
    csrJA.reserve( expectedEntries );

    IndexType numValues = 0;

    csrIA[0] = numValues;

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            // density determines the true ratio of random bool value

            bool takeIt = common::Math::randomBool( density );

            if ( takeIt )
            {
                csrJA.push_back( j );
                ++numValues;
            }
        }

        csrIA[i + 1] = numValues;
    }
}

template<typename ValueType>
void MatrixCreator::fillRandomImpl( Matrix<ValueType>& matrix, float density )
{
    dmemo::DistributionPtr saveColDist = matrix.getColDistributionPtr();

    const IndexType localRowSize = matrix.getRowDistribution().getLocalSize();
    const IndexType colSize      = matrix.getColDistribution().getGlobalSize();

    std::vector<IndexType> csrIA;
    std::vector<IndexType> csrJA;   

    randomCSRPattern( csrIA, csrJA, localRowSize, colSize, density );

    SCAI_LOG_INFO( logger, "generated csrIA( size = " << csrIA.size() << " ), csrJA( size = " << csrJA.size() << " )" )
 
    IndexType numValues = csrJA.size();    

    // now we draw the non-zero values

    hmemo::HArray<ValueType> values( numValues );

    utilskernel::HArrayUtils::setRandom( values, 1 );

    // ToDo: use move operations when available

    hmemo::HArrayRef<IndexType> ia( csrIA );
    hmemo::HArrayRef<IndexType> ja( csrJA );

    CSRStorage<ValueType> localStorage( localRowSize, colSize, ia, ja, values );

    // The new matrix data has the same row distribution as the input matrix

    matrix.assignLocal( localStorage, matrix.getRowDistributionPtr() );

    // take over the original column distribution to build halo

    matrix.redistribute( matrix.getRowDistributionPtr(), saveColDist );

    SCAI_LOG_INFO( logger, "matrix: " << matrix )
}

/* ------------------------------------------------------------------------- */

// Metaprogramming to translate assign( _Matrix ) to assignImpl( Matrix<ValueType> )

template<typename TList>
struct RandomWrapper;

template<>
struct RandomWrapper<common::mepr::NullType>
{
    static void fillRandom( _Matrix& matrix, float )
    {
        COMMON_THROWEXCEPTION( "Unsupported matrix value type: --> " << matrix )
    }
};

template<typename H, typename T>
struct RandomWrapper<common::mepr::TypeList<H, T> >
{
    static void fillRandom( _Matrix& matrix, float density )
    {
        if ( matrix.getValueType() == common::getScalarType<H>() )
        {
            MatrixCreator::fillRandomImpl( static_cast<Matrix<H>& >( matrix ), density );
        }
        else
        {
            RandomWrapper<T>::fillRandom( matrix, density );
        }
    }
};

void MatrixCreator::fillRandom( _Matrix& matrix, float density )
{
    // use meta programming to call the version

    RandomWrapper<SCAI_NUMERIC_TYPES_HOST_LIST>::fillRandom( matrix, density );
}

/* ------------------------------------------------------------------------- */

void MatrixCreator::buildRandom(
    _Matrix& matrix,
    const IndexType size,
    const float density )
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

    out.setCSRData( nRows * nRepeat, nCols * nRepeat, outIA, outJA, outValues );
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
                IndexType nRowValues = rIA[i + 1] - rIA[i];

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

    out.setCSRData( nRows * nRepeatRow, nCols * nRepeatCol, outIA, outJA, outValues );
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

    IndexType nChunks;  // will be number of chunks for this processor

    // this replication will never split any of the storages
    // bit tricky: use a cyclic( 1 ) distribution of nrepeat to get the local size of this processor

    {
        CyclicDistribution cdist( nRepeat, 1, comm );
        nChunks = cdist.getNumLocalChunks();
    }

    // we will take a general block distribution for the rows

    dmemo::DistributionPtr rowDist = genBlockDistribution( nChunks * storage.getNumRows(), comm );

    // we will take also a general block distribution for the columns to avoid the translation into global indexes

    dmemo::DistributionPtr colDist = genBlockDistribution( nChunks * storage.getNumColumns(), comm );

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

    dmemo::DistributionPtr rowDist = dmemo::genBlockDistribution( nLocal, comm );
    dmemo::DistributionPtr colDist = dmemo::noDistribution( storage.getNumColumns() * nRepeatCol  );

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

