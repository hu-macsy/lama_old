/**
 * @file GridVector.cpp
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
 * @brief Implementation of methods for grid vector
 * @author Thomas Brandes
 * @date 12.05.2017
 */

#include <scai/lama/GridVector.hpp>

// other SCAI libraries

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/SectionKernelTrait.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/io/MatlabIO.hpp>

#include <scai/common/macros/instantiate.hpp>
#include <scai/common/SCAITypes.hpp>

#include <scai/tracing.hpp>

namespace scai
{

using namespace hmemo;

namespace lama
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, GridVector<ValueType>::logger, "Vector.GridVector" )

/* ------------------------------------------------------------------------- */

template<typename ValueType>
GridVector<ValueType>::GridVector( hmemo::ContextPtr context ) : DenseVector<ValueType>( context )
{
    // give it a default grid distribution for consistency

    this->setDistributionPtr( dmemo::gridDistributionReplicated( common::Grid() ) );
}

template<typename ValueType>
GridVector<ValueType>::GridVector( const common::Grid& grid, const ValueType value, hmemo::ContextPtr context ) :

    DenseVector<ValueType>( dmemo::gridDistributionReplicated( grid ), value, context )
{
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
GridVector<ValueType>::GridVector( dmemo::DistributionPtr dist, const ValueType value, hmemo::ContextPtr context ) :

    DenseVector<ValueType>( dist, value, context )

{
     SCAI_ASSERT_ERROR( dynamic_cast<const dmemo::GridDistribution*>( dist.get() ), 
                        "Grid vector has no grid distribution, dist = " << *dist )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void GridVector<ValueType>::reduce( const GridVector<ValueType>& other, IndexType dim, const common::BinaryOp redOp )
{
    SCAI_ASSERT_VALID_INDEX_ERROR( dim, other.nDims(), "illegeal reduction dim on this grid " << other.globalGrid() )

    COMMON_THROWEXCEPTION( "reduction on grids not supported yet, reduction op = " << redOp )
}

template<typename ValueType>
void GridVector<ValueType>::gemm( const ValueType alpha, const GridVector<ValueType>& v1, const GridVector<ValueType>& v2 )
{
    const common::Grid& resGrid = this->globalGrid();
    const common::Grid& grid1 = v1.globalGrid();
    const common::Grid& grid2 = v2.globalGrid();

    SCAI_ASSERT_EQ_ERROR( 2, resGrid.nDims(), "gemm only on two-dimensional grids." )
    SCAI_ASSERT_EQ_ERROR( 2, grid1.nDims(), "gemm only on two-dimensional grids." )
    SCAI_ASSERT_EQ_ERROR( 2, grid2.nDims(), "gemm only on two-dimensional grids." )

    SCAI_ASSERT_EQ_ERROR( resGrid.size( 0 ), grid1.size( 0 ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( resGrid.size( 1 ), grid2.size( 1 ), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( grid1.size( 1 ), grid2.size( 0 ), "size mismatch" )

    // ToDo: not yet for distributed grids

    const IndexType m = resGrid.size( 0 );
    const IndexType n = resGrid.size( 1 );
    const IndexType k = grid1.size( 1 );

    int lda = grid1.size( 1 );
    int ldb = grid2.size( 1 );
    int ldc = resGrid.size( 1 );

    using common::MatrixOp;

    if ( lda != 0 && n != 0 && m != 0 )
    {
        static utilskernel::LAMAKernel<blaskernel::BLASKernelTrait::gemm<ValueType> > gemm;

        ContextPtr loc = this->getContextPtr();

        gemm.getSupportedContext( loc );

        GridWriteAccess<ValueType> wRes( *this, loc );
        GridReadAccess<ValueType> rA( v1, loc );
        GridReadAccess<ValueType> rB( v2, loc );
  
        SCAI_CONTEXT_ACCESS( loc )

        ValueType beta = 1;

        gemm[loc]( MatrixOp::NORMAL, MatrixOp::NORMAL, 
                   m, n, k, alpha, rA.get(), lda, rB.get(), ldb, beta,
                   wRes.get(), ldc );
    }
}

/* ========================================================================= */

template<typename ValueType>
void GridVector<ValueType>::setDiagonal( const GridSection<ValueType>& diagonal, const int diagonalNumber )
{
    SCAI_ASSERT_EQ_ERROR( 2, this->nDims(), "setDiagonal only on two-dimensional grids." )

    const IndexType numRows = this->size( 0 );
    const IndexType numCols = this->size( 1 );

    IndexType offsetTarget = 0;
    IndexType sizesTarget[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesTarget[SCAI_GRID_MAX_DIMENSION];

    IndexType expectedDiagSize = numRows;

    if ( diagonalNumber >= 0 )
    {
        const IndexType col = diagonalNumber;  // column where the diagonal starts
        SCAI_ASSERT_VALID_INDEX( col, numCols, "diagonal " << diagonalNumber << " out of range" )
        expectedDiagSize = common::Math::min( numCols - col, numRows );
        GridSection<ValueType> colSection = (*this)( Range( 0, expectedDiagSize ), col );
        colSection.getDopeVector( offsetTarget, sizesTarget, distancesTarget );
    }
    else
    {
        const IndexType row = -diagonalNumber;  // row where the diagonal starts
        SCAI_ASSERT_VALID_INDEX( row, numRows, "diagonal " << diagonalNumber << " out of range" );
        expectedDiagSize = common::Math::min( numRows - row, numCols );
        GridSection<ValueType> colSection = (*this)( Range( row, expectedDiagSize ), 0 );
        colSection.getDopeVector( offsetTarget, sizesTarget, distancesTarget );
    }

    // We cannot define a section for the diagonal directly, but we can do it via a column section
    // where we fake later the distances

    IndexType offsetSource = 0;
    IndexType sizesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesSource[SCAI_GRID_MAX_DIMENSION];

    IndexType dimsSource = diagonal.getDopeVector( offsetSource, sizesSource, distancesSource );

    SCAI_ASSERT_EQ_ERROR( 1, dimsSource, "setDiagonal, diagonal is not one-dimensional" )
    SCAI_ASSERT_EQ_ERROR( expectedDiagSize, sizesSource[0], "size of diagonal does not match" )

    distancesTarget[0] = distancesTarget[0] + 1;  // this gives the diagonal

    static utilskernel::LAMAKernel<utilskernel::SectionKernelTrait::assign<ValueType> > assign;

    hmemo::ContextPtr loc = this->getContextPtr();

    assign.getSupportedContext( loc );

    {
        GridReadAccess<ValueType> rSource( diagonal.mGridVector, loc );
        GridWriteAccess<ValueType> wTarget( *this, loc );

        SCAI_CONTEXT_ACCESS( loc )

        const ValueType* sourcePtr = rSource.get() + offsetSource;
        ValueType* targetPtr = wTarget.get() + offsetTarget;
    
        common::BinaryOp op = common::BinaryOp::COPY;

        bool swap = false;

        assign[loc]( targetPtr, dimsSource, sizesSource, distancesTarget, sourcePtr, distancesSource, op, swap );
    }
}

/* -- IO ------------------------------------------------------------------- */

template<typename ValueType>
void GridVector<ValueType>::writeLocalToFile( FileIO& file ) const
{
    SCAI_LOG_INFO( logger, "write local part to file " << file << ": " << *this )

    file.writeGridArray( this->mLocalValues, this->localGrid() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void GridVector<ValueType>::readFromFile( FileIO& file )
{
    auto comm = file.getCommunicatorPtr();

    if ( file.getDistributedIOMode() == DistributedIOMode::SINGLE )
    {
        const PartitionId MASTER = 0;
        const PartitionId myRank = comm->getRank();

        IndexType nDims;
        IndexType globalGridDims[SCAI_GRID_MAX_DIMENSION];

        common::Grid grid;

        if ( myRank == MASTER )
        {
            SCAI_LOG_INFO( logger, *comm << ": read local array from file " << file )
            file.readGridArray( this->mLocalValues, grid );
            nDims = grid.nDims();
            grid.getSizes( globalGridDims );
        }

        // single distribution for grid distribution not available, so replicate it

        comm->bcast( &nDims, 1, MASTER );
        comm->bcast( globalGridDims, nDims, MASTER );
        comm->bcastArray( this->mLocalValues, MASTER );

        common::Grid globalGrid( nDims, globalGridDims );

        this->setDistributionPtr( dmemo::gridDistributionReplicated( globalGrid ) );
    }
    else
    {
        // INDEPENDENT or COLLECTIVE 

        common::Grid localGrid;

        file.readGridArray( this->mLocalValues, localGrid );

        auto dist = dmemo::gridDistributionByLocalGrid( localGrid, comm );

        this->setDistributionPtr( dist );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void GridVector<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "GridVector<" << this->getValueType() << ">" << "( size = " << size() << ", local = " << this->mLocalValues.size()
           << ", dist = " << this->getDistribution() << ", loc  = " << *this->getContextPtr() << " )";
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void GridVector<ValueType>:: operator=( const GridSection<ValueType>& other )
{
    GridSection<ValueType> fullSection = *this;
    fullSection = other;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
GridVector<ValueType>* GridVector<ValueType>::copy() const
{
    // create a new grid vector with the copy constructor
    return new GridVector<ValueType>( *this );
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( GridVector, SCAI_ARRAY_TYPES_HOST )

/* ========================================================================= */
/*       Free functions instantiatons                                        */
/* ========================================================================= */

#define GRID_VECTOR_FUNCTIONS( ValueType )                  \
                                                            \
    template GridVector<ValueType> gridVectorUndefined(     \
        const common::Grid& grid,                           \
        hmemo::ContextPtr);                                 \
                                                            \
    template GridVector<ValueType> gridVectorUndefined(     \
        dmemo::DistributionPtr dist,                        \
        hmemo::ContextPtr);                                 \

    SCAI_COMMON_LOOP( GRID_VECTOR_FUNCTIONS, SCAI_ARRAY_TYPES_HOST )

#undef GRID_VECTOR_FUNCTIONS

} /* end namespace lama */


} /* end namespace scai */
