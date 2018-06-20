/**
 * @file StencilStorage.cpp
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
 * @brief Implementation and instantiation for template class StencilStorage.
 * @author Thomas Brandes
 * @date 04.06.2017
 */

// hpp
#include "StencilStorage.hpp"
#include <scai/sparsekernel/openmp/OpenMPStencilKernel.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/sparsekernel/StencilKernelTrait.hpp>

#include <scai/lama/storage/StorageMethods.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/Scalar.hpp>

// internal scai libraries
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/tasking/NoSyncToken.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>

namespace scai
{

using namespace hmemo;
using namespace dmemo;
using namespace utilskernel;
using namespace tasking;

using std::unique_ptr;
using std::shared_ptr;
using common::TypeTraits;
using common::BinaryOp;
using common::Grid;

namespace lama
{

template<typename ValueType>
StencilStorage<ValueType>::StencilStorage() :

    MatrixStorage<ValueType>( 0, 0, Context::getContextPtr() ),
    mGrid( common::Grid1D( 0 ) ),
    mStencil( common::Stencil1D<ValueType>() )
{
}

template<typename ValueType>
StencilStorage<ValueType>::StencilStorage( const common::Grid& grid, const common::Stencil<ValueType>&  stencil ) :

    MatrixStorage<ValueType>( grid.size(), grid.size(), Context::getContextPtr() ),
    mGrid( grid ),
    mStencil( stencil )
{
    // #dimension of grid and stencil must be equal

    SCAI_ASSERT_EQ_ERROR( grid.nDims(), stencil.nDims(), "dimensions of grid an stencil must be equal" )
}

template<typename ValueType>
StencilStorage<ValueType>::~StencilStorage()
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::conj()
{
    COMMON_THROWEXCEPTION( "conj unsupported" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> StencilStorage<ValueType>::l1Norm() const
{
    COMMON_THROWEXCEPTION( "l1Norm unsupported" )
    return 0;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> StencilStorage<ValueType>::l2Norm() const
{
    COMMON_THROWEXCEPTION( "l2Norm unsupported" )
    return 0;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> StencilStorage<ValueType>::maxNorm() const
{
    COMMON_THROWEXCEPTION( "maxNorm unsupported" )
    return 0;
}

template<typename ValueType>
void StencilStorage<ValueType>::setIdentity( const common::Grid& grid )
{
    mGrid = grid;

    switch ( grid.nDims() )
    {
        case 1 : 
        {
            common::Stencil1D<ValueType> stencil1D;
            stencil1D.addPoint( 0, 1 );
            mStencil = stencil1D;
            break;
        }
        case 2 : 
        {
            common::Stencil2D<ValueType> stencil2D;
            stencil2D.addPoint( 0, 0, 1 );
            mStencil = stencil2D;
            break;
        }
        case 3 : 
        {
            common::Stencil3D<ValueType> stencil3D;
            stencil3D.addPoint( 0, 0, 0, 1 );
            mStencil = stencil3D;
            break;
        }
        case 4 : 
        {
            common::Stencil4D<ValueType> stencil4D;
            stencil4D.addPoint( 0, 0, 0, 0, 1 );
            mStencil = stencil4D;
            break;
        }
        default:
            COMMON_THROWEXCEPTION( "unsupported stencil dim for grid " << grid )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t StencilStorage<ValueType>::getMemoryUsageImpl() const
{
    // just take the memory needed for the stencil data

    size_t nBytes1 = sizeof( IndexType ) * mStencil.nDims() * mStencil.nPoints();
    size_t nBytes2 = sizeof( ValueType ) * mStencil.nPoints();

    return nBytes1 + nBytes2;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
static const ValueType* getDiagonalPtr( const common::Stencil<ValueType>& stencil )
{
    const IndexType nDims = stencil.nDims();

    const int* stencilPositions = stencil.positions();
    const ValueType* stencilValues    = stencil.values();

    for ( IndexType k = 0; k < stencil.nPoints(); ++k )
    {
        bool isZero = true;

        for ( IndexType idim = 0; idim < nDims; ++idim )
        {
            isZero = isZero && stencilPositions[ k * nDims + idim ] == 0;
        }

        if ( isZero )
        {
            return stencilValues + k;
        }
    }

    return NULL;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::getDiagonal( HArray<ValueType>& array ) const
{
    ValueType diagonalValue = 0;   // stencil matrix has one single diagonal value

    const ValueType* diagonalPtr = getDiagonalPtr( mStencil );

    if ( diagonalPtr != NULL )
    {
        diagonalValue = *diagonalPtr;
    }

    SCAI_LOG_ERROR( logger, "diagonal value is = " << diagonalValue );

    HArrayUtils::setScalar( array, diagonalValue, common::BinaryOp::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::scale( const ValueType val )
{
    mStencil.scale( val );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::buildCSRSizes( HArray<IndexType>& sizeIA ) const
{
    ContextPtr ctx = Context::getHostPtr();    //  we build all data on the host

    IndexType n = this->getNumRows();

    std::unique_ptr<int[]> stencilOffsets( new int[ mStencil.nPoints() ] );

    IndexType gridDistances[SCAI_GRID_MAX_DIMENSION];

    mGrid.getDistances( gridDistances );

    mStencil.getLinearOffsets( stencilOffsets.get(), gridDistances );

    {
        WriteOnlyAccess<IndexType> sizes( sizeIA, n );
  
        sparsekernel::OpenMPStencilKernel::stencilLocalSizes( 
            sizes.get(), mStencil.nDims(), mGrid.sizes(), gridDistances, mGrid.borders(),
            mStencil.nPoints(), mStencil.positions() );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::buildCSRData( HArray<IndexType>& csrIA, HArray<IndexType>& csrJA, _HArray& csrValues ) const
{
    ContextPtr ctx = Context::getHostPtr();    //  we build all data on the host

    IndexType n = this->getNumRows();

    std::unique_ptr<int[]> stencilOffsets( new int[ mStencil.nPoints() ] );

    IndexType gridDistances[SCAI_GRID_MAX_DIMENSION];

    mGrid.getDistances( gridDistances );

    mStencil.getLinearOffsets( stencilOffsets.get(), gridDistances );

    for ( IndexType i = 0; i < mStencil.nPoints(); ++i )
    {
        SCAI_LOG_DEBUG( logger, "point = " << i << ", offset = " << stencilOffsets[i] << ", val = " << mStencil.values()[i] )
    }

    {
        WriteOnlyAccess<IndexType> sizes( csrIA, n );
  
        sparsekernel::OpenMPStencilKernel::stencilLocalSizes( 
            sizes.get(), mStencil.nDims(), mGrid.sizes(), gridDistances, mGrid.borders(),
            mStencil.nPoints(), mStencil.positions() );
    }

    // now build offset array

    IndexType nnz = HArrayUtils::scan1( csrIA );

    SCAI_LOG_DEBUG( logger, "Stencil -> CSR, nnz = " << nnz )

    // compute ja and values array

    HArray<ValueType> typedValues;
 
    {
        ReadAccess<IndexType> ia( csrIA );
        WriteOnlyAccess<IndexType> ja( csrJA, nnz );
        WriteOnlyAccess<ValueType> values( typedValues, nnz );
 
        sparsekernel::OpenMPStencilKernel::stencilLocalCSR( 
            ja.get(), values.get(), ia.get(),
            mStencil.nDims(), mGrid.sizes(), gridDistances, mGrid.borders(),
            mStencil.nPoints(), mStencil.positions(), mStencil.values(), stencilOffsets.get() );
    }

    if ( typedValues.getValueType() == csrValues.getValueType() )
    {
        typedValues.swap( static_cast<HArray<ValueType>&>( csrValues ) );
    }
    else
    {
        HArrayUtils::_assign( csrValues, typedValues );  // conversion required
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    IndexType curPos[SCAI_GRID_MAX_DIMENSION];   // grid position that corresponds to this row i
    IndexType newPos[SCAI_GRID_MAX_DIMENSION];   // for neighbored stencil points

    mGrid.gridPos( curPos, i );

    const IndexType nPoints = mStencil.nPoints();
    const IndexType nDims = mGrid.nDims();

    IndexType countValidPoints = 0;

    const int* stencilOffset = mStencil.positions();
    const ValueType* stencilValue = mStencil.values();

    {
        // we allocate the arrays with sufficient size but resize later

        WriteOnlyAccess<IndexType> wJA( jA, nPoints );
        WriteOnlyAccess<ValueType> wValues( values, nPoints );

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            // Build the new point
                
            for ( IndexType i = 0; i < nDims; ++i )
            {
                newPos[i] = curPos[i];
            }

            bool valid = mGrid.getOffsetPos( newPos, &stencilOffset[ p * nDims ] );

            if ( !valid ) 
            {
                continue;
            }

            IndexType col = mGrid.linearPos( newPos );
    
            // due to reflecting boundaries a col pos can appear twice

            bool found = false;
 
            for ( IndexType jj = 0; jj < countValidPoints; ++jj )
            {
                if ( wJA[jj] == col )
                {
                    found = true;
                    wValues[jj] += stencilValue[ p ];
                    break;
                }
            }

            if ( !found )
            {
                wJA[ countValidPoints ] = col;
                wValues[ countValidPoints ] = stencilValue[ p ];

                countValidPoints++;
            }
        }
    
        wJA.resize( countValidPoints );
        wValues.resize( countValidPoints );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::getRow( hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    values.resize( getNumColumns() );

    HArrayUtils::setScalar<ValueType>( values, 0, BinaryOp::COPY );

    HArray<IndexType> sparseIA;
    HArray<ValueType> sparseValues;

    getSparseRow( sparseIA, sparseValues, i );

    bool unique = true;  // sparseIA has only unique values

    HArrayUtils::scatter( values, sparseIA, unique, sparseValues, BinaryOp::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType StencilStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    HArray<IndexType> sparseIA;
    HArray<ValueType> sparseValues;

    getSparseRow( sparseIA, sparseValues, i );

    SCAI_ASSERT_EQ_DEBUG( sparseIA.size(), sparseValues.size(), "serious mismatch" )

    ValueType value = 0;

    {
        ReadAccess<IndexType> rIA( sparseIA );

        for ( IndexType jj = 0; jj < rIA.size(); ++jj )
        {
            if ( rIA[jj] == j )
            {
                value = sparseValues[jj];
                break;
            }
        }
    }

    return value;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* StencilStorage<ValueType>::incGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    LAMAKernel<sparsekernel::StencilKernelTrait::stencilGEMV<ValueType> > stencilGEMV;

    ContextPtr loc = this->getContextPtr();

    stencilGEMV.getSupportedContext( loc );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType gridDistances[SCAI_GRID_MAX_DIMENSION];

    IndexType width[2 * SCAI_GRID_MAX_DIMENSION];

    mGrid.getDistances( gridDistances );

    mStencil.getWidth( width );

    std::unique_ptr<int[]> stencilOffsets( new int[ mStencil.nPoints() ] );

    mStencil.getLinearOffsets( stencilOffsets.get(), gridDistances );

    SCAI_LOG_INFO ( logger, "incGEMV, grid = " << mGrid << ", stencil = " << mStencil << ", done at " << *loc )

    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc );

    stencilGEMV[loc]( wResult.get(), alpha, rX.get(), 
                      mGrid.nDims(), mGrid.sizes(), width, gridDistances, mGrid.borders(),
                      mStencil.nPoints(), mStencil.positions(), mStencil.values(),
                      stencilOffsets.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.CSR.jacobiIterate" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local stencil data, omega = " << omega )

    SCAI_ASSERT_EQ_DEBUG( getNumRows(), oldSolution.size(), "illegal size for old solution" )
    SCAI_ASSERT_EQ_DEBUG( getNumRows(), rhs.size(), "illegal size for rhs" )
    SCAI_ASSERT_EQ_DEBUG( getNumRows(), getNumColumns(), "jacobiIterate only on square matrices" )

    // solution = omega * ( rhs - B * oldSolution ) * dinv  + ( 1 - omega ) * oldSolution
    // solution = omega * rhs * dinv + ( 1 - omega ) * oldSolution;
    // solution -= B * oldSolution * ( omega * dinv )

    // okay that is tricky stuff, swap the diagonal entry with 0
    // would cause serious problems if stencil storage is used in other operation

    ValueType* diagonalPtr = const_cast<ValueType*>( getDiagonalPtr( mStencil ) );

    ValueType diagonalValue = *diagonalPtr;

    *diagonalPtr = ValueType( 0 );

    ValueType alpha = omega / diagonalValue;

    HArrayUtils::arrayPlusArray( solution, alpha, rhs, ValueType( 1 ) - omega, oldSolution, getContextPtr() );

    const bool async = false;

    incGEMV( solution, -alpha, oldSolution, async );

    *diagonalPtr = diagonalValue;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::matrixTimesVector(

    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    if ( common::isTranspose( op ) )
    {
        COMMON_THROWEXCEPTION( "transpose( A ) * x unsupported for stencil storage" )
    }

    if ( mGrid.size() == 0 )
    {
        return;
    }

    SCAI_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x 
                         << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.Stencil.matrixTimesVector" )

    if ( alpha == common::Constants::ZERO )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::compute( result, beta, BinaryOp::MULT, y, this->getContextPtr() );
        return;
    }

    ContextPtr loc = this->getContextPtr();

    // GEMV only implemented as y += A * x, so split

    // Step 1: result = beta * y

    if ( beta == common::Constants::ZERO )
    {
        result.clear();
        result.resize( getNumRows() );
        HArrayUtils::setScalar( result, ValueType( 0 ), BinaryOp::COPY, loc );
    }
    else
    {
        SCAI_ASSERT_EQUAL( y.size(), getNumRows(), "size mismatch y, beta = " << beta )
        HArrayUtils::compute( result, beta, BinaryOp::MULT, y, this->getContextPtr() );
    }

    bool async = false;

    SyncToken* token = incGEMV( result, alpha, x, async );

    SCAI_ASSERT( token == NULL, "syncrhonous execution cannot have token" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* StencilStorage<ValueType>::matrixTimesVectorAsync(

    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    if ( common::isTranspose( op ) )
    {
        COMMON_THROWEXCEPTION( "transpose( A ) * x unsupported for stencil storage" )
    }

    if ( mGrid.size() == 0 )
    {
        return NULL;
    }

    SCAI_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x 
                         << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.Stencil.matrixTimesVectorAsync" )

    if ( alpha == common::Constants::ZERO )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::compute( result, beta, BinaryOp::MULT, y, this->getContextPtr() );
        return NULL;
    }

    ContextPtr loc = this->getContextPtr();

    // GEMV only implemented as y += A * x, so split

    // Step 1: result = beta * y

    if ( beta == common::Constants::ZERO )
    {
        result.clear();
        result.resize( getNumRows() );
        HArrayUtils::setScalar( result, ValueType( 0 ), BinaryOp::COPY, loc );
    }
    else
    {
        SCAI_ASSERT_EQUAL( y.size(), getNumRows(), "size mismatch y, beta = " << beta )
        HArrayUtils::compute( result, beta, BinaryOp::MULT, y, this->getContextPtr() );
    }

    bool async = true;

    SyncToken* token = incGEMV( result, alpha, x, async );

    return token;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::assign( const _MatrixStorage& other )
{
    SCAI_ASSERT_EQ_ERROR( Format::STENCIL, other.getFormat(),
                          "stencil matrix/storage: assign only possible if other is also stencil" );

    SCAI_ASSERT_EQ_ERROR( this->getValueType(), other.getValueType(),
                          "stencil matrix/storage: assigon only supported here if other has same value type" );

    _MatrixStorage::setDimension( other.getNumRows(), other.getNumColumns() );

    const StencilStorage<ValueType> otherStencilStorage = static_cast<const StencilStorage<ValueType>&>( other );

    mGrid = otherStencilStorage.mGrid;

    mStencil = otherStencilStorage.mStencil;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::assignTranspose( const MatrixStorage<ValueType>& other )
{
    _MatrixStorage::setDimension( other.getNumColumns(), other.getNumRows() );

    SCAI_ASSERT_EQ_ERROR( Format::STENCIL, other.getFormat(),
                          "stencil matrix/storage: transposed argument must also be stencil matrix" );

    const StencilStorage<ValueType> otherStencilStorage = static_cast<const StencilStorage<ValueType>&>( other );

    mGrid = otherStencilStorage.mGrid;

    const common::Grid::BorderType* borders = mGrid.borders();

    for ( IndexType i = 0; i < mGrid.nDims(); ++i )
    {
        mGrid.setBorderType( i, borders[ 2 * i + 1 ], borders[ 2 * i ] );
    }

    mStencil.transpose( otherStencilStorage.mStencil );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::clear()
{
    mGrid = common::Grid1D( 0 );
    mStencil = common::Stencil1D<ValueType>( 1 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::purge() 
{
    clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
StencilStorage<ValueType>* StencilStorage<ValueType>::copy() const
{
    return new StencilStorage<ValueType>( *this );
}

/* ========================================================================= */
/*  Static fatory methods and related virtual methods                        */
/* ========================================================================= */

template<typename ValueType>
std::string StencilStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "StencilStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* StencilStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

template<typename ValueType>
const char* StencilStorage<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
MatrixStorageCreateKeyType StencilStorage<ValueType>::getCreateValue() const
{
    return MatrixStorageCreateKeyType( Format::STENCIL, common::getScalarType<ValueType>() );
}

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, StencilStorage<ValueType>::logger, "MatrixStorage.StencilStorage" )

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( StencilStorage, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
