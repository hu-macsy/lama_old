/**
 * @file GridSection.cpp
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
 * @brief Implementation of operations on distributed sections
 * @author Thomas Brandes
 * @date 10.05.2017
 */

#include <scai/lama/GridSection.hpp>
#include <scai/lama/GridVector.hpp>

#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>

#include <scai/utilskernel/SectionKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/loop.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
{

/* ---------------------------------------------------------------------------------------*/

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, GridSection<ValueType>::logger, "GridSection" )

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void GridSection<ValueType>::setDim( IndexType k, const Range& r )
{
   SecDimInfo& sec = mDimension[k];

   sec.globalRange[0] = r.mLB;
   sec.globalRange[1] = r.mUB;
   sec.globalRange[2] = r.mStride;

   sec.isLocalized = false;
   sec.isRange     = r.mStride != 0;

   if ( r.mUB == invalidIndex )
   {
       sec.globalRange[1] = mGlobalGrid.size( k );
   }

   if ( r.mStride != 0 )
   {
       sec.isRange = true;
       sec.globalRange[2] = r.mStride;
   }
   else
   {
       sec.isRange = false;
       sec.globalRange[2] = 1;
   }
}

/* ---------------------------------------------------------------------------------------*/

static IndexType correctUB( const IndexType lb, const IndexType ub, const IndexType str )
{
    if ( lb >= ub )
    {
        return ub;
    }
    else
    {
        return lb + ( ( ub - lb - 1 ) / str ) * str + 1;
    }
}

static int raiseLB ( const IndexType lb1, const IndexType lb2, const IndexType stride )
{ 
    if ( lb1 < lb2 )
    {
        return lb1 + ( ( lb2 - lb1 - 1 ) / stride + 1 ) * stride;
    }
    else
    {
        return lb1;
    }
}

/* ---------------------------------------------------------------------------------------*/

static void intersect( IndexType localRange[], const IndexType globalRange[], IndexType lb, IndexType size )
{
    // compute the local range of the global range for the owned values

    if ( globalRange[2] == 1 )
    {
        localRange[0] = common::Math::max( globalRange[0], lb );
        localRange[1] = common::Math::min( globalRange[1], lb + size );
        localRange[2] = 1;
    }
    else 
    {
        IndexType high = common::Math::min( correctUB( globalRange[0], globalRange[1], globalRange[2] ), lb + size );
        localRange[0] = raiseLB( globalRange[0], lb, globalRange[2] );
        localRange[1] = correctUB( localRange[0], high, globalRange[2] );
        localRange[2] = globalRange[2];
    }

    // move it back to origin

    // localRange[0] -= lb;
    // localRange[1] -= lb;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void GridSection<ValueType>::localize()
{
    const IndexType* localLB = mGridVector.localLB();

    // compute the local ranges of the global section

    for ( IndexType i = 0; i < mNDims; ++i )
    {
        SecDimInfo& sec = mDimension[i];

        if ( sec.isLocalized )
        {
            continue;
        }

        // localRange = intersect( globalRange, [lb,ub] )

        intersect( sec.localRange, sec.globalRange, localLB[i], mLocalGrid.size( i ) );

        if ( mGlobalGrid.size( i ) != mLocalGrid.size( i ) )
        {
            SCAI_LOG_INFO( logger, "dim " << i << ", (global) range = " << sec.globalRange[0] << ":" << sec.globalRange[1] 
                                   << ":" << sec.globalRange[2] << " is on " << localLB[i] << ":" << localLB[i] + mLocalGrid.size( i ) 
                                   << " (local) range = " << sec.localRange[0] << ":" << sec.localRange[1] << ":" << sec.localRange[2] )
        }
        else
        {
            SCAI_ASSERT_EQ_ERROR( sec.localRange[0], sec.globalRange[0], "mismatch" )
            // ub adapted SCAI_ASSERT_EQ_ERROR( sec.localRange[1], sec.globalRange[1], "mismatch" )
            SCAI_ASSERT_EQ_ERROR( sec.localRange[2], sec.globalRange[2], "mismatch" )
        }

        if ( sec.localRange[1] > sec.localRange[0] )
        {
            sec.localRange[0] -= localLB[i];
            sec.localRange[1] -= localLB[i];
        } 
        else
        { 
            sec.localRange[0] -= localLB[i];
            sec.localRange[1] = sec.localRange[0];
        }

        sec.isLocalized = true;
    }
}   

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>::GridSection( const GridVector<ValueType>& gridVector ) :

    mGridVector( const_cast<GridVector<ValueType>& >( gridVector ) ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )
{
    Range r;   // full range

    mNDims = mGlobalGrid.nDims();
   
    for ( IndexType idim = 0; idim < mNDims; idim++ )
    {
       setDim( idim, r );
    }
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>::GridSection( GridVector<ValueType>& gridVector, const Range& r0 ) :

    mGridVector( gridVector ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )
{
    mNDims = 1;

    setDim( 0, r0 );
}

template<typename ValueType>
GridSection<ValueType>::GridSection( GridVector<ValueType>& gridVector, const Range& r0, const Range& r1 ) :

    mGridVector( gridVector ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )
{
    mNDims = 2;

    setDim( 0, r0 );
    setDim( 1, r1 );
}

template<typename ValueType>
GridSection<ValueType>::GridSection( GridVector<ValueType>& gridVector, const Range& r0, const Range& r1, const Range& r2 ) :

    mGridVector( gridVector ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )
{
    mNDims = 3;

    setDim( 0, r0 );
    setDim( 1, r1 );
    setDim( 2, r2 );
}

template<typename ValueType>
GridSection<ValueType>::GridSection( GridVector<ValueType>& gridVector, const Range& r0, const Range& r1, const Range& r2, const Range& r3 ) :

    mGridVector( gridVector ),
    mGlobalGrid( gridVector.globalGrid() ),
    mLocalGrid( gridVector.localGrid() )
{
    mNDims = 4;

    setDim( 0, r0 );
    setDim( 1, r1 );
    setDim( 2, r2 );
    setDim( 3, r3 );
}

template<typename ValueType>
IndexType GridSection<ValueType>::getDopeVector( IndexType& offset, IndexType sizes[], IndexType distances[] ) const
{
    SCAI_LOG_DEBUG( logger, "get dope vector grid section, #dims = " << mNDims )

    GridSection<ValueType>& tmp = const_cast<GridSection<ValueType>&>( *this );
    tmp.localize();
 
    mLocalGrid.getDistances( distances );

    offset = 0;
    
    IndexType nSectionDims = 0;  // counts only ranges, not fixed elements

    for ( IndexType i = 0; i < mNDims; ++i )
    {
        const SecDimInfo& sec = mDimension[i];

        SCAI_LOG_DEBUG( logger, "secDim[" << i << "] : " << sec.globalRange[0] << ":" << sec.globalRange[1] << ":" << sec.globalRange[2] 
                                 << ", is range = " << sec.isRange )

        offset += sec.localRange[0] * distances[i];

        SCAI_LOG_DEBUG( logger, "offset = " << offset )

        if ( sec.isRange )
        {
            if ( sec.localRange[2] == 1 )
            {
                sizes[ nSectionDims ] = sec.localRange[1] - sec.localRange[0];
                distances[ nSectionDims] = distances[ i ];
            }
            else if ( sec.localRange[2] > 1 )
            {
                IndexType str = sec.localRange[2];
                sizes[ nSectionDims ] = ( sec.localRange[1] - sec.localRange[0] + str - 1 ) / str;
                distances[ nSectionDims] = distances[ i ] * str;
            }
            nSectionDims++;
        }
    }

    SCAI_LOG_DEBUG( logger, "section dims = " << nSectionDims << ", offset = " << offset );

    for ( IndexType k = 0; k < nSectionDims; ++k )
    {
        SCAI_LOG_DEBUG( logger, "size[" << k << "] = " << sizes[k] << ", distance = " << distances[k] )
    }

    return nSectionDims;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void GridSection<ValueType>::binOp ( const GridSection<ValueType>& other, common::BinaryOp op, bool swap )
{
    IndexType offsetSource = 0;
    IndexType offsetTarget = 0;

    IndexType sizesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType sizesTarget[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesTarget[SCAI_GRID_MAX_DIMENSION];

    IndexType dimsSource = other.getDopeVector( offsetSource, sizesSource, distancesSource );
    IndexType dimsTarget = getDopeVector( offsetTarget, sizesTarget, distancesTarget );

    SCAI_ASSERT_EQ_ERROR( dimsSource, dimsTarget, "section dimensions do not match" )

    for ( IndexType i = 0; i < dimsSource; ++i )
    {
        SCAI_ASSERT_EQ_ERROR( sizesSource[i], sizesTarget[i], "size mismatch for section dim = " << i )
    }

    static utilskernel::LAMAKernel<utilskernel::SectionKernelTrait::assign<ValueType> > assign;

    hmemo::ContextPtr loc = mGridVector.getContextPtr();

    assign.getSupportedContext( loc );

    {
        GridReadAccess<ValueType> rSource( other.mGridVector, loc );
        GridWriteAccess<ValueType> wTarget( mGridVector, loc );
    
        const ValueType* sourcePtr = rSource.get() + offsetSource;
        ValueType* targetPtr = wTarget.get() + offsetTarget;

        SCAI_CONTEXT_ACCESS( loc )

        assign[loc]( targetPtr, dimsSource, sizesSource, distancesTarget, sourcePtr, distancesSource, op, swap );
    }
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
template<typename OtherValueType>
void GridSection<ValueType>::unaryOp( const GridSection<OtherValueType>& other, common::UnaryOp op )
{
    IndexType offsetSource = 0;
    IndexType offsetTarget = 0;

    IndexType sizesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType sizesTarget[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesTarget[SCAI_GRID_MAX_DIMENSION];

    IndexType dimsSource = other.getDopeVector( offsetSource, sizesSource, distancesSource );
    IndexType dimsTarget = getDopeVector( offsetTarget, sizesTarget, distancesTarget );

    SCAI_ASSERT_EQ_ERROR( dimsSource, dimsTarget, "section dimensions do not match" )

    for ( IndexType i = 0; i < dimsSource; ++i )
    {
        SCAI_ASSERT_EQ_ERROR( sizesSource[i], sizesTarget[i], "size mismatch for section dim = " << i )
    }

    static utilskernel::LAMAKernel<utilskernel::SectionKernelTrait::unaryOp<ValueType, OtherValueType> > unaryOp;

    hmemo::ContextPtr loc = mGridVector.getContextPtr();

    unaryOp.getSupportedContext( loc );

    {
        GridReadAccess<OtherValueType> rSource( other.mGridVector, loc );
        GridWriteAccess<ValueType> wTarget( mGridVector, loc );

        const OtherValueType* sourcePtr = rSource.get() + offsetSource;
        ValueType* targetPtr = wTarget.get() + offsetTarget;

        SCAI_CONTEXT_ACCESS( loc )

        unaryOp[loc]( targetPtr, dimsSource, sizesSource, distancesTarget, sourcePtr, distancesSource, op );
    }
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType GridSection<ValueType>::get() const
{
    // ToDo: not yet for distributed grids

    SCAI_ASSERT_ERROR( mGridVector.getDistribution().isReplicated(), "get not yet supported for distributed vectors." )

    IndexType offset;
    IndexType sizes[SCAI_GRID_MAX_DIMENSION];
    IndexType distances[SCAI_GRID_MAX_DIMENSION];

    IndexType rank = getDopeVector( offset, sizes, distances );

    SCAI_ASSERT_EQ_ERROR( 0, rank, "section must have rank 0 to get the value" )

    const hmemo::HArray<ValueType>& vals = mGridVector.getLocalValues();

    return vals[ offset ];
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator= ( const GridSection<ValueType>& other )
{
    binOp( other, common::BinaryOp::COPY, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator+= ( const GridSection<ValueType>& other )
{
    binOp( other, common::BinaryOp::ADD, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator-= ( const GridSection<ValueType>& other )
{
    binOp( other, common::BinaryOp::SUB, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator*= ( const GridSection<ValueType>& other )
{
    binOp( other, common::BinaryOp::MULT, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator/= ( const GridSection<ValueType>& other )
{
    binOp( other, common::BinaryOp::DIVIDE, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void GridSection<ValueType>::binOpScalar ( const ValueType other, common::BinaryOp op, bool swap )
{
    IndexType sectionOffset = 0;
    IndexType sectionSizes[SCAI_GRID_MAX_DIMENSION];
    IndexType sectionDistances[SCAI_GRID_MAX_DIMENSION];

    IndexType sectionDims = getDopeVector( sectionOffset, sectionSizes, sectionDistances );

    SCAI_LOG_DEBUG( logger, "section dims = " << sectionDims << ", offset = " << sectionOffset );

    for ( IndexType k = 0; k < sectionDims; ++k )
    {
        SCAI_LOG_TRACE( logger, "size[" << k << "] = " << sectionSizes[k] << ", distance = " << sectionDistances[k] )
    }

    static utilskernel::LAMAKernel<utilskernel::SectionKernelTrait::assignScalar<ValueType> > assignScalar;
    hmemo::ContextPtr loc = mGridVector.getContextPtr();
    assignScalar.getSupportedContext( loc );

    {
        GridWriteAccess<ValueType> wGridVector( mGridVector, loc );
        ValueType* sectionPtr = wGridVector.get() + sectionOffset;

        SCAI_CONTEXT_ACCESS( loc )

        assignScalar[loc]( sectionPtr, sectionDims, sectionSizes, sectionDistances, other, op, swap );
    }
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator= ( const ValueType other )
{
    binOpScalar( other, common::BinaryOp::COPY, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator+= ( const ValueType other )
{
    binOpScalar( other, common::BinaryOp::ADD, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator-= ( const ValueType other )
{
    binOpScalar( other, common::BinaryOp::SUB, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator*= ( const ValueType other )
{
    binOpScalar( other, common::BinaryOp::MULT, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
GridSection<ValueType>& GridSection<ValueType>::operator/= ( const ValueType other )
{
    binOpScalar( other, common::BinaryOp::DIVIDE, false );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void GridSection<ValueType>::conj( const GridSection<ValueType>& other )
{
    unaryOp( other, common::UnaryOp::CONJ );
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void GridSection<ValueType>::assignTranspose( const GridSection<ValueType>& other, bool conjFlag )
{
    IndexType offsetSource = 0;
    IndexType offsetTarget = 0;

    IndexType sizesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType sizesTarget[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesSource[SCAI_GRID_MAX_DIMENSION];
    IndexType distancesTarget[SCAI_GRID_MAX_DIMENSION];

    IndexType dimsSource = other.getDopeVector( offsetSource, sizesSource, distancesSource );
    IndexType dimsTarget = getDopeVector( offsetTarget, sizesTarget, distancesTarget );

    SCAI_ASSERT_EQ_ERROR( 2, dimsSource, "transpose only for two-dimensional sections." )
    SCAI_ASSERT_EQ_ERROR( 2, dimsTarget, "transpose only for two-dimensional sections." )

    SCAI_ASSERT_EQ_ERROR( sizesSource[0], sizesTarget[1], "not matching for transpose" );
    SCAI_ASSERT_EQ_ERROR( sizesSource[1], sizesTarget[0], "not matching for transpose" );

    // section kernels work fine if we just switch distances of source section
 
    std::swap( distancesSource[0], distancesSource[1] );

    static utilskernel::LAMAKernel<utilskernel::SectionKernelTrait::unaryOp<ValueType, ValueType> > unaryOp;

    hmemo::ContextPtr loc = mGridVector.getContextPtr();

    unaryOp.getSupportedContext( loc );

    {
        GridReadAccess<ValueType> rSource( other.mGridVector, loc );
        GridWriteAccess<ValueType> wTarget( mGridVector, loc );

        const ValueType* sourcePtr = rSource.get() + offsetSource;
        ValueType* targetPtr = wTarget.get() + offsetTarget;

        common::UnaryOp op = conjFlag ? common::UnaryOp::CONJ : common::UnaryOp::COPY;

        SCAI_CONTEXT_ACCESS( loc )

        unaryOp[loc]( targetPtr, dimsSource, sizesSource, distancesTarget, sourcePtr, distancesSource, op );
    }
}

/* ---------------------------------------------------------------------------------------*/

// template GridSection<double>& GridSection<double>::operator= ( const GridSection<float>& other );

#define HARRAUTILS_SPECIFIER_LVL2( TargetType, SourceType )                          \
    template void GridSection<TargetType>::unaryOp(           \
            const GridSection<SourceType>& other, common::UnaryOp );           \

#define HARRAYUTILS_SPECIFIER( ValueType )                                      \
    SCAI_COMMON_LOOP_LVL2( ValueType, HARRAUTILS_SPECIFIER_LVL2, SCAI_ARRAY_TYPES_HOST )

SCAI_COMMON_LOOP( HARRAYUTILS_SPECIFIER, SCAI_ARRAY_TYPES_HOST )

#undef HARRAYUTILS_SPECIFIER
#undef HARRAYUTILS_SPECIFIER_LVL2

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( GridSection, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
