/**
 * @file HArrayUtils.cpp
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
 * @brief Implemenation of the utility routines for HArrays using corresponding kernel routines.
 * @author Thomas Brandes
 * @date 10.10.2011
 */

// hpp
#include <scai/utilskernel/HArrayUtils.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/SparseKernelTrait.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>

#include <scai/utilskernel/mepr/UtilsWrapper.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/ScalarType.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/loop.hpp>

// std
#include <iostream>

namespace scai
{

using namespace hmemo;

using common::BinaryOp;
using common::CompareOp;
using common::UnaryOp;
using common::isBinarySupported;
using common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( HArrayUtils::logger, "HArrayUtils" )

/* --------------------------------------------------------------------------- */

void HArrayUtils::_assign( _HArray& target, const _HArray& source, const ContextPtr prefLoc )
{
    if ( &target == &source )
    {
        SCAI_LOG_INFO( logger, "assign skipped, target and source are aliased" )
    }
    else
    {
        mepr::UtilsWrapperTT<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::setArray( target, source, common::BinaryOp::COPY, prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::assign( HArray<TargetValueType>& target, const HArray<SourceValueType>& source, const ContextPtr prefLoc )
{
    if ( static_cast<const _HArray*>( &target ) == static_cast<const _HArray*>( &source ) )
    {
        SCAI_LOG_INFO( logger, "assign skipped, target and source are aliased" )
        return;
    }

    // set should be available on interface for each loc

    static LAMAKernel<UtilKernelTrait::set<TargetValueType, SourceValueType> > set;

    ContextPtr loc = prefLoc;

    if ( !loc )
    {
        // if no context is given we assign where source has a valid copy available
        loc = source.getValidContext();
    }

    set.getSupportedContext( loc );

    const IndexType n = source.size();

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<SourceValueType> sourceVals( source, loc );
    WriteOnlyAccess<TargetValueType> targetVals( target, loc, n );
    set[loc]( targetVals.get(), sourceVals.get(), n, common::BinaryOp::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::assignResized( 
    HArray<TargetValueType>& target, 
    const IndexType newSize,
    const HArray<SourceValueType>& source, 
    const ContextPtr prefLoc )
{
    if ( static_cast<const _HArray*>( &target ) == static_cast<const _HArray*>( &source ) )
    {
        // target and source are same, so no assign required, but maybe resize and fill up

        IndexType n = source.size();

        target.resize( newSize );

        if ( newSize > n )
        {
            // fill up with zero

            fillArraySection( target, n, 1, TargetValueType( 0 ), newSize - n, common::BinaryOp::COPY, prefLoc );
        }

        return;
    }

    static LAMAKernel<UtilKernelTrait::set<TargetValueType, SourceValueType> > set;
    static LAMAKernel<UtilKernelTrait::setVal<TargetValueType>> setVal;

    ContextPtr loc = prefLoc;

    if ( !loc )
    {
        // if no context is given we assign where source has a valid copy available
        loc = source.getValidContext();
    }

    set.getSupportedContext( loc, setVal );

    const IndexType n = std::min( source.size(), newSize );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<SourceValueType> sourceVals( source, loc );
    WriteOnlyAccess<TargetValueType> targetVals( target, loc, newSize );
    set[loc]( targetVals.get(), sourceVals.get(), n, common::BinaryOp::COPY );

    if ( n < newSize )
    {
        setVal[loc]( targetVals.get() + n, newSize - n, TargetValueType( 0 ), common::BinaryOp::COPY );
    }
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::_setArray(
    _HArray& target,
    const _HArray& source,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    mepr::UtilsWrapperTT<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::setArray( target, source, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::setArray(
    HArray<TargetValueType>& target,
    const HArray<SourceValueType>& source,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    SCAI_LOG_DEBUG( logger, "setArray<" << common::TypeTraits<TargetValueType>::id() << 
                    ", " << common::TypeTraits<SourceValueType>::id() << ">, target[ " << target.size()
                    << "] <- " << op << " source[" << source.size() << "]" )

    if ( op == common::BinaryOp::COPY )
    {
        assign( target, source, prefLoc );
        return;
    }

    SCAI_ASSERT_ERROR( common::isBinarySupported<TargetValueType>( op ),
                       op << " not supported for " << TypeTraits<TargetValueType>::id() )

    SCAI_ASSERT_EQ_ERROR( target.size(), source.size(), "array size mismatch for set( op = " << op << " )" )

    static LAMAKernel<UtilKernelTrait::set<TargetValueType, SourceValueType> > set;

    ContextPtr loc = prefLoc;

    if ( !loc )
    {
        // if no context is given we assign where target has a valid copy available
        loc = target.getValidContext();
    }

    set.getSupportedContext( loc );  // just in case, but set should be available on interface for each loc

    const IndexType n = source.size();

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<SourceValueType> sourceVals( source, loc );
    WriteAccess<TargetValueType> targetVals( target, loc );
    set[loc]( targetVals.get(), sourceVals.get(), n, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::appendArray(
    hmemo::HArray<ValueType>& array1,
    const hmemo::HArray<ValueType>& array2,
    hmemo::ContextPtr context )
{
    const IndexType n1 = array1.size();
    const IndexType n2 = array2.size();

    static LAMAKernel<UtilKernelTrait::set<ValueType, ValueType> > set;

    ContextPtr loc = context;

    if ( !loc )
    {
        loc = array1.getValidContext();
    }

    set.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteAccess<ValueType> wArray1( array1, loc );
    wArray1.resize( n1 + n2 );
    ReadAccess<ValueType> rArray2( array2, loc );

    set[loc]( wArray1.get() + n1, rArray2.get(), n2, common::BinaryOp::COPY );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::_setArraySection(
    _HArray& target,
    const IndexType targetOffset,
    const IndexType targetStride,
    const _HArray& source,
    const IndexType sourceOffset,
    const IndexType sourceStride,
    const IndexType n,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    mepr::UtilsWrapperTT<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::
    setArraySection(
        target, targetOffset, targetStride,
        source, sourceOffset, sourceStride,
        n, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::setArraySection(
    HArray<TargetValueType>& target,
    const IndexType targetOffset,
    const IndexType targetStride,
    const HArray<SourceValueType>& source,
    const IndexType sourceOffset,
    const IndexType sourceStride,
    const IndexType n,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<TargetValueType>( op ),
                       op << " not supported for " << TypeTraits<TargetValueType>::id() )

    // in contrary to setArray we assume correctly allocated target

    static LAMAKernel<UtilKernelTrait::setSection<TargetValueType, SourceValueType> > setSection;

    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        // optional argument has not been set, take default

        loc = source.getValidContext();
    }

    setSection.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<SourceValueType> sourceVals( source, loc );
    WriteAccess<TargetValueType> targetVals( target, loc );

    setSection[loc]( targetVals.get() + targetOffset, targetStride,
                     sourceVals.get() + sourceOffset, sourceStride,
                     n, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::fillArraySection(
    HArray<ValueType>& array,
    const IndexType offset,
    const IndexType stride,
    const ValueType val,
    const IndexType n,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::fillSection<ValueType> > fillSection;

    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        // optional argument has not been set, take default

        loc = array.getValidContext();  // avoids moving of data
    }

    fillSection.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteAccess<ValueType> wValues( array, loc );

    fillSection[loc]( wValues.get() + offset, stride, val, n, op );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::_gather(
    _HArray& target,
    const _HArray& source,
    const HArray<IndexType>& indexes,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    // use metaprogramming to call the gather version with the correct value types for target and source
    mepr::UtilsWrapperTT<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::gather( target, source, indexes, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename SourceValueType>
void HArrayUtils::_sparseGather(
    _HArray& target,
    const SourceValueType sourceZeroValue,
    const HArray<SourceValueType>& sourceNonZeroValues,
    const HArray<IndexType>& sourceNonZeroIndexes,
    const HArray<IndexType>& indexes,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    // use metaprogramming to call sparseGather version with the correct value type for target

    mepr::UtilsWrapperT<SourceValueType, SCAI_ARRAY_TYPES_HOST_LIST>::
           sparseGather( target, sourceZeroValue, sourceNonZeroValues, sourceNonZeroIndexes, indexes, op, prefLoc );

}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::sparseGather(
    HArray<TargetValueType>& target,
    const SourceValueType sourceZeroValue,
    const HArray<SourceValueType>& sourceNonZeroValues,
    const HArray<IndexType>& sourceNonZeroIndexes,
    const HArray<IndexType>& indexes,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "sparseGather<" << TypeTraits<TargetValueType>::id()
                            << ", " << TypeTraits<SourceValueType>::id() << ">" )

    const IndexType n = indexes.size();   // number of values to access

    // choose location for the operation where source array is currently valid

    static LAMAKernel<UtilKernelTrait::setGatherSparse<TargetValueType, SourceValueType> > setGatherSparse;

    ContextPtr loc = prefLoc;

    // default: assume that we need the target values where the indexes are

    if ( loc == ContextPtr() )
    {
        loc = indexes.getValidContext();
    }

    setGatherSparse.getSupportedContext( loc );

    if ( op == BinaryOp::COPY )
    {
        target.clear();       // old values are no more needed
        target.resize( n );    
    }
    else
    {
        SCAI_ASSERT_EQ_ERROR( target.size(), indexes.size(), "size mismatch" )
    }

    SCAI_ASSERT_EQ_ERROR( sourceNonZeroValues.size(), sourceNonZeroIndexes.size(), "serious size mismatch" )

    const IndexType nnz = sourceNonZeroValues.size();

    WriteAccess<TargetValueType> wTarget( target, loc );
    ReadAccess<SourceValueType> rSourceVals( sourceNonZeroValues, loc );
    ReadAccess<IndexType> rSourceIndexes( sourceNonZeroIndexes, loc );
    ReadAccess<IndexType> rIndexes( indexes, loc );

    setGatherSparse[loc] ( wTarget.get(), sourceZeroValue, rSourceVals.get(), rSourceIndexes.get(), nnz, rIndexes.get(), op, n );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::_scatter(
    _HArray& target,
    const HArray<IndexType>& indexes,
    const bool unique,
    const _HArray& source,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    // use metaprogramming to call the scatter version with the correct value types for target and source
    mepr::UtilsWrapperTT<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::scatter( target, indexes, unique, source, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::gather(
    HArray<TargetValueType>& target,
    const HArray<SourceValueType>& source,
    const HArray<IndexType>& indexes,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<TargetValueType>( op ),
                       op << " not supported for " << TypeTraits<TargetValueType>::id() )

    // check alias, only possible for same value type

    if ( target.getValueType() == source.getValueType() )
    {
        const HArray<TargetValueType>& tsource = reinterpret_cast<const HArray<TargetValueType>&>( source );
        SCAI_ASSERT_NE_ERROR( &target, &tsource, "gather: alias of source and target unsupported" )
    }

    SCAI_REGION( "HArray.gather" )

    SCAI_LOG_INFO( logger, "target[ " << target.size() << " ] " << op << " = source [ indexes [ "
                   << indexes.size() << " ] : " << source.size() << " ]" )

    // choose location for the operation where source array is currently valid
    static LAMAKernel<UtilKernelTrait::setGather<TargetValueType, SourceValueType> > setGather;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = source.getValidContext();
    }

    setGather.getSupportedContext( loc );
    const IndexType n = indexes.size();

    if ( op == BinaryOp::COPY )
    {
        //  target[i] = source[ indexes[i] ]
        //  alias of target and indexes is supported

        SCAI_CONTEXT_ACCESS( loc )
        ReadAccess<SourceValueType> rSource( source, loc );
        ReadAccess<IndexType> rIndexes( indexes, loc );
        WriteOnlyAccess<TargetValueType> wTarget( target, loc, n );

        setGather[loc] ( wTarget.get(), rSource.get(), rIndexes.get(), op, n );
    }
    else
    {
        //  target[i] op= source[ indexes[i] ]

        SCAI_ASSERT_GE_ERROR( target.size(), n, "target array too small" )

        WriteAccess<TargetValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc )
        ReadAccess<SourceValueType> rSource( source, loc );
        ReadAccess<IndexType> rIndexes( indexes, loc );
        setGather[loc] ( wTarget.get(), rSource.get(), rIndexes.get(), op, n );
    }
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::scatter(
    HArray<TargetValueType>& target,
    const HArray<IndexType>& indexes,
    const bool unique,
    const HArray<SourceValueType>& source,
    const BinaryOp op,
    const ContextPtr prefLoc )
{
    if ( unique )
    {
        // each supported binary operator can be used

        SCAI_ASSERT_ERROR( isBinarySupported<TargetValueType>( op ),
                           op << " not supported for " << TypeTraits<TargetValueType>::id() )
    }
    else
    {
        // only binary operators that allow atomic updates

        SCAI_ASSERT(  op == BinaryOp::COPY || op == BinaryOp::ADD || op == BinaryOp::SUB  ,
                      "Unsupported reduction op " << op << " for not unique indexes" )
    }

    SCAI_REGION( "HArray.scatter" )

    SCAI_LOG_INFO( logger, "target[ indexes[ " << indexes.size() << " ] : " << target.size()
                   << " ] = source[ " << source.size() << " ]" )

    if ( target.size() == 0 )
    {
        SCAI_ASSERT_EQ_ERROR( 0, indexes.size(), "scatter on empty array" )
        return;
    }

    SCAI_ASSERT_DEBUG( HArrayUtils::validIndexes( indexes, target.size(), prefLoc ),
                       "illegal scatter index, target has size " << target.size() )

    static LAMAKernel<UtilKernelTrait::setScatter<TargetValueType, SourceValueType> > setScatter;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = source.getValidContext();
    }

    setScatter.getSupportedContext( loc );
    const IndexType n = indexes.size();
    WriteAccess<TargetValueType> wTarget( target, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<SourceValueType> rSource( source, loc );
    ReadAccess<IndexType> rIndexes( indexes, loc );
    //  target[ indexes[i] ] = source[i]
    setScatter[loc] ( wTarget.get(), rIndexes.get(), unique, rSource.get(), op, n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::transpose(
    hmemo::HArray<ValueType>& target,
    const IndexType n1,
    const IndexType n2,
    const HArray<ValueType>& source,
    const bool conj,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( source.size(), n1 * n2, "source array has not " << n1 << " x " << n2 << " elements." )

    static LAMAKernel<blaskernel::BLASKernelTrait::geam<ValueType> > geam;

    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        // default location: where we have valid copy of the source data
        loc = source.getValidContext();
    }

    geam.getSupportedContext( loc );

    if ( &target == &source )
    {
        HArray<ValueType> tmpTarget;
        transpose( tmpTarget, n1, n2, source, conj, loc );
        target = std::move( tmpTarget );
    }
    else
    {
        ReadAccess<ValueType> rSource( source, loc );
        WriteOnlyAccess<ValueType> wTarget( target, loc, n1 * n2 );

        SCAI_CONTEXT_ACCESS( loc )

        common::MatrixOp op = conj ? common::MatrixOp::CONJ_TRANSPOSE : common::MatrixOp::TRANSPOSE;

        geam[loc]( wTarget.get(), n2, n1, n2,
                   ValueType( 1 ), rSource.get(), n1, op,
                   ValueType( 0 ), rSource.get(), n2, common::MatrixOp::NORMAL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setScalar(
    HArray<ValueType>& target,
    const ValueType value,
    const BinaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<ValueType>( op ),
                       op << " not supported for " << TypeTraits<ValueType>::id() )

    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        // default location: where we have valid copy of the data
        loc = target.getValidContext();
    }

    setVal.getSupportedContext( loc );
    const IndexType n = target.size();
    SCAI_LOG_INFO( logger, target << " = " << value << ", to do at " << *loc << ", n = " << n )

    if ( op == BinaryOp::COPY )
    {
        // Note: very important is to specify the size n here as it might not have been allocated
        WriteOnlyAccess<ValueType> wTarget( target, loc, n );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( wTarget.get(), n, value, op );
    }
    else
    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( wTarget.get(), n, value, op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setSameValue(
    HArray<ValueType>& array,
    const IndexType n,
    const ValueType value,
    ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        // default location: where we have valid copy of the data
        loc = array.getValidContext();
    }

    setVal.getSupportedContext( loc );

    WriteOnlyAccess<ValueType> wTarget( array, loc, n );
    SCAI_CONTEXT_ACCESS( loc )
    setVal[loc]( wTarget.get(), n, value, common::BinaryOp::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::_setSameValue(
    hmemo::_HArray& array,
    const IndexType n,
    const ValueType val,
    hmemo::ContextPtr prefLoc )
{
    // use meta-programming to resolve the type of target array

    mepr::UtilsWrapperT<ValueType, SCAI_ARRAY_TYPES_HOST_LIST>::setSameValue( array, n, val, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setVal(
    HArray<ValueType>& target,
    const IndexType index,
    const ValueType val,
    const BinaryOp op )
{
    SCAI_ASSERT_ERROR( isBinarySupported<ValueType>( op ),
                       op << " not supported for " << TypeTraits<ValueType>::id() )

    // setting single value will directly copy to the device with the valid incarnation

    ContextPtr loc = target.getValidContext();

    if ( op == BinaryOp::COPY )
    {
        WriteAccess<ValueType> wTarget( target, loc );
        wTarget.setValue( val, index );
    }
    else
    {
        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
        setVal.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        WriteAccess<ValueType> wTarget( target, loc );
        setVal[loc]( wTarget.get() + index, 1, val, op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::getVal(
    const HArray<ValueType>& array,
    const IndexType index )
{
    SCAI_ASSERT_DEBUG( index < array.size(), "index = " << index << " out of range for array = " << array );
    // get the data from a valid context, avoids any memory copy.
    ContextPtr loc = array.getValidContext();
    ReadAccess<ValueType> rArray( array, loc );
    ValueType val;
    rArray.getValue( val, index );
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::reduce(
    const HArray<ValueType>& array,
    const BinaryOp redOp,
    const ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;
    // preferred location: where valid values of the array are available
    ContextPtr loc = array.getValidContext( prefLoc );
    reduce.getSupportedContext( loc );
    ReadAccess<ValueType> readArray( array, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType zero = common::zeroBinary<ValueType>( redOp );
    ValueType redVal = reduce[loc]( readArray.get(), readArray.size(), zero, redOp );
    return redVal;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> HArrayUtils::l1Norm( const HArray<ValueType>& array, const ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        return ValueType( 0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;
    // preferred location: where valid values of the array are available
    ContextPtr loc = array.getValidContext( prefLoc );
    asum.getSupportedContext( loc );
    ReadAccess<ValueType> readArray( array, loc );
    SCAI_CONTEXT_ACCESS( loc )
    RealType<ValueType> result = asum[loc]( n, readArray.get(), 1 );
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> HArrayUtils::l2Norm( const HArray<ValueType>& array, const ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        return ValueType( 0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::nrm2<ValueType> > nrm2;
    // preferred location: where valid values of the array are available
    ContextPtr loc = array.getValidContext( prefLoc );
    nrm2.getSupportedContext( loc );
    ReadAccess<ValueType> readArray( array, loc );
    SCAI_CONTEXT_ACCESS( loc )
    RealType<ValueType> result = nrm2[loc]( n, readArray.get(), 1 );
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> HArrayUtils::maxDiffNorm(
    const HArray<ValueType>& array1,
    const HArray<ValueType>& array2,
    ContextPtr prefLoc )
{
    // is just a special case of a reduction with two combined arrays

    return reduce2( array1, array2, BinaryOp::SUB, BinaryOp::ABS_MAX, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::reduce2(
    const HArray<ValueType>& array1,
    const HArray<ValueType>& array2,
    const BinaryOp binOp,
    const BinaryOp redOp,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<ValueType>( binOp ),
                       binOp << " not supported for " << TypeTraits<ValueType>::id() )

    ValueType zero = common::zeroBinary<ValueType>( redOp ); // checks also for valid op

    SCAI_ASSERT_EQUAL( array1.size(), array2.size(), "array size mismatch for reduce2" )

    const IndexType n = array1.size();

    static LAMAKernel<UtilKernelTrait::reduce2<ValueType> > reduce2;

    ContextPtr loc = prefLoc;

    // Rule for default location: where array1 has valid values

    if ( loc == ContextPtr() )
    {
        loc = array1.getValidContext();
    }

    reduce2.getSupportedContext( loc );

    ReadAccess<ValueType> readArray1( array1, loc );
    ReadAccess<ValueType> readArray2( array2, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType redVal = reduce2[loc]( readArray1.get(),
                                     readArray2.get(),
                                     n,
                                     binOp,
                                     zero,
                                     redOp );
    return redVal;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool HArrayUtils::all(
    const HArray<ValueType>& array1,
    const CompareOp op,
    const HArray<ValueType>& array2,
    ContextPtr prefLoc )
{
    const IndexType n = array1.size();

    static LAMAKernel<UtilKernelTrait::allCompare<ValueType> > allCompare;

    ContextPtr loc = prefLoc;

    // Rule for default location: where array1 has valid values

    if ( loc == ContextPtr() )
    {
        loc = array1.getValidContext();
    }

    allCompare.getSupportedContext( loc );

    ReadAccess<ValueType> readArray1( array1, loc );
    ReadAccess<ValueType> readArray2( array2, loc );

    SCAI_CONTEXT_ACCESS( loc )

    bool allVal = allCompare[loc]( readArray1.get(),
                                   readArray2.get(),
                                   n,
                                   op );
    return allVal;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool HArrayUtils::allScalar(
    const HArray<ValueType>& array,
    const CompareOp op,
    const ValueType val,
    ContextPtr prefLoc )
{
    const IndexType n = array.size();

    static LAMAKernel<UtilKernelTrait::allCompareScalar<ValueType> > allCompare;

    ContextPtr loc = prefLoc;

    // Rule for default location: where array1 has valid values

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    allCompare.getSupportedContext( loc );

    ReadAccess<ValueType> readArray( array, loc );

    SCAI_CONTEXT_ACCESS( loc )

    bool allVal = allCompare[loc]( readArray.get(),
                                   val,
                                   n,
                                   op );
    return allVal;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::dotProduct(
    const HArray<ValueType>& array1,
    const HArray<ValueType>& array2,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQUAL( array1.size(), array2.size(), "array size mismatch for dotproduct" )
    static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;
    ContextPtr loc = prefLoc;

    // Rule for default location: where array1 has valid values

    if ( loc == ContextPtr() )
    {
        loc = array1.getValidContext();
    }

    dot.getSupportedContext( loc );
    ReadAccess<ValueType> readArray1( array1, loc );
    ReadAccess<ValueType> readArray2( array2, loc );
    SCAI_CONTEXT_ACCESS( loc )
    const ValueType res = dot[loc]( readArray1.size(), readArray1.get(), 1, readArray2.get(), 1 );
    return res;

    // Note: other solution: reduce2( array1, array2, n, BinaryOp::MULT, 0, BinaryOp::ADD )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::axpy(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "result += " << alpha << " * x"
                   ", x = " << x << ", result = " << result )

    if ( alpha == common::Constants::ZERO )
    {
        return;
    }

    SCAI_ASSERT_EQUAL( result.size(), x.size(), "size mismatch" )
    IndexType n = x.size();
    static LAMAKernel<blaskernel::BLASKernelTrait::axpy<ValueType> > axpy;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = result.getValidContext();
    }

    axpy.getSupportedContext( loc );  // axpy must be available at loc
    ReadAccess<ValueType> xAccess( x, loc );
    WriteAccess<ValueType> resultAccess( result, loc, true );
    SCAI_CONTEXT_ACCESS( loc )
    axpy[loc]( n, alpha, xAccess.get(), 1, resultAccess.get(), 1 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::arrayPlusArray(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "result = " << alpha << " * x + " << beta << " * y" <<
                   ", x = " << x << ", y = " << y << ", result = " << result )

    // check for zero terms as we do not need read access and assert correct sizes
 
    bool noSwapScalar = false;

    if ( beta == common::Constants::ZERO )
    {
        // result = alpha * x [ + 0 * y ]

        if ( alpha == common::Constants::ZERO )
        {
            // result = 0 * x + 0 * y, be careful about size of result

            IndexType n = common::Math::max( x.size(), y.size() );
            result.clear();
            result.resize( n );
            setScalar( result, ValueType( 0 ), BinaryOp::COPY, prefLoc );
        }
        else
        {
            binaryOpScalar( result, x, alpha, BinaryOp::MULT, noSwapScalar, prefLoc );
        }

        return;
    }

    if ( alpha == common::Constants::ZERO )
    {
        // result = beta * y
        binaryOpScalar( result, y, beta, BinaryOp::MULT, noSwapScalar, prefLoc );
        return;
    }

    if ( &y == &result && beta == common::Constants::ONE )
    {
        axpy( result, alpha, x, prefLoc );
        return;
    }

    if ( &x == &result && alpha == common::Constants::ONE )
    {
        axpy( result, beta, y, prefLoc );
        return;
    }

    SCAI_ASSERT_EQUAL( x.size(), y.size(), "size mismatch" )
    const IndexType n = x.size();
    static LAMAKernel<blaskernel::BLASKernelTrait::sum<ValueType> > sum;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    sum.getSupportedContext( loc );
    // no alias checks are done here
    ReadAccess<ValueType> xAccess( x, loc );
    ReadAccess<ValueType> yAccess( y, loc );
    WriteOnlyAccess<ValueType> resultAccess( result, loc, x.size() );
    SCAI_CONTEXT_ACCESS( loc )
    sum[loc]( n, alpha, xAccess.get(), beta, yAccess.get(), resultAccess.get() );
}

template<typename ValueType>
void HArrayUtils::arrayPlusScalar(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    ContextPtr prefLoc )
{
    // operation is result = alpha * x + b (elementwise)
    SCAI_LOG_INFO( logger, "result = " << alpha << " * x + " << beta << ", x = " << x
                   << ", result = " << result )

    // check for zero terms as we do not need read access and assert correct sizes

    if ( alpha == common::Constants::ZERO ) // result = b
    {
        setSameValue( result, x.size(), beta, prefLoc );
        return;
    }

    if ( &x == &result && alpha == common::Constants::ONE ) // result += b (elementwise)
    {
        setScalar( result, beta, BinaryOp::ADD, prefLoc );
        return;
    }

    if ( alpha == common::Constants::ONE ) // result = x + b (elementwise)
    {
        binaryOpScalar( result, x, beta, BinaryOp::ADD, false, prefLoc );
        return;
    }

    if ( beta == common::Constants::ZERO )
    {
        if ( alpha == common::Constants::ZERO ) // result = 0
        {
            setSameValue<ValueType>( result, x.size(), 0, prefLoc );
        }
        else // result = alpha * x
        {
            binaryOpScalar( result, x, alpha, BinaryOp::MULT, false, prefLoc );
        }

        return;
    }

    // else full operation: result = alpha * x + b (elementwise)

    const IndexType n = x.size();

    static LAMAKernel<UtilKernelTrait::scaleVectorAddScalar<ValueType> > scaleVectorAddScalar;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    scaleVectorAddScalar.getSupportedContext( loc );

    // due to possible alias of result and x, write access must follow read

    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    SCAI_CONTEXT_ACCESS( loc )

    scaleVectorAddScalar[loc]( wResult.get(), rX.get(), n, alpha, beta );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::arrayTimesArray(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const HArray<ValueType>& y,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQUAL( x.size(), y.size(), "size mismatch" )
    const IndexType n = x.size();

    if ( alpha == common::Constants::ZERO )
    {
        // result = alpha * x * y -> result = 0

        setSameValue<ValueType>( result, n, 0, prefLoc );
        return;
    }

    static LAMAKernel<UtilKernelTrait::binaryOp<ValueType> > binaryOp;

    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    binaryOp.getSupportedContext( loc );

    // no alias checks are done here
    {
        ReadAccess<ValueType> xAccess( x, loc );
        ReadAccess<ValueType> yAccess( y, loc );
        WriteOnlyAccess<ValueType> resultAccess( result, loc, x.size() );
        SCAI_CONTEXT_ACCESS( loc )
        binaryOp[loc]( resultAccess.get(), xAccess.get(), yAccess.get(), n, BinaryOp::MULT );
    }

    if ( alpha != common::Constants::ONE )
    {
        setScalar( result, alpha, BinaryOp::MULT, prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::unaryOp(
    HArray<ValueType>& result,
    const HArray<ValueType>& x,
    const UnaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT( common::isUnarySupported<ValueType>( op ),
                 op << " not supported for " << TypeTraits<ValueType>::id() )

    const IndexType n = x.size();

    static LAMAKernel<UtilKernelTrait::unaryOp<ValueType> > unaryOpKernel;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    unaryOpKernel.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // due to possible alias of result and x, write access must follow read

    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    unaryOpKernel[loc]( wResult.get(), rX.get(), n, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::binaryOp(
    HArray<ValueType>& result,
    const HArray<ValueType>& x,
    const HArray<ValueType>& y,
    const BinaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT( isBinarySupported<ValueType>( op ),
                 op << " not supported for " << TypeTraits<ValueType>::id() )

    const IndexType n = x.size();

    SCAI_ASSERT_GE_ERROR( y.size(), n, "binaryOp: y has to small size" )

    static LAMAKernel<UtilKernelTrait::binaryOp<ValueType> > binaryOp;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    binaryOp.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // due to possible alias of result and x, write access must follow read

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    binaryOp[loc]( wResult.get(), rX.get(), rY.get(), n, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::binaryOpScalar(
    HArray<ValueType>& result,
    const HArray<ValueType>& x,
    const ValueType y,
    const BinaryOp op,
    const bool swapScalar,
    ContextPtr prefLoc )
{
    SCAI_ASSERT( isBinarySupported<ValueType>( op ),
                 op << " not supported for " << TypeTraits<ValueType>::id() )

    const IndexType n = x.size();

    static LAMAKernel<UtilKernelTrait::binaryOpScalar<ValueType> > binaryOpScalar;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    binaryOpScalar.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // due to possible alias of result and x, write access must follow read

    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    binaryOpScalar[loc]( wResult.get(), rX.get(), y, n, op, swapScalar );
}

/* --------------------------------------------------------------------------- */

bool HArrayUtils::validIndexes(
    const HArray<IndexType>& array,
    const IndexType size,
    const ContextPtr prefLoc )
{
    const IndexType n = array.size();
    static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;
    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    validIndexes.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> rValues( array, loc );
    bool valid = validIndexes[loc]( rValues.get(), n, size );
    return valid;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool HArrayUtils::isSorted(
    const HArray<ValueType>& array,
    const CompareOp op,
    ContextPtr prefLoc )
{
    const IndexType n = array.size();
    static LAMAKernel<UtilKernelTrait::isSorted<ValueType> > isSorted;
    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    isSorted.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rValues( array, loc );
    bool sorted = isSorted[loc]( rValues.get(), n, op );
    return sorted;
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::setOrder( HArray<IndexType>& array, IndexType n, ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();   // should be first context
    }

    setOrder.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteOnlyAccess<IndexType> wArray( array, loc, n );
    setOrder[loc]( wArray.get(), n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setSequence( HArray<ValueType>& array, ValueType startValue, ValueType inc, IndexType n, ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::setSequence<ValueType> > setSequence;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();   // should be first context
    }

    setSequence.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteOnlyAccess<ValueType> wArray( array, loc, n );
    setSequence[loc]( wArray.get(), startValue, inc, n );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::setRandom( _HArray& array,
                             const IndexType bound,
                             const ContextPtr prefLoc )
{
    // use meta-programming to call fillRandom<ValueType> with the type of array

    mepr::UtilsWrapper< SCAI_ARRAY_TYPES_HOST_LIST>::fillRandom( array, bound, 1.0f, prefLoc );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::setSparseRandom( _HArray& array,
                                   const float fillRate,
                                   const IndexType bound,
                                   const ContextPtr prefLoc )
{
    // use meta-programming to call fillRandom<ValueType> with the type of array

    mepr::UtilsWrapper< SCAI_ARRAY_TYPES_HOST_LIST>::fillRandom( array, bound, fillRate, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::fillRandom( HArray<ValueType>& array,
                              const IndexType bound,
                              const float fillRate,
                              const ContextPtr prefLoc )
{
    ContextPtr loc = Context::getHostPtr();   // currently only available on host

    const IndexType n = array.size();

    if ( fillRate >= 1.0f )
    {
        WriteOnlyAccess<ValueType> wArray( array, loc, n );

        for ( IndexType i = 0; i < n; ++i )
        {
            wArray[i] = common::Math::random<ValueType>( bound );
        }
    }
    else
    {
        if ( n > 0 && !array.isInitialized() )
        {
            COMMON_THROWEXCEPTION( "filling of uninitialized array: " << array )
        }

        WriteAccess<ValueType> wArray( array, loc );

        for ( IndexType i = 0; i < n; ++i )
        {
            bool fillIt = common::Math::randomBool( fillRate );
 
            if ( fillIt )
            {
                wArray[i] = common::Math::random<ValueType>( bound );
            }
        }
    }

    if ( prefLoc != ContextPtr() )
    {
        array.prefetch( prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::randomSparseIndexes( HArray<IndexType>& array, const IndexType n, const float probability )
{
    IndexType cnt = 0;

    WriteOnlyAccess<IndexType> wArray( array, n );

    for ( IndexType i = 0; i < n; ++i )
    {
        bool takeIt = common::Math::randomBool( probability );

        if ( takeIt )
        {
            wArray[ cnt++ ] = i;
        }
    }

    wArray.resize( cnt );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::scan1(
    HArray<ValueType>& array,
    ContextPtr prefLoc )
{
    const IndexType n = array.size();

    static LAMAKernel<UtilKernelTrait::scan<ValueType> > scan;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    scan.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteAccess<ValueType> wValues( array, loc );

    wValues.resize( n + 1 );   // One additional element will be added at the end

    bool exclusive = true;
    bool append = true;
    ValueType zero = 0;

    ValueType total = scan[loc]( wValues.get(), n, zero, exclusive, append );

    return total;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::scan(
    HArray<ValueType>& array,
    const ValueType first,
    const bool exclusive,
    ContextPtr prefLoc )
{
    const IndexType n = array.size();

    static LAMAKernel<UtilKernelTrait::scan<ValueType> > scan;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    scan.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteAccess<ValueType> wValues( array, loc );
 
    bool append = false;   // do not change the size

    ValueType total = scan[loc]( wValues.get(), n, first, exclusive, append );

    return total;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::unscan(
    HArray<ValueType>& array,
    ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n < 1 )
    {
        return ValueType( 0 );
    }

    static LAMAKernel<UtilKernelTrait::unscan<ValueType> > unscan;
    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    unscan.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wValues( array, loc );
    ValueType first = unscan[loc]( wValues.get(), n - 1 );

    // One additional element will be removed from end
    wValues.resize( n - 1 );

    return first;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::sort(
    HArray<IndexType>* perm,
    HArray<ValueType>* outValues,
    const HArray<ValueType>& inValues,
    const bool ascending,
    ContextPtr prefLoc )
{
    const IndexType n = inValues.size();

    if ( n == 0 )
    {
        if ( perm )
        {
            perm->clear();
        }

        return;
    }

    static LAMAKernel<UtilKernelTrait::sort<ValueType> > sort;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = inValues.getValidContext();
    }

    sort.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // ReadAccess for inValues is always before WriteAccess for outValues, so alias is okay

    ReadAccess<ValueType> rValues( inValues, loc );

    if ( perm == NULL )
    {
        WriteOnlyAccess<ValueType> wValues( *outValues, loc, n );
        sort[loc]( NULL, wValues.get(), rValues.get(), n, ascending );
    }
    else if ( outValues == NULL )
    {
        WriteOnlyAccess<IndexType> wPerm( *perm, loc, n );
        sort[loc]( wPerm.get(), NULL, rValues.get(), n, ascending );
    }
    else
    {
        WriteOnlyAccess<ValueType> wValues( *outValues, loc, n );
        WriteOnlyAccess<IndexType> wPerm( *perm, loc, n );
        sort[loc]( wPerm.get(), wValues.get(), rValues.get(), n, ascending );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::sortSparseEntries(
    HArray<IndexType>& indexes,
    HArray<ValueType>& values,
    const bool ascending,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( indexes.size(), values.size(), "keys and values must have same number of entries" )

    const IndexType n = indexes.size();

    if ( n < 2 )
    {
        return;
    }

    static LAMAKernel<UtilKernelTrait::sortInPlace<ValueType> > sortInPlace;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = indexes.getValidContext();
    }

    sortInPlace.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteAccess<IndexType> wIndexes( indexes, loc );
    WriteAccess<ValueType> wValues( values, loc );
    sortInPlace[loc]( wIndexes.get(), wValues.get(), n, ascending );
}

/* --------------------------------------------------------------------------- */

template<typename BucketType>
void HArrayUtils::bucketSort(
    HArray<IndexType>& offsets,
    HArray<IndexType>& perm,
    const HArray<BucketType>& array,
    const BucketType nb,
    ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        perm.clear();
        IndexType zero = 0;
        offsets.setSameValue( nb + 1, zero );  // offsets = { 0, 0, ..., 0 }
        return;
    }

    static LAMAKernel<UtilKernelTrait::countBuckets<BucketType> > countBuckets;
    static LAMAKernel<UtilKernelTrait::scan<IndexType> > scan;
    static LAMAKernel<UtilKernelTrait::sortInBuckets<BucketType> > sortInBuckets;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    countBuckets.getSupportedContext( loc, scan, sortInBuckets );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<BucketType> bucketMap( array, loc );
    WriteOnlyAccess<IndexType> sizes( offsets, loc, static_cast<IndexType>( nb + 1 ) );

    countBuckets[loc]( sizes.get(), nb, bucketMap, n );
    IndexType first = 0;
    IndexType total = scan[loc]( sizes.get(), nb, first, true, true );

    // Note: total can be < n if array contains values < 0 or >= nb

    WriteOnlyAccess<IndexType> wPerm( perm, loc, total );
    sortInBuckets[loc]( wPerm, sizes, nb, bucketMap, n );
}

/* --------------------------------------------------------------------------- */

template<typename BucketType>
void HArrayUtils::bucketCount(
    HArray<IndexType>& bucketSizes,
    const HArray<BucketType>& array,
    const BucketType nb,
    ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        IndexType zeroVal = 0;
        bucketSizes.setSameValue( nb, zeroVal );
        return;
    }

    static LAMAKernel<UtilKernelTrait::countBuckets<BucketType> > countBuckets;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    countBuckets.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<BucketType> bucketMap( array, loc );

    bucketSizes.reserve( loc, nb + 1 );

    WriteOnlyAccess<IndexType> wSizes( bucketSizes, loc, nb );

    countBuckets[loc]( wSizes, nb, bucketMap, n );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::inversePerm(
    HArray<IndexType>& invPerm,
    const HArray<IndexType>& perm,
    ContextPtr prefLoc )
{
    IndexType n = perm.size();

    if ( n == 0 )
    {
        invPerm.clear();  // reset it to zero, done
        return;
    }

    SCAI_ASSERT_NE_ERROR( &invPerm, &perm, "inversePerm: alias not supported" )

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = perm.getValidContext();
    }

    static LAMAKernel<UtilKernelTrait::setInversePerm> setInversePerm;
    static LAMAKernel<UtilKernelTrait::setVal<IndexType> > setVal;
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;

    setInversePerm.getSupportedContext( loc, setVal, reduce );

    // it is also checked that perm is really a permutation, otherwise invPerm would contain strange values

    {
        WriteOnlyAccess<IndexType> wInversePerm( invPerm, loc, n );
        ReadAccess<IndexType> rPerm( perm, loc );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( wInversePerm.get(), n, n, BinaryOp::COPY );
        setInversePerm[loc]( wInversePerm.get(), rPerm.get(), n );
        IndexType maxIndex = reduce[loc]( wInversePerm.get(), n, 0, BinaryOp::MAX );
        SCAI_ASSERT_ERROR( maxIndex < n, "Perm array does not cover all row indexes, #n = " << n );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::mergeSort(
    HArray<ValueType>& values,
    const HArray<IndexType>& offsets,
    bool ascending,
    ContextPtr prefLoc )
{
    mergeSortOptional( values, NULL, offsets, ascending, prefLoc );
}

template<typename ValueType>
void HArrayUtils::mergeSort(
    HArray<ValueType>& values,
    HArray<IndexType>& perm,
    const HArray<IndexType>& offsets,
    bool ascending,
    ContextPtr prefLoc )
{
    mergeSortOptional( values, &perm, offsets, ascending, prefLoc );
}

template<typename ValueType>
void HArrayUtils::mergeSortOptional(
    HArray<ValueType>& values,
    HArray<IndexType>* perm,
    const HArray<IndexType>& offsets,
    bool ascending,
    ContextPtr )
{
    IndexType n  = values.size();       // number of values to sort
    IndexType nb = offsets.size() - 1;  // number of sorted subarray

    SCAI_LOG_INFO( logger, "mergeSort of values[" << n << "] with " << nb << " sorted subarrays." )

    if ( perm )
    {
        SCAI_ASSERT_EQ_ERROR( values.size(), perm->size(), "array size mismatch" )
    }

    // we need a (writeable) copy of the offset array for running offsets

    HArray<IndexType> copyOffsets( offsets );

    // Temporary array required for the sorted values, in-place not supported

    HArray<ValueType> sorted;
    HArray<IndexType> newPerm;

    {
        WriteAccess<IndexType> wOffsets( copyOffsets );
        ReadAccess<IndexType> rOffsets( offsets );
        WriteOnlyAccess<ValueType> wSorted( sorted, n );
        ReadAccess<ValueType> rValues( values );

        std::unique_ptr<WriteOnlyAccess<IndexType> > wPerm;
        std::unique_ptr<ReadAccess<IndexType> > rPerm;

        if ( perm )
        {
            wPerm.reset( new WriteOnlyAccess<IndexType>( newPerm, n ) );
            rPerm.reset( new ReadAccess<IndexType>( *perm ) );
        }

        for ( IndexType i = 0; i < n; ++i )
        {
            // find the next minimal element

            IndexType nextIndex = invalidIndex;

            for ( IndexType k = 0; k < nb; ++k )
            {
                if ( wOffsets[k] == rOffsets[k + 1] )
                {
                    continue;  // no more values in subarray k
                }

                if ( nextIndex == invalidIndex )
                {
                    nextIndex = k;
                }
                else if ( ascending && rValues[wOffsets[k]] < rValues[wOffsets[nextIndex]] )
                {
                    nextIndex = k;
                }
                else if ( !ascending && rValues[wOffsets[k]] > rValues[wOffsets[nextIndex]] )
                {
                    nextIndex = k;
                }
            }

            SCAI_ASSERT_NE_ERROR( nextIndex, invalidIndex, "no more elements found" )

            IndexType& pos = wOffsets[nextIndex];

            SCAI_LOG_DEBUG( logger, "Next elem in subarray = " << nextIndex
                            << ", pos = " << pos << ", val = " << rValues[pos] );

            // sort in the next minimal element

            wSorted[i] = rValues[pos];

            if ( perm )
            {
                ( *wPerm )[i] = ( *rPerm )[pos];
            }

            pos++;
        }

        // Proof: wOffsets[k] == rOffsets[k+1]

        for ( IndexType k = 0; k < nb; ++k )
        {
            SCAI_ASSERT_EQ_ERROR( rOffsets[k + 1], wOffsets[k], "serious problem during mergesort" )
        }
    }

    // set the sorted values back in the input array, also for perm if available

    values.swap( sorted );

    if ( perm )
    {
        perm->swap( newPerm );
    }
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
void HArrayUtils::mergeSortOptional(
    HArray<ComplexFloat>& ,
    HArray<IndexType>*,
    const HArray<IndexType>&,
    bool,
    ContextPtr )
{
    COMMON_THROWEXCEPTION( "unsupported" )
}

template<>
void HArrayUtils::mergeSortOptional(
    HArray<ComplexDouble>& ,
    HArray<IndexType>*,
    const HArray<IndexType>&,
    bool,
    ContextPtr )
{
    COMMON_THROWEXCEPTION( "unsupported" )
}

template<>
void HArrayUtils::mergeSortOptional(
    HArray<ComplexLongDouble>& ,
    HArray<IndexType>*,
    const HArray<IndexType>&,
    bool,
    ContextPtr )
{
    COMMON_THROWEXCEPTION( "unsupported" )
}

#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::elimDoubles(
    HArray<IndexType>& indexes,
    HArray<ValueType>& values,
    const common::BinaryOp op )
{
    SCAI_ASSERT_EQ_ERROR( indexes.size(), values.size(), "serious size mismatch" )

    IndexType n = indexes.size();

    if ( n <= 1 )
    {
        return;
    }

    // mergeSparse elims also double values on sorted lists, just set one list empty

    static LAMAKernel<SparseKernelTrait::mergeSparse<ValueType> > mergeSparse;

    // default location for merge

    ContextPtr loc = indexes.getValidContext();

    mergeSparse.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    {
        WriteAccess<IndexType> wIndexes( indexes, loc );
        WriteAccess<ValueType> wValues( values, loc );

        IndexType nMerged = mergeSparse[loc]( wIndexes.get(), wValues.get(),
                                              wIndexes.get(), wValues.get(), n,
                                              NULL, NULL, 0, op );

        // nMerged can be smaller than n if one index list contained multiple values

        wIndexes.resize( nMerged );
        wValues.resize( nMerged );
    }
}

/* --------------------------------------------------------------------------- */

template<typename TargetType, typename SourceType>
void HArrayUtils::buildSparseArray(
    HArray<TargetType>& sparseArray,
    HArray<IndexType>& sparseIndexes,
    const HArray<SourceType>& denseArray,
    const SourceType zeroValue,
    ContextPtr prefLoc )
{
    const IndexType n = denseArray.size();
    static LAMAKernel<SparseKernelTrait::countNonZeros<SourceType> > countNonZeros;
    static LAMAKernel<SparseKernelTrait::compress<TargetType, SourceType> > compress;
    ContextPtr loc = prefLoc;

    // default location for conversion: where we have the dense values

    if ( loc == ContextPtr() )
    {
        loc = denseArray.getValidContext();
    }

    compress.getSupportedContext( loc, countNonZeros );
    SCAI_CONTEXT_ACCESS( loc )
    SourceType eps = TypeTraits<SourceType>::eps1();
    ReadAccess<SourceType> rDenseArray( denseArray, loc );
    // we count the non-zeros at first to have sizes for sparse data
    IndexType sparseN = countNonZeros[loc]( rDenseArray.get(), n, zeroValue, eps );
    SCAI_LOG_DEBUG( logger, "counted nonZeros for zero = " << zeroValue << ", dense " << n << " is now " << sparseN )
    WriteOnlyAccess<TargetType> wSparseArray( sparseArray, loc, sparseN );
    WriteOnlyAccess<IndexType> wSparseIndexes( sparseIndexes, loc, sparseN );
    IndexType sparseN1 = compress[loc]( wSparseArray.get(), wSparseIndexes.get(), rDenseArray.get(), n, zeroValue, eps );
    SCAI_ASSERT_EQ_ERROR( sparseN, sparseN1, "zero = " << zeroValue << ", counted nonzeros does not match with compessed size" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType HArrayUtils::insertSorted(
    HArray<ValueType>& array,
    const ValueType value, 
    ContextPtr )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    if ( common::TypeTraits<ValueType>::stype != common::TypeTraits<RealType>::stype )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
        return 0;
    }
    else
    {
        RealType value1 = value;

        IndexType n = array.size();
 
        HArray<RealType>& array1 = reinterpret_cast<HArray<RealType>&>( array );

        WriteAccess<RealType> wArray( array1 );

        wArray.resize( n + 1 );

        IndexType pos = n;

        // move up entries of index array until we have a smaller element

        for ( IndexType i = n;  i-- > 0;  )
        {
            if ( wArray[i] > value1 ) 
            {
                pos = i;
                wArray[i+1] = wArray[i];
            }
            else
            {
                break;
            }
        }

        wArray[pos] = value1;

        return pos;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::insertAtPos(
    HArray<ValueType>& array,
    const IndexType pos,
    const ValueType value,
    ContextPtr )
{
    IndexType n = array.size();

    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, n + 1, "illegal insert position" )

    WriteAccess<ValueType> wArray( array );

    wArray.resize( n + 1 );

    // set the value at pos and swap the remaining values up like a bubble

    ValueType swapValue = value;

    for ( IndexType i = pos; i <= n; ++i )
    {
        std::swap( swapValue, wArray[i] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::buildDenseArray(
    HArray<ValueType>& denseArray,
    const IndexType denseN,
    const HArray<ValueType>& sparseArray,
    const HArray<IndexType>& sparseIndexes,
    const ValueType zero,
    ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "buildDenseArray<" << common::TypeTraits<ValueType>::id() << ">, N = " << denseN 
                             << ", zero = " << zero << ", nnz = " << sparseArray.size() )

    const IndexType sparseN = sparseArray.size();

    SCAI_ASSERT_LE_ERROR( sparseN, denseN, "sparse array cannot have more entries than dense array." )

    SCAI_ASSERT_DEBUG( HArrayUtils::validIndexes( sparseIndexes, denseN, prefLoc ),
                       "Illegal index for sparse indexes, dense array size = " << denseN );

    setSameValue( denseArray, denseN, zero, prefLoc );

    bool unique = true;
    HArrayUtils::scatter( denseArray, sparseIndexes, unique, sparseArray, BinaryOp::COPY, prefLoc );
}

/* --------------------------------------------------------------------------- */

IndexType HArrayUtils::findPosInSortedIndexes( const HArray<IndexType>& indexes, const IndexType index )
{
    ReadAccess<IndexType> rIndexes( indexes );

    IndexType first = 0;
    IndexType last  = indexes.size();

    while ( first < last )
    {
        IndexType middle = first + ( last - first ) / 2;

        if ( rIndexes[middle] == index )
        {
            return middle;
        }
        else if ( rIndexes[middle] > index )
        {
            last = middle;
        }
        else
        {
            first = middle + 1;
        }
    }

    return invalidIndex;
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::findPosInSortedIndexesV( HArray<IndexType>& outPos,
                                           const HArray<IndexType>& indexes, 
                                           const HArray<IndexType> inPos )
{
    const IndexType n = inPos.size();
    const IndexType m = indexes.size();

    // Note: alias of outPos and inPos is safe !!

    ReadAccess<IndexType> rIn( inPos );
    ReadAccess<IndexType> rIndexes( indexes );
    WriteOnlyAccess<IndexType> wOut( outPos, n );
    OpenMPUtils::binarySearch( wOut.get(), rIndexes.get(), m, rIn.get(), n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::addSparse(
    HArray<IndexType>& resultIndexes,
    HArray<ValueType>& resultValues,
    const HArray<IndexType>& indexes1,
    const HArray<ValueType>& values1,
    const ValueType zero1,
    const ValueType alpha,
    const HArray<IndexType>& indexes2,
    const HArray<ValueType>& values2,
    const ValueType zero2,
    const ValueType beta,
    ContextPtr prefLoc )
{
    static LAMAKernel<SparseKernelTrait::countAddSparse> countAddSparse;
    static LAMAKernel<SparseKernelTrait::addSparse<ValueType> > addSparse;

    ContextPtr loc = prefLoc;

    // default location for conversion: where we have the dense values

    if ( loc == ContextPtr() )
    {
        loc = indexes1.getValidContext();
    }

    addSparse.getSupportedContext( loc, countAddSparse );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIndexes1( indexes1, loc );
    ReadAccess<IndexType> rIndexes2( indexes2, loc );

    IndexType n1 = indexes1.size();
    IndexType n2 = indexes2.size();

    IndexType n = countAddSparse[loc]( rIndexes1.get(), n1, rIndexes2.get(), n2 );

    if ( n == n1 && n == n2 )
    {
        SCAI_LOG_DEBUG( logger, "addSparse, same pattern for both operands" )

        rIndexes1.release();
        rIndexes2.release();

        if ( ( &resultIndexes != &indexes1 ) && ( &resultIndexes != &indexes2 ) )
        {
            // no alias, so we have to copy the indexes

            resultIndexes = indexes1;
        }

        arrayPlusArray( resultValues, alpha, values1, beta, values2, loc );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "addSparse, different pattern, n = " << n << ", n1 = " << n1 << ", n2 = " << n2 )

        WriteOnlyAccess<IndexType> wIndexes( resultIndexes, loc, n );
        WriteOnlyAccess<ValueType> wValues( resultValues, loc, n );

        ReadAccess<ValueType> rValues1( values1, loc );
        ReadAccess<ValueType> rValues2( values2, loc );

        addSparse[loc]( wIndexes.get(), wValues.get(),
                        rIndexes1.get(), rValues1.get(), zero1, n1, alpha,
                        rIndexes2.get(), rValues2.get(), zero2, n2, beta );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::binaryOpSparse(
    HArray<IndexType>& resultIndexes,
    HArray<ValueType>& resultValues,
    const HArray<IndexType>& indexes1,
    const HArray<ValueType>& values1,
    const ValueType zero1,
    const HArray<IndexType>& indexes2,
    const HArray<ValueType>& values2,
    const ValueType zero2,
    const BinaryOp op,
    ContextPtr prefLoc )
{
    static LAMAKernel<SparseKernelTrait::countAddSparse> countAddSparse;
    static LAMAKernel<SparseKernelTrait::binopSparse<ValueType> > binopSparse;

    ContextPtr loc = prefLoc;

    // default location for conversion: where we have the dense values

    if ( loc == ContextPtr() )
    {
        loc = indexes1.getValidContext();
    }

    binopSparse.getSupportedContext( loc, countAddSparse );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIndexes1( indexes1, loc );
    ReadAccess<IndexType> rIndexes2( indexes2, loc );

    IndexType n1 = indexes1.size();
    IndexType n2 = indexes2.size();

    IndexType n = countAddSparse[loc]( rIndexes1.get(), n1, rIndexes2.get(), n2 );

    if ( n == n1 && n == n2 )
    {
        // sparse input arrays have same pattern, so we can call binary op for values

        rIndexes1.release();
        rIndexes2.release();

         
        resultIndexes = indexes1;
        binaryOp( resultValues, values1, values2, op, prefLoc );
    }
    else
    {
        WriteOnlyAccess<IndexType> wIndexes( resultIndexes, loc, n );
        WriteOnlyAccess<ValueType> wValues( resultValues, loc, n );

        ReadAccess<ValueType> rValues1( values1, loc );
        ReadAccess<ValueType> rValues2( values2, loc );

        binopSparse[loc]( wIndexes.get(), wValues.get(),
                          rIndexes1.get(), rValues1.get(), zero1, n1, 
                          rIndexes2.get(), rValues2.get(), zero2, n2, op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType HArrayUtils::allSparse(
    bool& allFlag,
    const HArray<IndexType>& indexes1,
    const HArray<ValueType>& values1,
    const ValueType zero1,
    const HArray<IndexType>& indexes2,
    const HArray<ValueType>& values2,
    const ValueType zero2,
    const CompareOp op,
    ContextPtr prefLoc )
{
    static LAMAKernel<SparseKernelTrait::allCompareSparse<ValueType> > allCompareSparse;

    ContextPtr loc = prefLoc;

    // default location for conversion: where we have the dense values

    if ( loc == ContextPtr() )
    {
        loc = indexes1.getValidContext();
    }

    allCompareSparse.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIndexes1( indexes1, loc );
    ReadAccess<IndexType> rIndexes2( indexes2, loc );

    IndexType n1 = indexes1.size();
    IndexType n2 = indexes2.size();

    ReadAccess<ValueType> rValues1( values1, loc );
    ReadAccess<ValueType> rValues2( values2, loc );

    return allCompareSparse[loc]( allFlag,
                                  rIndexes1.get(), rValues1.get(), zero1, n1, 
                                  rIndexes2.get(), rValues2.get(), zero2, n2, op );

}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::mergeSparse(
    HArray<IndexType>& resultIndexes,
    HArray<ValueType>& resultValues,
    const HArray<IndexType>& indexes1,
    const HArray<ValueType>& values1,
    const HArray<IndexType>& indexes2,
    const HArray<ValueType>& values2,
    const BinaryOp op,
    ContextPtr prefLoc )
{
    static LAMAKernel<SparseKernelTrait::countAddSparse> countAddSparse;
    static LAMAKernel<SparseKernelTrait::mergeSparse<ValueType> > mergeSparse;

    ContextPtr loc = prefLoc;

    // default location for conversion: where we have the dense values

    if ( loc == ContextPtr() )
    {
        loc = indexes1.getValidContext();
    }

    mergeSparse.getSupportedContext( loc, countAddSparse );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIndexes1( indexes1, loc );
    ReadAccess<IndexType> rIndexes2( indexes2, loc );

    IndexType n1 = indexes1.size();
    IndexType n2 = indexes2.size();

    IndexType n = countAddSparse[loc]( rIndexes1.get(), n1, rIndexes2.get(), n2 );

    SCAI_LOG_INFO( logger, "counted new sparse n = " << n << " from " << n1 << " + " << n2 )

    {
        WriteOnlyAccess<IndexType> wIndexes( resultIndexes, loc, n );
        WriteOnlyAccess<ValueType> wValues( resultValues, loc, n );

        ReadAccess<ValueType> rValues1( values1, loc );
        ReadAccess<ValueType> rValues2( values2, loc );

        IndexType nMerged = mergeSparse[loc]( wIndexes.get(), wValues.get(),
                                              rIndexes1.get(), rValues1.get(), n1, 
                                              rIndexes2.get(), rValues2.get(), n2, op );

        // nMerged can be smaller than n if one index list contained multiple values

        resultIndexes.resize( nMerged );
        resultValues.resize( nMerged );
    }
}

template <typename ValueType>
void HArrayUtils::mergeAndMap(
    hmemo::HArray<ValueType> & result,
    hmemo::HArray<IndexType> & xMap,
    hmemo::HArray<IndexType> & yMap,
    const hmemo::HArray<ValueType> & x,
    const hmemo::HArray<ValueType> & y,
    const common::CompareOp comparator,
    hmemo::ContextPtr prefLoc )
{
    // Currently we do all operations on the host, so ignore prefLoc
    (void) prefLoc;

    const auto stype = common::TypeTraits<ValueType>::stype;

    SCAI_ASSERT_ERROR(!common::isComplex(stype),
                      "merge is not applicable to complex numbers. ValueType == " << stype);
    SCAI_ASSERT_DEBUG(HArrayUtils::isSorted(x, comparator), "x must be sorted according to the specified comparator " << comparator);
    SCAI_ASSERT_DEBUG(HArrayUtils::isSorted(y, comparator), "y must be sorted according to the specified comparator " << comparator);

    using hmemo::ReadAccess;
    using hmemo::WriteAccess;
    using hmemo::WriteOnlyAccess;

    ReadAccess<ValueType> rx(x);
    ReadAccess<ValueType> ry(y);
    WriteOnlyAccess<IndexType> wxMap(xMap, x.size());
    WriteOnlyAccess<IndexType> wyMap(yMap, y.size());

    auto & z = result;
    const auto zSize = rx.size() + ry.size();
    WriteOnlyAccess<ValueType> wz(z, zSize);

    IndexType xIndex = 0;
    IndexType yIndex = 0;
    IndexType zIndex = 0;

    while ( xIndex < rx.size() && yIndex < ry.size() )
    {
        SCAI_ASSERT_DEBUG(zIndex < z.size(), "zIndex should never go out of bounds");

        const auto & xElement = rx[xIndex];
        const auto & yElement = ry[yIndex];

        if (common::compare(xElement, comparator, yElement))
        {
            wxMap[xIndex] = zIndex;
            wz[zIndex] = xElement;
            ++xIndex;
        }
        else
        {
            wyMap[yIndex] = zIndex;
            wz[zIndex] = yElement;
            ++yIndex;
        }
        ++zIndex;
    }

    while ( xIndex < rx.size() )
    {
        wxMap[xIndex] = zIndex;
        wz[zIndex++] = rx[xIndex++];
    }

    while ( yIndex < ry.size() )
    {
        wyMap[yIndex] = zIndex;
        wz[zIndex++] = ry[yIndex++];
    }

    SCAI_ASSERT(xIndex == x.size(), "must have read all elements in x.");
    SCAI_ASSERT(yIndex == y.size(), "must have read all elements in y.");
    SCAI_ASSERT(zIndex == z.size(), "must have written to all elements in z.");
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::selectComplexPart(
    HArray<RealType<ValueType> >& realValues,
    const HArray<ValueType>& complexValues,
    const common::ComplexPart kind,
    ContextPtr )
{
    // ToDo: define and use kernel routine for any device

    IndexType n = complexValues.size();

    ReadAccess<ValueType> rComplex( complexValues );
    WriteOnlyAccess<RealType<ValueType> > wReal( realValues, n );

    if ( kind == common::ComplexPart::REAL )
    {
        for ( IndexType i = 0; i < n; ++i )
        {
            wReal[i] = common::Math::real( rComplex[i] );
        }
    }
    else
    {
        for ( IndexType i = 0; i < n; ++i )
        {
            wReal[i] = common::Math::imag( rComplex[i] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::buildComplex(
    HArray<ValueType>& complexValues,
    const HArray<RealType<ValueType> >& realValues,
    const HArray<RealType<ValueType> >& imagValues,
    ContextPtr )
{
    // ToDo: define and use kernel routine for any context

    SCAI_ASSERT_EQ_ERROR( realValues.size(), imagValues.size(), "serious size mismatch" )

    IndexType n = realValues.size();

    ReadAccess<RealType<ValueType> > rReal( realValues );
    ReadAccess<RealType<ValueType> > rImag( imagValues );
    WriteOnlyAccess<ValueType> wComplex( complexValues, n );

    ValueType i = common::TypeTraits<ValueType>::imaginaryUnit();

    for ( IndexType k = 0; k < n; ++k )
    {
        wComplex[k] = rReal[k] + i * rImag[k];
    }
}

/* --------------------------------------------------------------------------- */

/** Makro for the instantiation of routines with two template arguments for source and target type. */

#define HARRAUTILS_SPECIFIER_LVL2( TargetType, SourceType )                          \
    template void HArrayUtils::gather<TargetType, SourceType>(                       \
            hmemo::HArray<TargetType>&,                                              \
            const hmemo::HArray<SourceType>&,                                        \
            const hmemo::HArray<IndexType>&,                                         \
            const BinaryOp,                                                          \
            const hmemo::ContextPtr );                                               \
    template void HArrayUtils::sparseGather<TargetType, SourceType>(                 \
            hmemo::HArray<TargetType>&,                                              \
            const SourceType,                                                        \
            const hmemo::HArray<SourceType>&,                                        \
            const hmemo::HArray<IndexType>&,                                         \
            const hmemo::HArray<IndexType>&,                                         \
            const BinaryOp,                                                          \
            const hmemo::ContextPtr );                                               \
    template void HArrayUtils::assign<TargetType, SourceType>(                       \
            hmemo::HArray<TargetType>&,                                              \
            const hmemo::HArray<SourceType>&,                                        \
            hmemo::ContextPtr );                                                     \
    template void HArrayUtils::assignResized<TargetType, SourceType>(                \
            hmemo::HArray<TargetType>&,                                              \
            const IndexType newSize,                                                 \
            const hmemo::HArray<SourceType>&,                                        \
            hmemo::ContextPtr );                                                     \
    template void HArrayUtils::setArray<TargetType, SourceType>(                     \
            hmemo::HArray<TargetType>&,                                              \
            const hmemo::HArray<SourceType>&,                                        \
            const BinaryOp,                                                          \
            hmemo::ContextPtr );                                                     \
    template void HArrayUtils::setArraySection<TargetType, SourceType>(              \
            hmemo::HArray<TargetType>&, const IndexType, const IndexType,            \
            const hmemo::HArray<SourceType>&, const IndexType, const IndexType,      \
            const IndexType,                                                         \
            const BinaryOp,                                                          \
            hmemo::ContextPtr );                                                     \
    template void HArrayUtils::scatter<TargetType, SourceType>(                      \
            hmemo::HArray<TargetType>&,                                              \
            const hmemo::HArray<IndexType>&,                                         \
            const bool,                                                              \
            const hmemo::HArray<SourceType>&,                                        \
            const BinaryOp,                                                          \
            const hmemo::ContextPtr );                                               \
    template void HArrayUtils::buildSparseArray<TargetType, SourceType>(             \
            hmemo::HArray<TargetType>&,                                              \
            hmemo::HArray<IndexType>&,                                               \
            const hmemo::HArray<SourceType>&,                                        \
            const SourceType,                                                        \
            hmemo::ContextPtr );

/** Makro for the instantiation of routines with one template argument for the value type. */

#define HARRAYUTILS_SPECIFIER( ValueType )                              \
    template void HArrayUtils::_sparseGather<ValueType>(                \
            hmemo::_HArray&,                                            \
            const ValueType,                                            \
            const hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<IndexType>&,                            \
            const BinaryOp,                                             \
            const hmemo::ContextPtr );                                  \
    template void HArrayUtils::setVal<ValueType>(                       \
            hmemo::HArray<ValueType>&,                                  \
            const IndexType,                                            \
            const ValueType,                                            \
            const BinaryOp );                                           \
    template ValueType HArrayUtils::getVal<ValueType>(                  \
            const hmemo::HArray<ValueType>&,                            \
            const IndexType );                                          \
    template void HArrayUtils::transpose<ValueType>(                    \
            hmemo::HArray<ValueType>&,                                  \
            const IndexType,                                            \
            const IndexType,                                            \
            const hmemo::HArray<ValueType>&,                            \
            const bool,                                                 \
            hmemo::ContextPtr);                                         \
    template void HArrayUtils::setScalar<ValueType>(                    \
            hmemo::HArray<ValueType>&,                                  \
            const ValueType,                                            \
            const BinaryOp,                                             \
            hmemo::ContextPtr);                                         \
    template void HArrayUtils::setSameValue<ValueType>(                 \
            hmemo::HArray<ValueType>&,                                  \
            const IndexType,                                            \
            const ValueType,                                            \
            hmemo::ContextPtr);                                         \
    template void HArrayUtils::_setSameValue<ValueType>(                \
            hmemo::_HArray&,                                            \
            const IndexType,                                            \
            const ValueType,                                            \
            hmemo::ContextPtr);                                         \
    template void HArrayUtils::fillArraySection<ValueType>(             \
            HArray<ValueType>& array,                                   \
            const IndexType,                                            \
            const IndexType,                                            \
            const ValueType,                                            \
            const IndexType,                                            \
            const BinaryOp,                                             \
            ContextPtr prefLoc );                                       \
    template void HArrayUtils::appendArray<ValueType>(                  \
            hmemo::HArray<ValueType>&,                                  \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template ValueType HArrayUtils::reduce<ValueType>(                  \
            const hmemo::HArray<ValueType>&,                            \
            const BinaryOp,                                             \
            hmemo::ContextPtr );                                        \
    template ValueType HArrayUtils::reduce2<ValueType>(                 \
            const hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const BinaryOp,                                             \
            const BinaryOp,                                             \
            hmemo::ContextPtr );                                        \
    template bool HArrayUtils::all<ValueType>(                          \
            const hmemo::HArray<ValueType>&,                            \
            const CompareOp,                                            \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template bool HArrayUtils::allScalar<ValueType>(                    \
            const hmemo::HArray<ValueType>&,                            \
            const CompareOp,                                            \
            const ValueType,                                            \
            hmemo::ContextPtr );                                        \
    template RealType<ValueType> HArrayUtils::maxDiffNorm<ValueType>(   \
            const hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::axpy<ValueType>(                         \
            hmemo::HArray<ValueType>&,                                  \
            const ValueType,                                            \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::arrayPlusArray<ValueType>(               \
            hmemo::HArray<ValueType>&,                                  \
            const ValueType,                                            \
            const hmemo::HArray<ValueType>&,                            \
            const ValueType,                                            \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::arrayPlusScalar<ValueType>(              \
            hmemo::HArray<ValueType>&,                                  \
            const ValueType,                                            \
            const hmemo::HArray<ValueType>&,                            \
            const ValueType, hmemo::ContextPtr );                       \
    template void HArrayUtils::arrayTimesArray<ValueType>(              \
            hmemo::HArray<ValueType>&,                                  \
            const ValueType,                                            \
            const hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template ValueType HArrayUtils::dotProduct<ValueType>(              \
            const hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template RealType<ValueType> HArrayUtils::l1Norm<ValueType>(        \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template RealType<ValueType> HArrayUtils::l2Norm<ValueType>(        \
            const hmemo::HArray<ValueType>&,                            \
            hmemo::ContextPtr );                                        \
    template bool HArrayUtils::isSorted<ValueType>(                     \
            const hmemo::HArray<ValueType>&,                            \
            const CompareOp,                                            \
            hmemo::ContextPtr );                                        \
    template ValueType HArrayUtils::scan1<ValueType>(                   \
            hmemo::HArray<ValueType>&,                                  \
            hmemo::ContextPtr );                                        \
    template ValueType HArrayUtils::scan<ValueType>(                    \
            hmemo::HArray<ValueType>&,                                  \
            const ValueType first,                                      \
            const bool exclusive,                                       \
            hmemo::ContextPtr );                                        \
    template ValueType HArrayUtils::unscan<ValueType>(                  \
            hmemo::HArray<ValueType>&, hmemo::ContextPtr );             \
    template void HArrayUtils::sort<ValueType>(                         \
            hmemo::HArray<IndexType>*,                                  \
            hmemo::HArray<ValueType>*,                                  \
            const hmemo::HArray<ValueType>&,                            \
            const bool,                                                 \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::sortSparseEntries<ValueType>(            \
            hmemo::HArray<IndexType>&,                                  \
            hmemo::HArray<ValueType>&,                                  \
            const bool,                                                 \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::mergeSort<ValueType>(                    \
            hmemo::HArray<ValueType>&,                                  \
            const hmemo::HArray<IndexType>&,                            \
            const bool,                                                 \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::mergeSort<ValueType>(                    \
            hmemo::HArray<ValueType>&,                                  \
            hmemo::HArray<IndexType>&,                                  \
            const hmemo::HArray<IndexType>&,                            \
            const bool,                                                 \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::elimDoubles<ValueType>(                  \
            hmemo::HArray<IndexType>&,                                  \
            hmemo::HArray<ValueType>&,                                  \
            const BinaryOp );                                           \
    template void HArrayUtils::setSequence<ValueType>(                  \
            hmemo::HArray<ValueType>&,                                  \
            ValueType,                                                  \
            ValueType,                                                  \
            IndexType,                                                  \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::fillRandom<ValueType>(                   \
            hmemo::HArray<ValueType>&,                                  \
            IndexType,                                                  \
            float,                                                      \
            hmemo::ContextPtr );                                        \
    template IndexType HArrayUtils::insertSorted(                       \
            hmemo::HArray<ValueType>& array,                            \
            const ValueType value,                                      \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::insertAtPos(                             \
            hmemo::HArray<ValueType>& array,                            \
            const IndexType pos,                                        \
            const ValueType val,                                        \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::buildDenseArray(                         \
            hmemo::HArray<ValueType>&,                                  \
            const IndexType,                                            \
            const hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<IndexType>&,                            \
            const ValueType,                                            \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::addSparse(                               \
            hmemo::HArray<IndexType>&,                                  \
            hmemo::HArray<ValueType>&,                                  \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const ValueType,                                            \
            const ValueType,                                            \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const ValueType,                                            \
            const ValueType,                                            \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::binaryOpSparse(                          \
            hmemo::HArray<IndexType>&,                                  \
            hmemo::HArray<ValueType>&,                                  \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const ValueType,                                            \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const ValueType,                                            \
            const BinaryOp,                                             \
            hmemo::ContextPtr );                                        \
    template IndexType HArrayUtils::allSparse(                          \
            bool&,                                                      \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const ValueType,                                            \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const ValueType,                                            \
            const CompareOp,                                            \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::mergeSparse(                             \
            hmemo::HArray<IndexType>&,                                  \
            hmemo::HArray<ValueType>&,                                  \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<IndexType>&,                            \
            const hmemo::HArray<ValueType>&,                            \
            const BinaryOp,                                             \
            hmemo::ContextPtr );                                        \
    template void HArrayUtils::mergeAndMap(                             \
        hmemo::HArray<ValueType> & result,                              \
        hmemo::HArray<IndexType> & xMap,                                \
        hmemo::HArray<IndexType> & yMap,                                \
        const hmemo::HArray<ValueType> & x,                             \
        const hmemo::HArray<ValueType> & y,                             \
        const common::CompareOp comparator,                             \
        hmemo::ContextPtr prefLoc );                                    \
    SCAI_COMMON_LOOP_LVL2( ValueType, HARRAUTILS_SPECIFIER_LVL2, SCAI_ARRAY_TYPES_HOST )

SCAI_COMMON_LOOP( HARRAYUTILS_SPECIFIER, SCAI_ARRAY_TYPES_HOST )

#undef HARRAYUTILS_SPECIFIER
#undef HARRAUTILS_SPECIFIER_LVL2

// The UnaryOp/binary operations are instantiated for all array types
// But keep in mind that some operations are not supported for IndexType
// and will therefore cause runtime errors, i.e. throw an exception

#define HARRAYUTILS_SPECIFIER( ValueType )                        \
    template void HArrayUtils::unaryOp<ValueType>(                \
            hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                      \
            const UnaryOp, hmemo::ContextPtr);                    \
    template void HArrayUtils::binaryOp<ValueType>(               \
            hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                      \
            const hmemo::HArray<ValueType>&,                      \
            const BinaryOp,                                       \
            hmemo::ContextPtr);                                   \
    template void HArrayUtils::binaryOpScalar<ValueType>(         \
            hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                      \
            const ValueType,                                      \
            const BinaryOp,                                       \
            const bool,                                           \
            hmemo::ContextPtr);

// Note: it is now safe to instantiate UnaryOp/binary kernels also for IndexType

SCAI_COMMON_LOOP( HARRAYUTILS_SPECIFIER, SCAI_ARRAY_TYPES_HOST )

#undef HARRAYUTILS_SPECIFIER

#define HARRAYUTILS_SPECIFIER( ValueType )                        \
    template void HArrayUtils::selectComplexPart(                 \
            hmemo::HArray<RealType<ValueType> >&,                 \
            const hmemo::HArray<ValueType>&,                      \
            const common::ComplexPart,                            \
            hmemo::ContextPtr );                                  \
    template void HArrayUtils::buildComplex(                      \
            hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<RealType<ValueType> >&,           \
            const hmemo::HArray<RealType<ValueType> >&,           \
            hmemo::ContextPtr );                                  \

// selectComplexPart uses Math::real and Math::imag that is not defined for IndexType

SCAI_COMMON_LOOP( HARRAYUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef HARRAYUTILS_SPECIFIER

template void HArrayUtils::bucketSort(
    hmemo::HArray<IndexType>& offsets,
    hmemo::HArray<IndexType>& perm,
    const hmemo::HArray<IndexType>& array,
    const IndexType nb,
    hmemo::ContextPtr prefLoc );

template void HArrayUtils::bucketCount(
    hmemo::HArray<IndexType>& sizes,
    const hmemo::HArray<IndexType>& array,
    const IndexType nb,
    hmemo::ContextPtr prefLoc );

/* --------------------------------------------------------------------------- */

} /* end namespace utilskernel */

} /* end namespace scai */
