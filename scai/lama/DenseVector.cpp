/**
 * @file DenseVector.cpp
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
 * @brief Implementations and instantiations for class DenseVector.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/DenseVector.hpp>

// local library
#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/expression/Expression.hpp>

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/PartitionIO.hpp>

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>
#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>

// std
#include <ostream>

namespace scai
{

using common::scoped_array;
using common::TypeTraits;
using utilskernel::HArrayUtils;
using utilskernel::LArray;

using namespace hmemo;
using namespace dmemo;

namespace lama
{

// common logger for all types of DenseVector

SCAI_LOG_DEF_LOGGER( _DenseVector::logger, "Vector.DenseVector" )

_DenseVector::_DenseVector( const IndexType n ) :

    Vector( n )
{
}

_DenseVector::_DenseVector( const IndexType n, ContextPtr context ) :

    Vector( n, context )
{
}

_DenseVector::_DenseVector( DistributionPtr dist ) :

    Vector( dist )
{
}

_DenseVector::_DenseVector( DistributionPtr dist, ContextPtr context ) :

    Vector( dist, context )
{
}

_DenseVector::_DenseVector( const _DenseVector& other ) :

    Vector( other )
{
}

_DenseVector::_DenseVector( const Vector& other ) :

    Vector( other )
{
}

_DenseVector* _DenseVector::create( common::scalar::ScalarType type )
{
    // There is only one factor for all vectors

    Vector* v = Vector::create( VectorCreateKeyType( Vector::DENSE, type ) );

    // reinterpret cast is safe

    return reinterpret_cast<_DenseVector*>( v );
}

/* ------------------------------------------------------------------------- */
/*  Constructors of DenseVector<ValueType>                                   */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector() :
    _DenseVector( 0 ),
    mLocalValues()
{
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( ContextPtr context ) :

    _DenseVector( 0, context ),
    mLocalValues()
{
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const IndexType size ) :

    _DenseVector( size ), 
    mLocalValues( size )

{
    SCAI_LOG_INFO( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">( size = " 
                            << size << " ), undefined values" )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution ) :

    _DenseVector( distribution ), 
    mLocalValues( distribution->getLocalSize() )
{
    SCAI_LOG_INFO( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">( dist = " 
                   << *distribution << "), undefined values" )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const IndexType size, ContextPtr context ) : 

    _DenseVector( size, context ), 
    mLocalValues( size )

{
    SCAI_LOG_INFO( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">( size = " 
                            << size << " ), undefined values, @ctx = " << *mContext )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, ContextPtr context ) : 

    _DenseVector( distribution, context ), 
    mLocalValues( distribution->getLocalSize() )

{
    SCAI_LOG_INFO( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">( dist = " 
                   << *distribution << "), undefined values, @ctx = " << *mContext )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const IndexType size, const ValueType value, ContextPtr context ) :

    _DenseVector( size, context ), 
    mLocalValues( size, value, context )

{
    SCAI_LOG_INFO( logger, "Construct dense vector, size = " << size << ", init =" << value )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, const ValueType value, ContextPtr context ) :

    _DenseVector( distribution, context ), 
    mLocalValues( distribution->getLocalSize(), value )

{
    SCAI_LOG_INFO( logger,
                   "Construct dense vector, size = " << distribution->getGlobalSize() << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", value = " << value )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const IndexType size, const ValueType startValue, const ValueType inc, ContextPtr context ) :

    _DenseVector( size, context ), 
    mLocalValues( size, startValue, inc, context )

{
    SCAI_LOG_INFO( logger, "Construct dense vector, size = " << size << ", startValue =" << startValue << ", inc=" << inc )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, const ValueType startValue, const ValueType inc, ContextPtr context ) :

    _DenseVector( distribution, context ),
    mLocalValues( distribution->getLocalSize(), context )

{
    allocate();    // correct allocation of mLocalValues

    SCAI_LOG_INFO( logger,
                   "Construct dense vector, size = " << distribution->getGlobalSize() << ", distribution = " << *distribution
                   << ", local size = " << distribution->getLocalSize() << ", startValue = " << startValue << ", inc=" << inc )

    // get my owned indexes

    HArray<IndexType> myGlobalIndexes( context );

    // mult with inc and add startValue

    distribution->getOwnedIndexes( myGlobalIndexes );

    // localValues[] =  indexes[] * inc + startValue

    HArrayUtils::assign( mLocalValues, myGlobalIndexes, context );
    HArrayUtils::assignScalar( mLocalValues, inc, utilskernel::binary::MULT, context );
    HArrayUtils::assignScalar( mLocalValues, startValue, utilskernel::binary::ADD, context );
}

template <typename ValueType>
void DenseVector<ValueType>::setSequence( const Scalar startValue, const Scalar inc, const IndexType n )
{
    setDistributionPtr( DistributionPtr( new NoDistribution( n ) ) );

    HArrayUtils::setSequence( mLocalValues, startValue.getValue<ValueType>(), inc.getValue<ValueType>(), n, getContextPtr() );
}

template <typename ValueType>
void DenseVector<ValueType>::setSequence( const Scalar startValue, const Scalar inc, DistributionPtr distribution )
{
    setDistributionPtr( distribution );

    if ( distribution->isReplicated() )
    {
        SCAI_ASSERT_EQ_DEBUG( distribution->getGlobalSize(), distribution->getLocalSize(), *distribution << " not replicated" );

        HArrayUtils::setSequence( mLocalValues, startValue.getValue<ValueType>(), inc.getValue<ValueType>(), distribution->getGlobalSize() );
        return;
    }

    ContextPtr context = getContextPtr();

    // get my owned indexes

    HArray<IndexType> myGlobalIndexes( context );

    // mult with inc and add startValue

    distribution->getOwnedIndexes( myGlobalIndexes );

    // localValues[] =  indexes[] * inc + startValue

    HArrayUtils::assign( mLocalValues, myGlobalIndexes, context );
    HArrayUtils::assignScalar( mLocalValues, inc.getValue<ValueType>(), utilskernel::binary::MULT, context );
    HArrayUtils::assignScalar( mLocalValues, startValue.getValue<ValueType>(), utilskernel::binary::ADD, context );
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Vector& other ) :

    _DenseVector( other )

{
    allocate();
    assign( other );
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Vector& other, DistributionPtr distribution ) :

    _DenseVector( other )

{
    allocate();       // make sure that local values array fits distribution, context
    assign( other );
    redistribute( distribution );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const std::string& filename ) : 

    _DenseVector( 0 )

{
    SCAI_LOG_INFO( logger, "Construct dense vector from file " << filename )
    readFromFile( filename );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::setRandom( dmemo::DistributionPtr distribution, const float fillRate )
{
    allocate( distribution );
    mLocalValues.setRandom( mLocalValues.size(), fillRate, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const _HArray& localValues, DistributionPtr distribution ) :

    _DenseVector( distribution ),
    mLocalValues( localValues )
{
    SCAI_ASSERT_EQ_ERROR( localValues.size(), distribution->getLocalSize(), "size mismatch" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const _HArray& values ) :

    _DenseVector( DistributionPtr( new NoDistribution( values.size() ) ) ),
    mLocalValues( values )

{
    SCAI_LOG_DEBUG( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() 
                             << ">( HArray<" << values.getValueType() << ">, size = " << values.size() << " )" )
}

/* ------------------------------------------------------------------------- */

/*
 * Constructors with Expressions as arguments
 */

// linear algebra expression: a*x
template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SV& expression ):

    _DenseVector( expression.getArg2() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a+x/x+a
template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SV_S& expression ) :

    _DenseVector( expression.getArg1().getArg2() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x + beta)" )
    Vector::operator=( expression );
}

// linear algebra expression: x*y

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_VV& expression ) :

    _DenseVector( expression.getArg1() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( x * y )" )
    Expression_SVV tmpExp( Scalar( 1 ), expression );
    Vector::operator=( tmpExp );
}

// linear algebra expression: s*x*y
template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SVV& expression ) :

    _DenseVector( expression.getArg2().getArg1() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x+b*y, inherit distribution/context from vector x

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SV_SV& expression ) :

    _DenseVector( expression.getArg1().getArg2() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x + beta * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SMV_SV& expression ) :

    _DenseVector( expression.getArg1().getArg2().getArg1().getRowDistributionPtr(),
                  expression.getArg1().getArg2().getArg1().getContextPtr() )
{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SVM_SV& expression ) :

    _DenseVector( expression.getArg1().getArg2().getArg2().getColDistributionPtr(),
                  expression.getArg1().getArg2().getArg2().getContextPtr() )
{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SMV& expression )

    : _DenseVector( expression.getArg2().getArg1().getRowDistributionPtr(),
                    expression.getArg2().getArg1().getContextPtr() )
{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x*A, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SVM& expression ) :

    _DenseVector( expression.getArg2().getArg2().getColDistributionPtr(),
                  expression.getArg2().getArg2().getContextPtr() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A )" )
    Vector::operator=( expression );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::~DenseVector()
{
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>& DenseVector<ValueType>::operator=( const DenseVector<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "DenseVector<" << TypeTraits<ValueType>::id() << "> = " <<
                   "DenseVector<" << TypeTraits<ValueType>::id() << ">" )

    assign( other );
    return *this;
}

template<typename ValueType>
DenseVector<ValueType>& DenseVector<ValueType>::operator=( const Scalar value )
{
    SCAI_LOG_INFO( logger, "DenseVector<" << TypeTraits<ValueType>::id() << "> = " << value )

    assign( value );
    return *this;
}

/* ------------------------------------------------------------------------- */

/** Determine splitting values for sorting distributed values.
 *
 *  A value v belongs to partition p if splitValues[p] <= v < splitValues[p+1]
 */
template<typename ValueType>
static void getSplitValues(
    ValueType splitValues[],
    const Communicator& comm,
    const ValueType sortedValues[],
    const IndexType n,
    const bool ascending )
{
    const PartitionId nPartitions = comm.getSize();

    if ( ascending )
    {
        ValueType minV = n > 0 ? sortedValues[0] : TypeTraits<ValueType>::getMax();
        ValueType maxV = n > 0 ? sortedValues[n - 1] : TypeTraits<ValueType>::getMin();

        splitValues[0]           = comm.min( minV );
        splitValues[nPartitions] = comm.max( maxV );
    }
    else
    {
        ValueType maxV = n > 0 ? sortedValues[0] : TypeTraits<ValueType>::getMin();
        ValueType minV = n > 0 ? sortedValues[n - 1] : TypeTraits<ValueType>::getMax();

        splitValues[0]           = comm.max( maxV );
        splitValues[nPartitions] = comm.min( minV );
    }

    // fill intermediate values by uniform distribution of range splitValues[0] .. splitValues[nPartitions]

    for ( PartitionId p = 1; p < nPartitions; ++p )
    {
        splitValues[p] = splitValues[0] + ( splitValues[nPartitions] - splitValues[0] ) * ValueType( p ) / ValueType( nPartitions );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::sort( bool ascending )
{
    sortImpl( NULL, this, *this, ascending );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::sort( DenseVector<IndexType>& perm, bool ascending )
{
    sortImpl( &perm, this, *this, ascending );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool DenseVector<ValueType>::isSorted( bool ascending ) const
{
    const dmemo::Distribution& distribution = getDistribution();
    const dmemo::Communicator& comm         = distribution.getCommunicator();

    PartitionId rank = comm.getRank();
    PartitionId size = comm.getSize();

    IndexType blockSize = distribution.getBlockDistributionSize();

    SCAI_ASSERT_NE_ERROR( blockSize, nIndex, "isSorted only possible on block distributed vectors" )

    const HArray<ValueType>& localValues = getLocalValues();

    bool is = HArrayUtils::isSorted( localValues, ascending );

    if ( size == 1 )
    {
        return is;    // we are done
    }

    SCAI_LOG_DEBUG( logger, "local array is sorted = " << is )

    const IndexType n = localValues.size();

    ValueType tmp;

    if ( n == 0 )
    {
        comm.shift( &tmp, IndexType( 1 ), &tmp, IndexType( 0 ), 1 );
        comm.shift( &tmp, IndexType( 1 ), &tmp, IndexType( 0 ), -1 );

        is = comm.all( is );

        return is;
    }

    ReadAccess<ValueType> rLocalValues( localValues );

    // send right value to right neighbor

    IndexType nt = comm.shift( &tmp, IndexType( 1 ), &rLocalValues[n - 1], IndexType( 1 ), 1 );

    if ( nt && rank > 0 )
    {
        if ( ascending )
        {
            if ( tmp > rLocalValues[0] )
            {
                is = false;
                SCAI_LOG_INFO( logger, comm << ": value of left neighbor " << tmp << " greater than my lowest value " << rLocalValues[0] )
            }
        }
        else
        {
            if ( tmp < rLocalValues[0] )
            {
                is = false;
                SCAI_LOG_INFO( logger, comm << ": value of left neighbor " << tmp << " less than my greatest value " << rLocalValues[0] )
            }
        }
    }

    // send left value to left neighbor

    nt = comm.shift( &tmp, IndexType( 1 ), &rLocalValues[0], IndexType( 1 ), -1 );

    if ( nt && rank < size - 1 )
    {
        if ( ascending )
        {
            if ( tmp < rLocalValues[n - 1] )
            {
                is = false;
                SCAI_LOG_INFO( logger, comm << ": value of right neighbor " << tmp << " less than my highest value " << rLocalValues[n - 1] )
            }
        }
        else
        {
            if ( tmp > rLocalValues[n - 1] )
            {
                is = false;
                SCAI_LOG_INFO( logger, comm << ": value of right neighbor " << tmp << " greater than my highest value " << rLocalValues[n - 1] )
            }
        }
    }

    is = comm.all( is );

    return is;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::sortImpl(
    DenseVector<IndexType>* perm,
    DenseVector<ValueType>* out,
    DenseVector<ValueType>& in,
    bool ascending )
{
    const dmemo::Distribution& distribution = in.getDistribution();

    const dmemo::Communicator& comm = distribution.getCommunicator();

    IndexType blockSize = distribution.getBlockDistributionSize();

    SCAI_ASSERT_NE_ERROR( blockSize, nIndex, "sort only possible on block distributed vectors" )

    // Now sort the local values

    LArray<IndexType>* localPerm = NULL;

    if ( perm )
    {
        localPerm = & perm->getLocalValues();
    }

    const HArray<ValueType>&  inValues = in.getLocalValues();

    SCAI_ASSERT_ERROR( out, "Optional argument out, must be present" )

    HArray<ValueType>& sortedValues = out->getLocalValues();

    utilskernel::HArrayUtils::sort( localPerm, &sortedValues, inValues, ascending );

    if ( localPerm )
    {
        const IndexType nLocalPerm = localPerm->size();
        SCAI_ASSERT_EQ_ERROR( nLocalPerm, inValues.size(), "size mismatch for perm array from sort" );

        // the local indexes of permutation must be translated to global indexes

        if ( nLocalPerm > 0 )
        {
            // due to block distribution we need only global index of first one

            *localPerm += distribution.local2global( 0 );
        }

        ReadAccess<IndexType> rPerm( *localPerm );

        for ( IndexType i = 0; i < nLocalPerm; ++i )
        {
            SCAI_LOG_TRACE( logger, "localPerm[" << i << "] = " << rPerm[i] )
        }
    }

    // Determine the splitting values

    PartitionId nPartitions = comm.getSize();

    common::scoped_array<ValueType> splitValues( new ValueType[ nPartitions + 1 ] );

    {
        ReadAccess<ValueType> rSortedValues( sortedValues );
        getSplitValues( splitValues.get(), comm, rSortedValues.get(), sortedValues.size(), ascending );
    }

    common::scoped_array<IndexType> quantities( new IndexType[ nPartitions ] );

    // Determine quantities for each processor

    for ( PartitionId ip = 0; ip < nPartitions; ++ip )
    {
        quantities[ip] = 0;
    }

    PartitionId p = 0;

    {
        ReadAccess<ValueType> rValues( sortedValues );

        for ( IndexType i = 0; i < blockSize; ++i )
        {
            if ( ascending )
            {
                while ( rValues[i] > splitValues[p + 1] )
                {
                    p++;
                    SCAI_ASSERT_LT_DEBUG( p, nPartitions, "Illegal split values, mLocalValues[" << i << "] = " << rValues[i] );
                }
            }
            else
            {
                while ( rValues[i] < splitValues[p + 1] )
                {
                    p++;
                    SCAI_ASSERT_LT_DEBUG( p, nPartitions, "Illegal split values, mLocalValues[" << i << "] = " << rValues[i] );
                }
            }

            quantities[p]++;
        }
    }

    // make communication plans for sending data and receiving data

    dmemo::CommunicationPlan sendPlan;

    sendPlan.allocate( quantities.get(), nPartitions );

    dmemo::CommunicationPlan recvPlan;

    recvPlan.allocateTranspose( sendPlan, comm );

    SCAI_LOG_INFO( logger, comm << ": send plan: " << sendPlan << ", rev plan: " << recvPlan );

    LArray<ValueType> newValues;

    IndexType newLocalSize = recvPlan.totalQuantity();

    {
        WriteOnlyAccess<ValueType> recvVals( newValues, newLocalSize );
        ReadAccess<ValueType> sendVals( sortedValues );
        comm.exchangeByPlan( recvVals.get(), recvPlan, sendVals.get(), sendPlan );
    }

    // Also communicate the original index positions of the array values

    if ( perm )
    {
        LArray<IndexType> newPerm;

        {
            WriteOnlyAccess<IndexType> recvVals( newPerm, newLocalSize );
            ReadAccess<IndexType> sendVals( *localPerm );
            comm.exchangeByPlan( recvVals.get(), recvPlan, sendVals.get(), sendPlan );
        }

        // accesses must be released before swapping of arrays

        localPerm->swap( newPerm );

        ReadAccess<IndexType> rPerm( *localPerm );

        for ( IndexType i = 0; i < localPerm->size(); ++i )
        {
            SCAI_LOG_TRACE( logger, comm << ": newPerm[" << i << "] = " << rPerm[i] )
        }
    }

    // Merge the values received from other processors

    HArray<IndexType> mergeOffsets;   // offsets for sorted subarrays in mergesort

    {
        IndexType nOffsets = recvPlan.size();

        WriteOnlyAccess<IndexType> wOffsets( mergeOffsets, nOffsets + 1 );

        wOffsets[0] = 0;

        for ( IndexType k = 0; k < nOffsets; ++k )
        {
            wOffsets[k + 1] = wOffsets[k] + recvPlan[k].quantity;
        }
    }

    if ( perm )
    {
        utilskernel::HArrayUtils::mergeSort( newValues, *localPerm, mergeOffsets, ascending );
    }
    else
    {
        utilskernel::HArrayUtils::mergeSort( newValues, mergeOffsets, ascending );
    }

    // Create the new general block distribution

    DistributionPtr newDist( new dmemo::GenBlockDistribution( distribution.getGlobalSize(), newLocalSize, distribution.getCommunicatorPtr() ) );

    if ( out )
    {
        out->swap( newValues, newDist );
    }

    if ( perm )
    {
        perm->swap( *localPerm, newDist );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::swap( HArray<ValueType>& newValues, DistributionPtr newDist )
{
    SCAI_LOG_DEBUG( logger, *this << ": swap with new local values = " << newValues << ", new dist = " << *newDist )

    SCAI_ASSERT_EQ_ERROR( newValues.size(), newDist->getLocalSize(), "serious mismatch" )

    mLocalValues.swap( newValues );
    setDistributionPtr( newDist );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::scalar::ScalarType DenseVector<ValueType>::getValueType() const
{
    return TypeTraits<ValueType>::stype;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::setDenseValues( const _HArray& values )
{
    const IndexType localSize = getDistribution().getLocalSize();

    SCAI_ASSERT_EQ_ERROR( localSize, values.size(), "size of array with local values does not match local size of distribution" )

    HArrayUtils::assign( mLocalValues, values, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::setSparseValues( const HArray<IndexType>& nonZeroIndexes, const _HArray& nonZeroValues )
{
    const IndexType localSize = getDistribution().getLocalSize();

    // Note: buildDenseArray checks for legal indexes, order of indexes does not matter

    HArrayUtils::buildDenseArray( mLocalValues, localSize, nonZeroValues, nonZeroIndexes, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>* DenseVector<ValueType>::copy() const
{
    // create a new dense vector with the copy constructor
    return new DenseVector<ValueType>( *this );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>* DenseVector<ValueType>::newVector() const
{
    common::unique_ptr<DenseVector<ValueType> > vector( new DenseVector<ValueType>() );
    vector->setContextPtr( this->getContextPtr() );
    return vector.release();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::getValue( IndexType globalIndex ) const
{
    ValueType myValue = 0;

    const IndexType localIndex = getDistribution().global2local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": getValue( globalIndex = " << globalIndex << " ) -> local : " << localIndex )

    if ( localIndex != nIndex )
    {
        myValue = mLocalValues[localIndex];
    }

    ValueType allValue = getDistribution().getCommunicator().sum( myValue );

    // works also fine for replicated distributions with NoCommunicator

    SCAI_LOG_TRACE( logger, "myValue = " << myValue << ", allValue = " << allValue )

    return Scalar( allValue );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::setValue( const IndexType globalIndex, const Scalar value )
{
    SCAI_LOG_TRACE( logger, *this << ": setValue( globalIndex = " << globalIndex << " ) = " <<  value )

    const IndexType localIndex = getDistribution().global2local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": set @g " << globalIndex << " is @l " << localIndex << " : " << value )

    if ( localIndex != nIndex )
    {
        mLocalValues[localIndex] = value.getValue<ValueType>();
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::min() const
{
    // Note: min returns the maximal representation value on zero-sized vectors, TypeTraits<ValueType>::getMax()
    ValueType localMin = mLocalValues.min();
    return Scalar( getDistribution().getCommunicator().min( localMin ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::max() const
{
    // Note: max returns the minimal representation value on zero-sized vectors
    ValueType localMax = mLocalValues.max();
    return Scalar( getDistribution().getCommunicator().max( localMax ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::l1Norm() const
{
    ValueType localL1Norm = mLocalValues.l1Norm();
    return Scalar( getDistribution().getCommunicator().sum( localL1Norm ) );
}

/*---------------------------------------------------------------------------*/
template<typename ValueType>
Scalar DenseVector<ValueType>::sum() const
{
    ValueType localsum = mLocalValues.sum();
    return Scalar( getDistribution().getCommunicator().sum( localsum ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::l2Norm() const
{
    // Note: we do not call l2Norm here for mLocalValues to avoid sqrt
    ValueType localDotProduct = mLocalValues.dotProduct( mLocalValues );
    ValueType globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    return Scalar( common::Math::sqrt( globalDotProduct ) );
}

template<>
Scalar DenseVector<IndexType>::l2Norm() const
{
    // Note: we do not call l2Norm here for mLocalValues to avoid sqrt

    ScalarRepType localDotProduct = mLocalValues.dotProduct( mLocalValues );
    ScalarRepType globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    return Scalar( common::Math::sqrt( globalDotProduct ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::maxNorm() const
{
    ValueType localMaxNorm = mLocalValues.maxNorm();
    const Communicator& comm = getDistribution().getCommunicator();
    ValueType globalMaxNorm = comm.max( localMaxNorm );
    SCAI_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm: " << localMaxNorm
                   << ", max norm global = " << globalMaxNorm )
    return Scalar( globalMaxNorm );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::swap( Vector& other )
{
    SCAI_LOG_DEBUG( logger, "swap:" << *this << " with " << other )
    DenseVector* otherPtr = dynamic_cast<DenseVector*>( &other );

    if ( !otherPtr )
    {
        COMMON_THROWEXCEPTION( "Tried to swap with a Vector of a different type." )
    }

    Vector::swapVector( other );
    mLocalValues.swap( otherPtr->mLocalValues );
    // mHaloValues.swap( otherPtr->mHaloValues );
}

template<typename ValueType>
void DenseVector<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "DenseVector<" << getValueType() << ">" << "( size = " << size() << ", local = " << mLocalValues.size()
           << ", dist = " << getDistribution() << ", loc  = " << *getContextPtr() << " )";
}

template<typename ValueType>
void DenseVector<ValueType>::assignScaledVector( const Scalar& alpha, const Vector& x )
{
    SCAI_LOG_INFO( logger, "assignScaledVector: this = " << alpha << " * x, with x = " << x << ", this = " << *this )

    const ValueType alphaV = alpha.getValue<ValueType>();

    if ( x.getVectorKind() != DENSE || x.getValueType() != getValueType() )
    {
        SCAI_LOG_INFO( logger, "this = " << alpha << " * x -> this = x;  this *= " << alpha )
        assign( x );
        HArrayUtils::setScalar( mLocalValues, alphaV, utilskernel::binary::MULT, mContext );
        return;
    }

    // make sure that in case of alias x == *this we do not reallocate

    SCAI_LOG_DEBUG( logger, "distribution X = " << x.getDistribution() << ", distribution this = " << getDistribution() )

    if ( x.getDistribution() != getDistribution() )
    {
        allocate( x.getDistributionPtr() );
    }

    const DenseVector<ValueType>& denseX = dynamic_cast<const DenseVector<ValueType>&>( x );

    SCAI_LOG_INFO( logger, "x = " << denseX << ", this = " << *this )

    SCAI_ASSERT_EQ_ERROR( mLocalValues.size(), denseX.mLocalValues.size(), "local size mismatch" )

    SCAI_LOG_DEBUG( logger, "call arrayPlusArray" )

    ValueType betaV = 0;

    utilskernel::HArrayUtils::arrayPlusArray( mLocalValues, alphaV, denseX.mLocalValues, betaV, denseX.mLocalValues, mContext );
}

template<typename ValueType>
void DenseVector<ValueType>::vectorPlusVector( const Scalar& alpha, const Vector& x, const Scalar& beta, const Vector& y )
{
    SCAI_LOG_INFO( logger, "vectorPlusVector: z = " << alpha << " * x + " << beta << " * y, with  x = " << x << ", y = " << y << ", z = " << *this )

    // Query for alpha == 0 or beta == 0 as x and y might be undefined 

    if ( alpha == Scalar( 0 ) )
    {
        assignScaledVector( beta, y );
        return;
    }

    if ( beta == Scalar( 0 ) )
    {
        assignScaledVector( alpha, x );
        return;
    }

    if ( x.getVectorKind() != DENSE || x.getValueType() != getValueType() )
    {
        SCAI_LOG_WARN( logger, "vectorPlusVector: use temporary DenseVector<" << getValueType() << "> for x = " << x )
        DenseVector<ValueType> xTmp( x );
        vectorPlusVector( alpha, xTmp, beta, y );
        return;
    }

    if ( y.getVectorKind() != DENSE || y.getValueType() != getValueType() )
    {
        SCAI_LOG_WARN( logger, "vectorPlusVector: use temporary DenseVector<" << getValueType() << "> for y = " << y )
        DenseVector<ValueType> yTmp( y );
        vectorPlusVector( alpha, x, beta, yTmp );
        return;
    }

    SCAI_REGION( "Vector.Dense.VplusV" )

    const ValueType alphaV = alpha.getValue<ValueType>();
    const ValueType betaV  = beta.getValue<ValueType>();

    if ( x.getDistribution() != y.getDistribution() )
    {
        COMMON_THROWEXCEPTION(
            "distribution do not match for z = alpha * x + beta * y, z = " << *this << " , x = " << x << " , y = " << y )
    }

    if ( x.getDistribution() != getDistribution() || x.size() != size() )
    {
        allocate( x.getDistributionPtr() );
    }

    const DenseVector<ValueType>& denseX = dynamic_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>& denseY = dynamic_cast<const DenseVector<ValueType>&>( y );

    if ( mLocalValues.size() != denseX.mLocalValues.size() )
    {
        SCAI_LOG_DEBUG( logger, "resize local values of z = this" )
        mLocalValues.clear();
        mLocalValues.resize( denseX.mLocalValues.size() );
    }

    SCAI_LOG_DEBUG( logger, "call arrayPlusArray" )
    utilskernel::HArrayUtils::arrayPlusArray( mLocalValues, alphaV, denseX.mLocalValues, betaV, denseY.mLocalValues, mContext );
}

template<typename ValueType>
void DenseVector<ValueType>::vectorTimesVector( const Scalar& alpha, const Vector& x, const Vector& y )
{
    SCAI_LOG_INFO( logger, "z = x * y, z = " << *this << " , x = " << x << " , y = " << y )
    SCAI_LOG_DEBUG( logger, "dist of x = " << x.getDistribution() )
    SCAI_LOG_DEBUG( logger, "dist of y = " << y.getDistribution() )

    if ( x.getDistribution() != y.getDistribution() )
    {
        COMMON_THROWEXCEPTION(
            "distribution do not match for z = x * y, z = " << *this << " , x = " << x << " , y = " << y )
    }

    if ( x.getDistribution() != getDistribution() || x.size() != size() )
    {
        allocate( x.getDistributionPtr() );
    }

    if ( x.getVectorKind() != DENSE || x.getValueType() != getValueType() )
    {
        DenseVector<ValueType> xTmp( x );
        vectorTimesVector( alpha, xTmp, y );
        return;
    }

    if ( y.getVectorKind() != DENSE || y.getValueType() != getValueType() )
    {
        DenseVector<ValueType> yTmp( y );
        vectorTimesVector( alpha, x, yTmp );
        return;
    }

    const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>& denseY = reinterpret_cast<const DenseVector<ValueType>&>( y );

    mLocalValues.resize( denseX.mLocalValues.size() );

    SCAI_LOG_DEBUG( logger, "call arrayTimesArray" )

    ValueType alphaV = alpha.getValue<ValueType>();

    utilskernel::HArrayUtils::arrayTimesArray( mLocalValues, alphaV, denseX.mLocalValues, denseY.mLocalValues, mContext );
}

template<typename ValueType>
void DenseVector<ValueType>::vectorPlusScalar( const Scalar& alpha, const Vector& x, const Scalar& beta )
{
    if ( x.getVectorKind() != DENSE || x.getValueType() != getValueType() )
    {
        SCAI_LOG_WARN( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">::vectorAddScalar, uses tmp for x" )

        DenseVector<ValueType> xTmp( x );
        vectorPlusScalar( alpha, xTmp, beta );
        return;
    }

    const DenseVector<ValueType>& denseX = dynamic_cast<const DenseVector<ValueType>&>( x );

    const ValueType alphaV = alpha.getValue<ValueType>();
    const ValueType betaV = beta.getValue<ValueType>();

    SCAI_LOG_INFO( logger, "z = alpha * x + beta, z = " << *this << ", alpha=  " << alpha
                   << " , x = " << x << " , beta = " << beta )
    SCAI_LOG_DEBUG( logger, "dist of x = " << x.getDistribution() )

    if ( x.getDistribution() != getDistribution() || x.size() != size() )
    {
        allocate( x.getDistributionPtr() );
    }

    if ( mLocalValues.size() != denseX.mLocalValues.size() )
    {
        mLocalValues.clear();
        mLocalValues.resize( denseX.mLocalValues.size() );
    }

    SCAI_LOG_DEBUG( logger, "call arrayPlusScalar" )

    utilskernel::HArrayUtils::arrayPlusScalar( mLocalValues, alphaV, denseX.mLocalValues, betaV, mContext );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::gather(
    const DenseVector<ValueType>& source,
    const DenseVector<IndexType>& index,
    const utilskernel::binary::BinaryOp op )
{
    if ( op != utilskernel::binary::COPY )
    {
        SCAI_ASSERT_EQ_ERROR( getDistribution(), index.getDistribution(), "both vectors must have same distribution" )
    }

    const Distribution& sourceDistribution = source.getDistribution();

    if ( sourceDistribution.isReplicated() )
    {
        // we can just do it locally and are done

        HArrayUtils::gather( mLocalValues, source.getLocalValues(), index.getLocalValues(), op );

        assign( mLocalValues, index.getDistributionPtr() );

        return;
    }

    const Communicator& comm = sourceDistribution.getCommunicator();

    // otherwise we have to set up a communication plan

    HArray<PartitionId> owners;

    sourceDistribution.computeOwners( owners, index.getLocalValues() );

    // set up required values by sorting the indexes corresponding to the owners via bucketsort

    const PartitionId size = comm.getSize();

    HArray<IndexType> offsets;  // used to allocate the communication plan

    HArray<IndexType> perm;     // used to sort required indexes and to scatter the gathered values

    HArrayUtils::bucketSort( offsets, perm, owners, size );

    SCAI_ASSERT_EQ_DEBUG( offsets.size(), size + 1, "wrong offsets" )
    SCAI_ASSERT_EQ_DEBUG( perm.size(), owners.size(), "illegal perm" )

    HArray<IndexType> requiredIndexes;  // local index values sorted by owner

    HArrayUtils::gather( requiredIndexes, index.getLocalValues(), perm, utilskernel::binary::COPY );

    // exchange communication plans

    dmemo::CommunicationPlan recvPlan;

    dmemo::CommunicationPlan sendPlan;

    {
        hmemo::ReadAccess<IndexType> rOffsets( offsets );
        recvPlan.allocateByOffsets( rOffsets.get(), size );
    }

    sendPlan.allocateTranspose( recvPlan, comm );

    SCAI_LOG_DEBUG( logger, comm << ": recvPlan = " << recvPlan << ", sendPlan = " << sendPlan )

    HArray<IndexType> sendIndexes;

    comm.exchangeByPlan( sendIndexes, sendPlan, requiredIndexes, recvPlan );

    // translate global sendIndexes to local indexes, all must be local

    {
        WriteAccess<IndexType> wSendIndexes( sendIndexes );

        for ( IndexType i = 0; i < sendIndexes.size(); ++i )
        {
            IndexType localIndex = sourceDistribution.global2local( wSendIndexes[i] );
            SCAI_ASSERT_NE_DEBUG( localIndex, nIndex, "got required index " << wSendIndexes[i] << " but I'm not owner" )
            wSendIndexes[i] = localIndex;
        }
    }

    // exchange communication plan

    HArray<ValueType> sendValues;  // values to send from my source values

    HArrayUtils::gather( sendValues, source.getLocalValues(), sendIndexes, utilskernel::binary::COPY );

    // send via communication plan

    HArray<ValueType> recvValues;

    comm.exchangeByPlan( recvValues, recvPlan, sendValues, sendPlan );

    if ( op == utilskernel::binary::COPY )
    {
        mLocalValues.resize( perm.size() );
    }

    // required indexes were sorted according to perm, using inverse perm here via scatter

    bool isUnique = false;

    HArrayUtils::scatter( mLocalValues, perm, isUnique, recvValues, op );

    assign( mLocalValues, index.getDistributionPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::scatter(
    const DenseVector<IndexType>& index,
    const DenseVector<ValueType>& source,
    const utilskernel::binary::BinaryOp op )
{
    bool hasUniqueIndexes = false;

    SCAI_ASSERT_EQ_ERROR( source.getDistribution(), index.getDistribution(), "both vectors must have same distribution" )

    // get the owners of my local index values, relevant is distribution of this target vector

    const Distribution& targetDistribution = getDistribution();

    SCAI_ASSERT_EQ_ERROR( source.getDistribution().isReplicated(), targetDistribution.isReplicated(),
                          "scatter: either all arrays are replicated or all are distributed" )

    if ( targetDistribution.isReplicated() )
    {
        // so all involved arrays are replicated, we can do it just locally

        HArrayUtils::scatter( mLocalValues, index.getLocalValues(), hasUniqueIndexes, source.getLocalValues(), op );

        return;
    }

    const Communicator& comm = targetDistribution.getCommunicator();

    HArray<PartitionId> owners;

    targetDistribution.computeOwners( owners, index.getLocalValues() );

    // set up provided values by sorting the indexes corresponding to the owners via bucketsort

    const PartitionId size = comm.getSize();

    HArray<IndexType> offsets;  // used to allocate the communication plan

    HArray<IndexType> perm;     // used to sort required indexes and to scatter the gathered values

    HArrayUtils::bucketSort( offsets, perm, owners, size );

    SCAI_ASSERT_EQ_DEBUG( offsets.size(), size + 1, "Internal error: wrong offsets" )
    SCAI_ASSERT_EQ_ERROR( perm.size(), owners.size(), "Illegal size for permutation, most likely due to out-of-range index values" )

    HArray<IndexType> sendIndexes;  // local index values sorted by owner

    HArray<ValueType> sendValues;   // local source values sorted same as index values

    HArrayUtils::gather( sendIndexes, index.getLocalValues(), perm, utilskernel::binary::COPY );

    HArrayUtils::gather( sendValues, source.getLocalValues(), perm, utilskernel::binary::COPY );

    // exchange communication plans

    dmemo::CommunicationPlan sendPlan;   // will be allocated by offsets from bucket sort

    dmemo::CommunicationPlan recvPlan;   // is the transposed send plan

    {
        hmemo::ReadAccess<IndexType> rOffsets( offsets );
        sendPlan.allocateByOffsets( rOffsets.get(), size );
    }

    recvPlan.allocateTranspose( sendPlan, comm );

    SCAI_LOG_DEBUG( logger, comm << ": sendPlan = " << sendPlan << ", recvPlan = " << recvPlan )

    HArray<IndexType> recvIndexes;
    HArray<ValueType> recvValues;

    comm.exchangeByPlan( recvIndexes, recvPlan, sendIndexes, sendPlan );
    comm.exchangeByPlan( recvValues, recvPlan, sendValues, sendPlan );

    // translate global recvIndexes to local indexes, all must be local

    {
        WriteAccess<IndexType> wRecvIndexes( recvIndexes );

        for ( IndexType i = 0; i < recvIndexes.size(); ++i )
        {
            IndexType localIndex = targetDistribution.global2local( wRecvIndexes[i] );
            SCAI_ASSERT_NE_DEBUG( localIndex, nIndex, "got required index " << wRecvIndexes[i] << " but I'm not owner" )
            wRecvIndexes[i] = localIndex;
        }
    }

    // Now scatter all received values

    HArrayUtils::scatter( mLocalValues, recvIndexes, hasUniqueIndexes, recvValues, op, source.getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::dotProduct( const Vector& other ) const
{
    SCAI_REGION( "Vector.Dense.dotP" )
    SCAI_LOG_INFO( logger, "Calculating dot product for " << *this << " * " << other )

    // add other->getVectorKind() == DENSE, if sparse is also supported

    SCAI_ASSERT_EQ_ERROR( getValueType(), other.getValueType(),
                          "dotProduct not supported for different value types. "
                          << *this << " x " << other )

    SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(),
                          "dotProduct not supported for vectors with different distributions. "
                          << *this  << " x " << other )

    const DenseVector<ValueType>* denseOther = dynamic_cast<const DenseVector<ValueType>*>( &other );

    SCAI_ASSERT_ERROR( denseOther, "dynamic_cast failed for other = " << other )

    SCAI_LOG_DEBUG( logger, "Calculating local dot product at " << *mContext )
    const IndexType localSize = mLocalValues.size();
    SCAI_ASSERT_EQ_DEBUG( localSize, getDistribution().getLocalSize(), "size mismatch" )
    const ValueType localDotProduct = mLocalValues.dotProduct( denseOther->mLocalValues );
    SCAI_LOG_DEBUG( logger, "Calculating global dot product form local dot product = " << localDotProduct )
    ValueType dotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    SCAI_LOG_DEBUG( logger, "Global dot product = " << dotProduct )

    return Scalar( dotProduct );
}

template<typename ValueType>
void DenseVector<ValueType>::scale( const Vector& other )
{
    SCAI_REGION( "Vector.Dense.scale" )
    SCAI_LOG_INFO( logger, "Scale " << *this << " with " << other )

    if ( getDistribution() != other.getDistribution() )
    {
        COMMON_THROWEXCEPTION( "distribution do not match for this * other, this = " << *this << " , other = " << other )
    }

    other.buildLocalValues( mLocalValues, utilskernel::binary::MULT, getContextPtr() );
}

template<typename ValueType>
void DenseVector<ValueType>::allocate()
{
    // allocate local dense values with correct local size, touch on its context
    WriteOnlyAccess<ValueType> dummyWAccess( mLocalValues, mContext, getDistribution().getLocalSize() );
}

template<typename ValueType>
void DenseVector<ValueType>::allocate( DistributionPtr distribution )
{
    setDistributionPtr( distribution );
    allocate();
}

template<typename ValueType>
void DenseVector<ValueType>::allocate( const IndexType n )
{
    setDistributionPtr( DistributionPtr( new NoDistribution( n ) ) );
    allocate();
}

template<typename ValueType>
void DenseVector<ValueType>::assign( const Scalar value )
{
    SCAI_LOG_DEBUG( logger, *this << ": assign " << value )
    // assign the scalar value on the home of this dense vector.
    HArrayUtils::setScalar( mLocalValues, value.getValue<ValueType>(), utilskernel::binary::COPY, mContext );
}

template<typename ValueType>
void DenseVector<ValueType>::add( const Scalar value )
{
    SCAI_LOG_DEBUG( logger, *this << ": add " << value )
    // assign the scalar value on the home of this dense vector.
    HArrayUtils::setScalar( mLocalValues, value.getValue<ValueType>(), utilskernel::binary::ADD, mContext );
}

template<typename ValueType>
void DenseVector<ValueType>::buildLocalValues(

     _HArray& localValues, 
     const utilskernel::binary::BinaryOp op,
     ContextPtr prefLoc ) const

{
    if ( op == utilskernel::binary::COPY )
    {
        HArrayUtils::assign( localValues, mLocalValues, prefLoc );
    }
    else
    {
        SCAI_ASSERT_EQ_ERROR( localValues.size(), mLocalValues.size(), "size mismatch" )
        HArrayUtils::setArray( localValues, mLocalValues, op, prefLoc );
    }
}

template<typename ValueType>
void DenseVector<ValueType>::prefetch( const ContextPtr location ) const
{
    mLocalValues.prefetch( location );
}

template<typename ValueType>
void DenseVector<ValueType>::wait() const
{
    mLocalValues.wait();
}

template<typename ValueType>
void DenseVector<ValueType>::invert()
{
    mLocalValues.invert();
}

template<typename ValueType>
void DenseVector<ValueType>::conj()
{
    mLocalValues.conj();
}

template<typename ValueType>
void DenseVector<ValueType>::exp()
{
    mLocalValues.exp();
}

template<typename ValueType>
void DenseVector<ValueType>::log()
{
    mLocalValues.log();
}

template<typename ValueType>
void DenseVector<ValueType>::floor()
{
    mLocalValues.floor();
}

template<typename ValueType>
void DenseVector<ValueType>::ceil()
{
    mLocalValues.ceil();
}

template<typename ValueType>
void DenseVector<ValueType>::sqrt()
{
    mLocalValues.sqrt();
}

template<typename ValueType>
void DenseVector<ValueType>::sin()
{
    mLocalValues.sin();
}

template<typename ValueType>
void DenseVector<ValueType>::cos()
{
    mLocalValues.cos();
}

template<typename ValueType>
void DenseVector<ValueType>::tan()
{
    mLocalValues.tan();
}

template<typename ValueType>
void DenseVector<ValueType>::atan()
{
    mLocalValues.atan();
}

template<typename ValueType>
void DenseVector<ValueType>::powExp( const Vector& other )
{
    const DenseVector<ValueType>& denseOther = dynamic_cast<const DenseVector<ValueType>&>( other );
    mLocalValues.powExp( denseOther.mLocalValues );
}

template<typename ValueType>
void DenseVector<ValueType>::powBase( const Vector& other )
{
    const DenseVector<ValueType>& denseOther = dynamic_cast<const DenseVector<ValueType>&>( other );
    mLocalValues.powBase( denseOther.mLocalValues );
}

template<typename ValueType>
void DenseVector<ValueType>::powBase( Scalar base )
{
    mLocalValues.powBase( base.getValue<ValueType>() );
}

template<typename ValueType>
void DenseVector<ValueType>::powExp( Scalar exp )
{
    mLocalValues.powExp( exp.getValue<ValueType>() );
}

template<typename ValueType>
size_t DenseVector<ValueType>::getMemoryUsage() const
{
    // Note: memory of mHaloValues is not counted, is just a temporary

    IndexType localSize = mLocalValues.size();

    // Note: do sum with IndexType, as size_t is not yet handled by TypeTraits

    IndexType globalSize = getDistribution().getCommunicator().sum( localSize );

    return sizeof( ValueType ) * globalSize;
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseVector<ValueType>::redistribute( DistributionPtr distribution )
{
    SCAI_LOG_INFO( logger, "redistribute " << *this << ", new dist = " << *distribution )

    SCAI_ASSERT_EQ_ERROR( size(), distribution->getGlobalSize(), "global size mismatch between old/new distribution" )

    if ( getDistribution() == *distribution )
    {
        SCAI_LOG_DEBUG( logger, *this << " redistribute to same distribution " << *distribution )
        // we can keep local/global values, but just set dist pointer
        setDistributionPtr( distribution );
    }
    else if ( getDistribution().isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": replicated vector" << " will be localized to " << *distribution )
        HArray<ValueType> newLocalValues;
        ContextPtr hostContext = Context::getHostPtr();
        {
            const IndexType newSize = distribution->getLocalSize();
            ReadAccess<ValueType> rLocalValues( mLocalValues, hostContext );
            WriteOnlyAccess<ValueType> wNewLocalValues( newLocalValues, hostContext, newSize );
            #pragma omp parallel for

            for ( IndexType i = 0; i < size(); ++i )
            {
                if ( distribution->isLocal( i ) )
                {
                    const IndexType iLocal = distribution->global2local( i );
                    SCAI_ASSERT_DEBUG( iLocal < newSize, "illegal index " << iLocal )
                    wNewLocalValues[iLocal] = rLocalValues[i];
                }
            }
        }
        mLocalValues.swap( newLocalValues );
        setDistributionPtr( distribution );
    }
    else if ( distribution->isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, *this << " will be replicated" )
        // replicate a distributed vector
        HArray<ValueType> globalValues;
        SCAI_LOG_DEBUG( logger, "globalValues" )
        ContextPtr hostContext = Context::getHostPtr();
        {
            ReadAccess<ValueType> localData( mLocalValues, hostContext );
            WriteOnlyAccess<ValueType> globalData( globalValues, hostContext, size() );
            getDistribution().replicate( globalData.get(), localData.get() );
        }

        mLocalValues.swap( globalValues );
        setDistributionPtr( distribution );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, *this << " will be redistributed to " << *distribution )
        // so we have now really a redistibution, build a Redistributor
        HArray<ValueType> newLocalValues( distribution->getLocalSize() );
        Redistributor redistributor( distribution, getDistributionPtr() ); // target, source distributions
        redistributor.redistribute( newLocalValues, mLocalValues );
        mLocalValues.swap( newLocalValues );
        setDistributionPtr( distribution );
    }
}

/* -- IO ------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::writeLocalToFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::scalar::ScalarType dataType,
    const FileIO::FileMode fileMode
) const
{
    std::string suffix = fileType;

    if ( suffix == "" )
    {
        suffix = FileIO::getSuffix( fileName );
    }

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        common::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        if ( dataType != common::scalar::UNKNOWN )
        {
            // overwrite the default settings

            fileIO->setDataType( dataType );
        }

        if ( fileMode != FileIO::DEFAULT_MODE )
        {
            // overwrite the default settings

            fileIO->setMode( fileMode );
        }

        fileIO->writeArray( mLocalValues, fileName );
    }
    else
    {
        COMMON_THROWEXCEPTION( "File : " << fileName << ", unknown suffix" )
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType DenseVector<ValueType>::readLocalFromFile( const std::string& fileName, const IndexType first, const IndexType n )
{
    SCAI_LOG_INFO( logger, "read local array from file " << fileName )

    FileIO::read( mLocalValues, fileName, common::scalar::INTERNAL, first, n );

    return mLocalValues.size();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void DenseVector<ValueType>::clearValues()
{
    mLocalValues.clear();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Vector* DenseVector<ValueType>::create()
{
    return new DenseVector<ValueType>();
}

template<typename ValueType>
VectorCreateKeyType DenseVector<ValueType>::createValue()
{
    return VectorCreateKeyType( Vector::DENSE, common::getScalarType<ValueType>() );
}

template<typename ValueType>
VectorCreateKeyType DenseVector<ValueType>::getCreateValue() const
{
    return createValue();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const DenseVector<ValueType>& other ) :

    _DenseVector( other )

{
    // implementation here can be simpler as DenseVector( const Vector& other )
    SCAI_LOG_INFO( logger,
                   "Copy of vector of global size " << size() << ", local size " << getDistribution().getLocalSize() )
    mLocalValues = other.getLocalValues();
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DenseVector, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
