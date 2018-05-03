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
#include <scai/lama/SparseVector.hpp>

// local library
#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/expression/Expression.hpp>

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/PartitionIO.hpp>

#include <scai/lama/mepr/VectorWrapper.hpp>

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>
#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/unsupported.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/mepr/TypeList.hpp>

// std
#include <ostream>
#include <memory>

namespace scai
{

using common::TypeTraits;
using common::BinaryOp;
using utilskernel::HArrayUtils;

using namespace hmemo;
using namespace dmemo;

namespace lama
{

/* ------------------------------------------------------------------------- */
/*  Constructors of DenseVector<ValueType>                                   */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( ContextPtr context ) :

    Vector<ValueType>( 0, context ),
    mLocalValues()
{
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const IndexType size, const ValueType value, ContextPtr context ) :

    Vector<ValueType>( size, context ), 
    mLocalValues( size, value, context )

{
    SCAI_LOG_INFO( logger, "Construct dense vector, size = " << size << ", init =" << value )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, const ValueType value, ContextPtr context ) :

    Vector<ValueType>( distribution, context ), 
    mLocalValues( distribution->getLocalSize(), value )

{
    SCAI_LOG_INFO( logger,
                   "Construct dense vector, size = " << distribution->getGlobalSize() 
                    << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", value = " << value )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::splitUp( DistributionPtr& dist, hmemo::HArray<ValueType>& localValues )
{
    dist        = getDistributionPtr();
    localValues = std::move( mLocalValues );

    setDistributionPtr( std::make_shared<NoDistribution>( 0 ) );
}

/* ------------------------------------------------------------------------- */

template <typename ValueType>
void DenseVector<ValueType>::fillLinearValues( const ValueType startValue, const ValueType inc )
{
    const Distribution& dist = getDistribution();

    if ( dist.isReplicated() )
    {
        SCAI_ASSERT_EQ_DEBUG( dist.getGlobalSize(), dist.getLocalSize(), dist << " not replicated" );

        HArrayUtils::setSequence( mLocalValues, startValue, inc, dist.getGlobalSize() );
        return;
    }

    ContextPtr context = getContextPtr();

    // get my owned indexes

    HArray<IndexType> myGlobalIndexes( context );

    // mult with inc and add startValue

    dist.getOwnedIndexes( myGlobalIndexes );

    // localValues[] =  indexes[] * inc + startValue

    HArrayUtils::assign( mLocalValues, myGlobalIndexes, context );
    HArrayUtils::setScalar( mLocalValues, inc, BinaryOp::MULT, context );
    HArrayUtils::setScalar( mLocalValues, startValue, BinaryOp::ADD, context );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::fillRandom( const IndexType bound )
{
    HArrayUtils::fillRandom( mLocalValues, bound, 1.0f, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::fillSparseRandom( const float fillRate, const IndexType bound )
{
    HArrayUtils::fillRandom( mLocalValues, bound, fillRate, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, HArray<ValueType> localValues, ContextPtr context ) :

    Vector<ValueType>( distribution, context ),
    mLocalValues( std::move( localValues ) )
{
    SCAI_ASSERT_EQ_ERROR( mLocalValues.size(), distribution->getLocalSize(), "size mismatch" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( HArray<ValueType> values, ContextPtr context ) :

    Vector<ValueType>( std::make_shared<NoDistribution>( values.size() ), context ),
    mLocalValues( std::move( values ) )

{
    SCAI_LOG_DEBUG( logger, "DenseVector<" << common::TypeTraits<ValueType>::id()
                             << ">( moved HArray<" << values.getValueType() 
                             << ">, size = " << values.size() << " )" )
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

    setDistributionPtr( other.getDistributionPtr() );
    mLocalValues = other.mLocalValues;

    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>& DenseVector<ValueType>::operator=( DenseVector<ValueType>&& other ) noexcept
{
    setDistributionPtr( other.getDistributionPtr() );
    mLocalValues = std::move( other.mLocalValues );
    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool DenseVector<ValueType>::isConsistent() const
{
    // for a dense vector we just have to check that the locally allocated
    // data has the same size as the local size of the distribution

    const Distribution& dist = getDistribution();

    IndexType consistencyErrors = 0;

    const IndexType localSize = dist.getLocalSize();

    if ( getLocalValues().size() != localSize )
    {
        consistencyErrors++;
    }

    // use communicator for global reduction to make sure that all processors return same value.

    consistencyErrors = dist.getCommunicator().sum( consistencyErrors );

    return 0 == consistencyErrors;
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
    const PartitionId numPartitions = comm.getSize();

    if ( ascending )
    {
        ValueType minV = n > 0 ? sortedValues[0] : TypeTraits<ValueType>::getMax();
        ValueType maxV = n > 0 ? sortedValues[n - 1] : TypeTraits<ValueType>::getMin();

        splitValues[0]           = comm.min( minV );
        splitValues[numPartitions] = comm.max( maxV );
    }
    else
    {
        ValueType maxV = n > 0 ? sortedValues[0] : TypeTraits<ValueType>::getMin();
        ValueType minV = n > 0 ? sortedValues[n - 1] : TypeTraits<ValueType>::getMax();

        splitValues[0]           = comm.max( maxV );
        splitValues[numPartitions] = comm.min( minV );
    }

    // fill intermediate values by uniform distribution of range splitValues[0] .. splitValues[numPartitions]

    for ( PartitionId p = 1; p < numPartitions; ++p )
    {
        splitValues[p] = splitValues[0] + ( splitValues[numPartitions] - splitValues[0] ) * ValueType( p ) / ValueType( numPartitions );
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

    SCAI_ASSERT_NE_ERROR( blockSize, invalidIndex, "isSorted only possible on block distributed vectors" )

    const HArray<ValueType>& localValues = getLocalValues();

    bool is = HArrayUtils::isSorted( localValues, ascending ? common::CompareOp::LE : common::CompareOp::GE );

    if ( size == 1 )
    {
        return is;    // we are done
    }

    SCAI_LOG_DEBUG( logger, "local array is sorted = " << is )

    const IndexType n = localValues.size();

    ValueType tmp;

    const int toRight = 1;
    const int toLeft  = -1;

    if ( n == 0 )
    {
        comm.shift( &tmp, IndexType( 1 ), &tmp, IndexType( 0 ), toRight );
        comm.shift( &tmp, IndexType( 1 ), &tmp, IndexType( 0 ), toLeft );

        is = comm.all( is );

        return is;
    }

    ReadAccess<ValueType> rLocalValues( localValues );

    // send right value to right neighbor

    IndexType nt = comm.shift( &tmp, IndexType( 1 ), &rLocalValues[n - 1], IndexType( 1 ), toRight );

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

    nt = comm.shift( &tmp, IndexType( 1 ), &rLocalValues[0], IndexType( 1 ), toLeft );

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

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
bool DenseVector<ComplexFloat>::isSorted( bool ) const
{
    COMMON_THROWEXCEPTION( "isSorted unsupported for complex vectors." )
}

template<>
bool DenseVector<ComplexDouble>::isSorted( bool ) const
{
    COMMON_THROWEXCEPTION( "isSorted unsupported for complex vectors." )
}

template<>
bool DenseVector<ComplexLongDouble>::isSorted( bool ) const
{
    COMMON_THROWEXCEPTION( "isSorted unsupported for complex vectors." )
}

#endif

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

    SCAI_ASSERT_NE_ERROR( blockSize, invalidIndex, "sort only possible on block distributed vectors" )

    // Now sort the local values

    HArray<IndexType>* localPerm = NULL;

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

            const IndexType globalLow = distribution.local2global( 0 );
            HArrayUtils::setScalar( *localPerm, globalLow, common::BinaryOp::ADD );
        }

        ReadAccess<IndexType> rPerm( *localPerm );

        for ( IndexType i = 0; i < nLocalPerm; ++i )
        {
            SCAI_LOG_TRACE( logger, "localPerm[" << i << "] = " << rPerm[i] )
        }
    }

    // Determine the splitting values

    PartitionId numPartitions = comm.getSize();

    std::unique_ptr<ValueType[]> splitValues( new ValueType[ numPartitions + 1 ] );

    {
        ReadAccess<ValueType> rSortedValues( sortedValues );
        getSplitValues( splitValues.get(), comm, rSortedValues.get(), sortedValues.size(), ascending );
    }

    std::unique_ptr<IndexType[]> quantities( new IndexType[ numPartitions ] );

    // Determine quantities for each processor

    for ( PartitionId ip = 0; ip < numPartitions; ++ip )
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
                    SCAI_ASSERT_LT_DEBUG( p, numPartitions, "Illegal split values, mLocalValues[" << i << "] = " << rValues[i] );
                }
            }
            else
            {
                while ( rValues[i] < splitValues[p + 1] )
                {
                    p++;
                    SCAI_ASSERT_LT_DEBUG( p, numPartitions, "Illegal split values, mLocalValues[" << i << "] = " << rValues[i] );
                }
            }

            quantities[p]++;
        }
    }

    // make communication plans for sending data and receiving data

    auto sendPlan = dmemo::CommunicationPlan::buildBySizes( quantities.get(), numPartitions );
    auto recvPlan = sendPlan.transpose( comm );

    SCAI_LOG_INFO( logger, comm << ": send plan: " << sendPlan << ", rev plan: " << recvPlan );

    HArray<ValueType> newValues;

    IndexType newLocalSize = recvPlan.totalQuantity();

    {
        WriteOnlyAccess<ValueType> recvVals( newValues, newLocalSize );
        ReadAccess<ValueType> sendVals( sortedValues );
        comm.exchangeByPlan( recvVals.get(), recvPlan, sendVals.get(), sendPlan );
    }

    // Also communicate the original index positions of the array values

    if ( perm )
    {
        HArray<IndexType> newPerm;

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

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
void DenseVector<ComplexFloat>::sortImpl( DenseVector<IndexType>*, DenseVector<ComplexFloat>*, DenseVector<ComplexFloat>&, bool )
{
    COMMON_THROWEXCEPTION( "no sort on complex vector" )
}

template<>
void DenseVector<ComplexDouble>::sortImpl( DenseVector<IndexType>*, DenseVector<ComplexDouble>*, DenseVector<ComplexDouble>&, bool )
{
    COMMON_THROWEXCEPTION( "no sort on complex vector" )
}

template<>
void DenseVector<ComplexLongDouble>::sortImpl( DenseVector<IndexType>*, DenseVector<ComplexLongDouble>*, DenseVector<ComplexLongDouble>&, bool )
{
    COMMON_THROWEXCEPTION( "no sort on complex vector" )
}

#endif

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::scan()
{
	// first, check that the input is some block distribution

    SCAI_ASSERT_NE_ERROR( getDistribution().getBlockDistributionSize(), invalidIndex,
                          "scan only supported for block distribution" )

    const Communicator& comm = getDistribution().getCommunicator();
 
    HArray<ValueType> prefixValues;
    
    ValueType val = HArrayUtils::sum( mLocalValues );

    ValueType scanVal = comm.scan( val ); // is inclusve scan

    scanVal -= val;  // exclusive value is needed

    SCAI_LOG_INFO( logger, comm << ": local sum = " << val << ", scan = " << scanVal )

    // now do the correct scan, start with val

    bool exclusive = false;

    HArrayUtils::scan( mLocalValues, scanVal, exclusive, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::scan( const DenseVector<ValueType>& other )
{
    // first, check that the input is some block distribution

    SCAI_ASSERT_NE_ERROR( other.getDistribution().getBlockDistributionSize(), invalidIndex,
                          "scan only supported for block distribution" )

    // currently scan is only supported by in-place array operations

    assign( other );
    scan();
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
void DenseVector<ValueType>::setDenseValues( const _HArray& values )
{
    const IndexType localSize = getDistribution().getLocalSize();

    SCAI_ASSERT_EQ_ERROR( localSize, values.size(), "size of array with local values does not match local size of distribution" )

    HArrayUtils::_assign( mLocalValues, values, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::fillSparseData( 
    const HArray<IndexType>& nonZeroIndexes, 
    const _HArray& nonZeroValues,
    const BinaryOp op )
{
    // Note: scatter checks for legal indexes

    HArrayUtils::_scatter( mLocalValues, nonZeroIndexes, false, nonZeroValues, op, getContextPtr() );
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
    std::unique_ptr<DenseVector<ValueType> > vector( new DenseVector<ValueType>( getContextPtr() ) );
    vector->setSameValue( this->getDistributionPtr(), 0 );
    return vector.release();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseVector<ValueType>::getValue( IndexType globalIndex ) const
{
    ValueType myValue = 0;

    const IndexType localIndex = getDistribution().global2local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": getValue( globalIndex = " << globalIndex << " ) -> local : " << localIndex )

    if ( localIndex != invalidIndex )
    {
        myValue = mLocalValues[localIndex];
    }

    ValueType allValue = getDistribution().getCommunicator().sum( myValue );

    // works also fine for replicated distributions with NoCommunicator

    SCAI_LOG_TRACE( logger, "myValue = " << myValue << ", allValue = " << allValue )

    return allValue;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::setValue( const IndexType globalIndex, const ValueType value )
{
    SCAI_ASSERT_VALID_INDEX_ERROR( globalIndex, size(), "out of range index" )

    SCAI_LOG_TRACE( logger, *this << ": setValue( globalIndex = " << globalIndex << " ) = " <<  value )

    const IndexType localIndex = getDistribution().global2local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": set @g " << globalIndex << " is @l " << localIndex << " : " << value )

    if ( localIndex != invalidIndex )
    {
        mLocalValues[localIndex] = value;
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseVector<ValueType>::min() const
{
    ValueType localMin = HArrayUtils::min( mLocalValues );
    return getDistribution().getCommunicator().min( localMin );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseVector<ValueType>::max() const
{
    // Note: max returns the minimal representation value on zero-sized vectors
    ValueType localMax = HArrayUtils::max( mLocalValues );
    return getDistribution().getCommunicator().max( localMax );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseVector<ValueType>::l1Norm() const
{
    auto localL1Norm = HArrayUtils::l1Norm( mLocalValues );
    return getDistribution().getCommunicator().sum( localL1Norm );
}

/*---------------------------------------------------------------------------*/
template<typename ValueType>
ValueType DenseVector<ValueType>::sum() const
{
    auto localsum = HArrayUtils::sum( mLocalValues );
    return getDistribution().getCommunicator().sum( localsum );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseVector<ValueType>::l2Norm() const
{
    // Note: we do not call l2Norm here for mLocalValues to avoid sqrt
    RealType<ValueType> localDotProduct = HArrayUtils::dotProduct( mLocalValues, mLocalValues );
    RealType<ValueType> globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    return common::Math::sqrt( globalDotProduct );
}

template<>
IndexType DenseVector<IndexType>::l2Norm() const
{
    // Note: we do not call l2Norm here for mLocalValues to avoid sqrt

    double localDotProduct = static_cast<double>( HArrayUtils::dotProduct( mLocalValues, mLocalValues ) );
    double globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    return IndexType( common::Math::sqrt( globalDotProduct ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseVector<ValueType>::maxNorm() const
{
    RealType<ValueType> localMaxNorm = HArrayUtils::maxNorm( mLocalValues );
    const Communicator& comm = getDistribution().getCommunicator();
    RealType<ValueType> globalMaxNorm = comm.max( localMaxNorm );
    SCAI_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm: " << localMaxNorm
                   << ", max norm global = " << globalMaxNorm )
    return globalMaxNorm;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseVector<ValueType>::maxDiffNorm( const Vector<ValueType>& other ) const
{
    bool temporaryNeeded = false;

    if ( other.getDistribution() != getDistribution() )
    {
        SCAI_UNSUPPORTED( "maxDiffNorm( x ), distribution mismatch implies temporary" )
        temporaryNeeded = true;
    }

    if ( other.getVectorKind() != VectorKind::DENSE )
    {
        // ToDo: should be supported maxDiff for SPARSE-DENSE 
        SCAI_UNSUPPORTED( "maxDiffNorm( x ), x not DENSE, uses temporary" )
        temporaryNeeded = true;
    }

    if ( temporaryNeeded )
    {
        // convert/distribute other vector so that it fits
        return maxDiffNorm( distribute<DenseVector<ValueType>>( other, getDistributionPtr() ) );
    }

    SCAI_ASSERT_DEBUG( dynamic_cast<const DenseVector<ValueType>*>( &other ), "wrong cast" )

    const DenseVector<ValueType>& denseOther = static_cast<const DenseVector<ValueType>&>( other );

    RealType<ValueType> localMaxNorm = HArrayUtils::maxDiffNorm( mLocalValues, denseOther.getLocalValues() );

    const Communicator& comm = getDistribution().getCommunicator();

    RealType<ValueType> globalMaxNorm = comm.max( localMaxNorm );

    SCAI_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm: " << localMaxNorm
                   << ", max norm global = " << globalMaxNorm )

    return globalMaxNorm;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool DenseVector<ValueType>::all( const common::CompareOp op, const ValueType value ) const
{
    bool localAll = HArrayUtils::allScalar( getLocalValues(), op, value );

    bool globalAll = getDistribution().getCommunicator().all( localAll );

    return globalAll;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool DenseVector<ValueType>::all( const common::CompareOp op, const Vector<ValueType>& other ) const
{
    SCAI_ASSERT_EQ_ERROR( other.getDistribution(), getDistribution(), "distribution mismatch for all compare, op = " << op )

    bool localAll;

    if ( other.getVectorKind() == VectorKind::DENSE )
    {
        const DenseVector<ValueType>& denseOther = static_cast<const DenseVector<ValueType>&>( other );
        localAll = HArrayUtils::all( getLocalValues(), op, denseOther.getLocalValues() );
    }
    else
    {
        HArray<ValueType> otherLocalValues;
        other.buildLocalValues( otherLocalValues );
        localAll = HArrayUtils::all( getLocalValues(), op, otherLocalValues );
    }

    bool globalAll = getDistribution().getCommunicator().all( localAll );

    return globalAll;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::swap( _Vector& other )
{
    SCAI_LOG_DEBUG( logger, "swap:" << *this << " with " << other )
    DenseVector* otherPtr = dynamic_cast<DenseVector*>( &other );

    if ( !otherPtr )
    {
        COMMON_THROWEXCEPTION( "Tried to swap with a Vector of a different type." )
    }

    _Vector::swapVector( other );
    mLocalValues.swap( otherPtr->mLocalValues );
    // mHaloValues.swap( otherPtr->mHaloValues );
}

template<typename ValueType>
void DenseVector<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "DenseVector<" << getValueType() << ">" << "( size = " << size() << ", local = " << mLocalValues.size()
           << ", dist = " << getDistribution() << ", loc  = " << *getContextPtr() << " )";
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::setScalar( const ValueType& alpha )
{
    SCAI_LOG_DEBUG( logger, *this << ": set scalar " << alpha )

    // assign the scalar alpha on the home of this dense vector.
    HArrayUtils::setScalar( mLocalValues, alpha, common::BinaryOp::COPY, this->getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::axpy( const ValueType& alpha, const Vector<ValueType>& x )
{
    SCAI_LOG_INFO( logger, "axpy: this += " << alpha << " * x, with x = " << x << ", this = " << *this )

    SCAI_ASSERT_EQ_ERROR( x.getDistribution(), getDistribution(), "distribution mismatch for axpy vector" )

    if ( x.getVectorKind() == VectorKind::SPARSE )
    {
        SCAI_REGION( "Vector.Dense.axpySparse" )

        SCAI_ASSERT_DEBUG( dynamic_cast<const SparseVector<ValueType>*>( &x ), "illegal cast" );

        const SparseVector<ValueType>& sparseX = static_cast<const SparseVector<ValueType>&>( x );

        const HArray<IndexType>& nonZeroIndexes = sparseX.getNonZeroIndexes();
        const HArray<ValueType>& nonZeroValues  = sparseX.getNonZeroValues();

        const ValueType xZero = sparseX.getZero();

        const bool unique = true;  // non-zero indexes in sparse vectors are always unique

        if ( xZero != common::Constants::ZERO )
        {
            // we have also to add the zero values

            HArray<ValueType> xDenseValues;
            sparseX.buildLocalValues( xDenseValues);
            utilskernel::HArrayUtils::axpy( mLocalValues, alpha, xDenseValues, getContextPtr() );
        }
        else if ( alpha == common::Constants::ONE )
        {
            HArrayUtils::scatter( mLocalValues, nonZeroIndexes, unique, nonZeroValues, BinaryOp::ADD, getContextPtr());
        }
        else
        {
            // ToDo: if there is scatter with scaling factor for the scattered values, employ it here

            HArray<ValueType> tmpValues;
            HArrayUtils::compute( tmpValues, nonZeroValues, BinaryOp::MULT, alpha, getContextPtr() );
            HArrayUtils::scatter( mLocalValues, nonZeroIndexes, unique, tmpValues, BinaryOp::ADD, getContextPtr() );
        }
    }
    else
    {
        SCAI_REGION( "Vector.Dense.axpyDense" )

        const DenseVector<ValueType>& denseX = static_cast<const DenseVector<ValueType>&>( x );

        const HArray<ValueType>& xValues = denseX.getLocalValues();

        utilskernel::HArrayUtils::axpy( mLocalValues, alpha, xValues, getContextPtr() );
    }
}

/* ------------------------------------------------------------------------- */
/*   this = alpha * x + beta * y                                             */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::vectorPlusVector( const ValueType& alpha, const Vector<ValueType>& x, 
                                               const ValueType& beta, const Vector<ValueType>& y )
{
    SCAI_LOG_INFO( logger, "vectorPlusVector: z = " << alpha << " * x + " << beta << " * y"
                             << ", with  x = " << x << ", y = " << y << ", z = " << *this )

    // Query for alpha == 0 or beta == 0 as x and y might be undefined 

    if ( alpha == ValueType( 0 ) )
    {
        // vector x is completely ignored, there are even no space checks
        binaryOpScalar( y, beta, BinaryOp::MULT, true );
        return;
    }

    if ( beta == ValueType( 0 ) )
    {
        // vector y is completely ignored, there are even no space checks
        binaryOpScalar( x, alpha, BinaryOp::MULT, true );
        return;
    }

    if ( alpha == ValueType( 1 ) && ( &x == this ) )
    {
        // take advantage of alias
        axpy( beta, y );
        return;
    }

    if ( beta == ValueType( 1 ) && ( &y == this ) )
    {
        // take advantage of alias
        axpy( alpha, x );
        return;
    }

    SCAI_ASSERT_EQ_ERROR( x.getDistribution(), y.getDistribution(),
                          "size/distribution mismatch of operands in alpha * x + beta * y" )

    if ( x.getVectorKind() == VectorKind::SPARSE && y.getVectorKind() == VectorKind::SPARSE )
    {
        // use a temporary sparse vector for addition 

        SCAI_UNSUPPORTED( "vectorPlusVector: use temporary SparseVector<" << getValueType() << ">" )

        SparseVector<ValueType> sTmp;
        sTmp.vectorPlusVector( alpha, x, beta, y );
        this->assign( sTmp );
    }

    if ( x.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_UNSUPPORTED( "vectorPlusVector: use temporary DenseVector<" << getValueType() << "> for x = " << x )
        vectorPlusVector( alpha, convert<DenseVector<ValueType>>( x ), beta, y );
        return;
    }

    if ( y.getVectorKind() != VectorKind::DENSE  )
    {
        SCAI_UNSUPPORTED( "vectorPlusVector: use temporary DenseVector<" << getValueType() << "> for y = " << y )
        vectorPlusVector( alpha, x, beta, convert<DenseVector<ValueType>>( y ) );
        return;
    }

    SCAI_REGION( "Vector.Dense.VplusV" )

    if ( x.getDistribution() != getDistribution() )
    {
        allocate( x.getDistributionPtr() );   // result inherits distribution of operands
    }

    const DenseVector<ValueType>& denseX = static_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>& denseY = static_cast<const DenseVector<ValueType>&>( y );

    SCAI_ASSERT_EQ_DEBUG( mLocalValues.size(), denseX.mLocalValues.size(), "serious space mismatch" )
    SCAI_ASSERT_EQ_DEBUG( mLocalValues.size(), denseY.mLocalValues.size(), "serious space mismatch" )

    SCAI_LOG_DEBUG( logger, "call arrayPlusArray" )

    utilskernel::HArrayUtils::arrayPlusArray( mLocalValues, alpha, denseX.mLocalValues, 
                                                            beta, denseY.mLocalValues, getContextPtr() );
}

/* ------------------------------------------------------------------------- */
/*   this = alpha * x * y   - elementwise multiplication                     */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::vectorTimesVector( 
    const ValueType& alpha, 
    const Vector<ValueType>& x, 
    const Vector<ValueType>& y )
{
    SCAI_LOG_INFO( logger, "z = x * y, z = " << *this << " , x = " << x << " , y = " << y )

    SCAI_LOG_DEBUG( logger, "dist of x = " << x.getDistribution() )
    SCAI_LOG_DEBUG( logger, "dist of y = " << y.getDistribution() )

    SCAI_ASSERT_EQ_ERROR( x.getDistribution(), y.getDistribution(),
                          "size/distribution mismatch of operands in alpha * x * y" )

    if ( x.getVectorKind() != VectorKind::DENSE || x.getValueType() != getValueType() )
    {
        vectorTimesVector( alpha, convert<DenseVector<ValueType>>( x ), y );
        return;
    }

    if ( y.getVectorKind() != VectorKind::DENSE || y.getValueType() != getValueType() )
    {
        vectorTimesVector( alpha, x, convert<DenseVector<ValueType>>( y ) );
        return;
    }

    setDistributionPtr( x.getDistributionPtr() );   // result inherits (same) space of operands, mostly its same

    const DenseVector<ValueType>& denseX = static_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>& denseY = static_cast<const DenseVector<ValueType>&>( y );

    // denseX.mLocalValues.size() == denseY.mLocalValues.size() is already verified
    // mLocalValues will be allocated with the correct size, alias are handled correctly.

    SCAI_LOG_DEBUG( logger, "call arrayTimesArray" )

    utilskernel::HArrayUtils::arrayTimesArray( mLocalValues, alpha, denseX.mLocalValues, denseY.mLocalValues, getContextPtr() );
}

/* ------------------------------------------------------------------------- */
/*   this = alpha * x + beta                                                 */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::vectorPlusScalar( const ValueType& alpha, const Vector<ValueType>& x, const ValueType& beta )
{
    if ( x.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_LOG_WARN( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">::vectorAddScalar, uses tmp for x" )

        vectorPlusScalar( alpha, convert<DenseVector<ValueType>>( x ), beta );
        return;
    }

    const DenseVector<ValueType>& denseX = dynamic_cast<const DenseVector<ValueType>&>( x );

    SCAI_LOG_INFO( logger, "z = alpha * x + beta, z = " << *this << ", alpha=  " << alpha
                   << " , x = " << x << " , beta = " << beta )
    SCAI_LOG_DEBUG( logger, "dist of x = " << x.getDistribution() )

    setDistributionPtr( x.getDistributionPtr() );

    SCAI_LOG_DEBUG( logger, "call arrayPlusScalar" )

    utilskernel::HArrayUtils::arrayPlusScalar( mLocalValues, alpha, denseX.mLocalValues, beta, getContextPtr() );

    SCAI_ASSERT_EQ_DEBUG( mLocalValues.size(), denseX.mLocalValues.size(), "serious mismatch" )
}

/* ----------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::selectComplexPart( Vector<RealType<ValueType> >& x, const common::ComplexPart part ) const
{
    HArray<RealType<ValueType> > localX;

    utilskernel::HArrayUtils::selectComplexPart( localX, getLocalValues(), part );

    x.assign( localX, getDistributionPtr() );
}

template<>
void DenseVector<IndexType>::selectComplexPart( Vector<IndexType>& x, const common::ComplexPart kind ) const
{
    if ( kind == common::ComplexPart::REAL )
    {
        x = *this;
    }
    else
    {
        x.setSameValue( getDistributionPtr(), IndexType( 0 ) );
    }
}

/* ----------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::buildComplex( const Vector<RealType<ValueType> >& x, const Vector<RealType<ValueType> >& y )
{
    SCAI_ASSERT_EQ_ERROR( x.getDistribution(), y.getDistribution(), "vectors must have same distribution" ) 

    typedef RealType<ValueType> real;

    if ( x.getVectorKind() != VectorKind::DENSE )
    {
        buildComplex( convert<DenseVector<real>>( x ), y );
    }
    else if ( y.getVectorKind() != VectorKind::DENSE )
    {
        buildComplex( x, convert<DenseVector<real>>( y ) );
    }
    else
    {
        const DenseVector<real>& xD = static_cast<const DenseVector<real>&>( x );
        const DenseVector<real>& yD = static_cast<const DenseVector<real>&>( y );

        setDistributionPtr( x.getDistributionPtr() );
        
        utilskernel::HArrayUtils::buildComplex( mLocalValues, xD.getLocalValues(), yD.getLocalValues(), getContextPtr() );
    }
}

/* ------------------------------------------------------------------------- */

template<>
void DenseVector<IndexType>::buildComplex( const Vector<IndexType>& x, const Vector<IndexType>& )
{
    *this = x;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::gather(
    const DenseVector<ValueType>& source,
    const DenseVector<IndexType>& index,
    const BinaryOp op )
{
    if ( op != BinaryOp::COPY )
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

    HArrayUtils::gather( requiredIndexes, index.getLocalValues(), perm, BinaryOp::COPY );

    // exchange communication plans

    auto recvPlan = dmemo::CommunicationPlan::buildByOffsets( hostReadAccess( offsets ).get(), size );
    auto sendPlan = recvPlan.transpose( comm );

    SCAI_LOG_DEBUG( logger, comm << ": recvPlan = " << recvPlan << ", sendPlan = " << sendPlan )

    HArray<IndexType> sendIndexes;

    comm.exchangeByPlan( sendIndexes, sendPlan, requiredIndexes, recvPlan );

    // translate global sendIndexes to local indexes, all must be local

    {
        WriteAccess<IndexType> wSendIndexes( sendIndexes );

        for ( IndexType i = 0; i < sendIndexes.size(); ++i )
        {
            IndexType localIndex = sourceDistribution.global2local( wSendIndexes[i] );
            SCAI_ASSERT_NE_DEBUG( localIndex, invalidIndex, "got required index " << wSendIndexes[i] << " but I'm not owner" )
            wSendIndexes[i] = localIndex;
        }
    }

    // exchange communication plan

    HArray<ValueType> sendValues;  // values to send from my source values

    HArrayUtils::gather( sendValues, source.getLocalValues(), sendIndexes, BinaryOp::COPY );

    // send via communication plan

    HArray<ValueType> recvValues;

    comm.exchangeByPlan( recvValues, recvPlan, sendValues, sendPlan );

    if ( op == BinaryOp::COPY )
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
    const BinaryOp op )
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

    HArrayUtils::gather( sendIndexes, index.getLocalValues(), perm, BinaryOp::COPY );

    HArrayUtils::gather( sendValues, source.getLocalValues(), perm, BinaryOp::COPY );

    // exchange communication plans

    auto sendPlan = dmemo::CommunicationPlan::buildByOffsets( hostReadAccess( offsets ).get(), size );
    auto recvPlan = sendPlan.transpose( comm );

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
            SCAI_ASSERT_NE_DEBUG( localIndex, invalidIndex, "got required index " << wRecvIndexes[i] << " but I'm not owner" )
            wRecvIndexes[i] = localIndex;
        }
    }

    // Now scatter all received values

    HArrayUtils::scatter( mLocalValues, recvIndexes, hasUniqueIndexes, recvValues, op, source.getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseVector<ValueType>::dotProduct( const Vector<ValueType>& other ) const
{
    SCAI_REGION( "Vector.Dense.dotP" )
    SCAI_LOG_INFO( logger, "Calculating dot product for " << *this << " * " << other )

    // add other->getVectorKind() == VectorKind::DENSE, if sparse is also supported

    SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(),
                          "dotProduct not supported for vectors with different distributions. "
                          << *this  << " x " << other )

    ValueType localDotProduct;

    if ( other.getVectorKind() == VectorKind::DENSE )
    {
        SCAI_ASSERT_DEBUG( dynamic_cast<const DenseVector<ValueType>*>( &other ), "dynamic cast failed, other = " << other )

        const DenseVector<ValueType>& denseOther = static_cast<const DenseVector<ValueType>&>( other );

        localDotProduct = HArrayUtils::dotProduct( mLocalValues, denseOther.getLocalValues() );
    }
    else
    {
        SCAI_ASSERT_DEBUG( dynamic_cast<const SparseVector<ValueType>*>( &other ), "dynamic cast failed, other = " << other )

        const SparseVector<ValueType>& sparseOther = static_cast<const SparseVector<ValueType>&>( other );
    
        HArray<ValueType> myValues;   // build values at same position as sparse vector

        gatherLocalValues( myValues, sparseOther.getNonZeroIndexes(), BinaryOp::COPY, getContextPtr() );

        localDotProduct = HArrayUtils::dotProduct( myValues, sparseOther.getNonZeroValues() );
    }

    SCAI_LOG_DEBUG( logger, "Calculating global dot product form local dot product = " << localDotProduct )
    ValueType dotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    SCAI_LOG_DEBUG( logger, "Global dot product = " << dotProduct )

    return dotProduct;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::setVector( const _Vector& other, BinaryOp op, const bool swapArgs )
{
    SCAI_REGION( "Vector.Dense.setVector" )

    SCAI_LOG_INFO( logger, "set " << *this << " with " << other << ", op = " << op )

    if ( op == common::BinaryOp::COPY )
    {
        SCAI_ASSERT_ERROR( !swapArgs, "swapping for binary COPY operator not allowed" )

        if ( &other == this )
        {
            return;
        }

        allocate( other.getDistributionPtr() );
    } 
    else 
    {
        SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(), "space mismatch" )
    }

    if ( !swapArgs )
    {
        other.buildLocalValues( mLocalValues, op, getContextPtr() );
    }
    else
    {
        COMMON_THROWEXCEPTION( "swapArgs not supported yet" );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::allocate()
{
    // allocate local dense values with correct local size, touch on its context
    WriteOnlyAccess<ValueType> dummyWAccess( mLocalValues, getContextPtr(), getDistribution().getLocalSize() );
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

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::assign( const _Vector& other )
{
    // translate virtual call to specific template call via wrapper

    mepr::VectorWrapper<DenseVector, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( this, other );
}

template<typename ValueType>
template<typename OtherValueType>
void DenseVector<ValueType>::assignImpl( const Vector<OtherValueType>& other )
{
    if ( other.getVectorKind() == VectorKind::SPARSE )
    {
        assignSparse( static_cast<const SparseVector<OtherValueType>&>( other ) );
    }
    else if ( other.getVectorKind() == VectorKind::DENSE )
    {
        assignDense( static_cast<const DenseVector<OtherValueType>&>( other ) );
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported vector kind for assign to sparse vector" )
    }
}

template<typename ValueType>
template<typename OtherValueType>
void DenseVector<ValueType>::assignSparse( const SparseVector<OtherValueType>& other )
{
    allocate( other.getDistributionPtr() );
    setScalar( static_cast<ValueType>( other.getZero() ) );
    fillSparseData( other.getNonZeroIndexes(), other.getNonZeroValues(), BinaryOp::COPY );
}

template<typename ValueType>
template<typename OtherValueType>
void DenseVector<ValueType>::assignDense( const DenseVector<OtherValueType>& other )
{

    setDistributionPtr( other.getDistributionPtr() );
    setDenseValues( other.getLocalValues() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::buildLocalValues(

     _HArray& localValues, 
     const BinaryOp op,
     ContextPtr prefLoc ) const

{
    HArrayUtils::_setArray( localValues, mLocalValues, op, prefLoc );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::gatherLocalValues(

     _HArray& localValues,
     const HArray<IndexType>& indexes,
     const BinaryOp op,
     ContextPtr prefLoc ) const
{
    HArrayUtils::_gather( localValues, mLocalValues, indexes, op, prefLoc );
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

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseVector<ValueType>::unaryOp( const Vector<ValueType>& x, common::UnaryOp op )
{
    if ( x.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_UNSUPPORTED( "denseVector = " << op << "( scalarVector ): tmpDenseVector = scalarVector added" )
        DenseVector<ValueType> tmpX;
        tmpX.assign( x );            // convert it, same distribution
        unaryOp( tmpX, op );
        return;
    }

    setDistributionPtr( x.getDistributionPtr() );     // do not allocate as it might invalidate x in case of alias

    const DenseVector<ValueType>& denseX = static_cast<const DenseVector<ValueType>&>( x );

    SCAI_LOG_INFO( logger, "this = " << op << " ( denseX ), denseX = " << x 
                   << ", mLocalValues = " << mLocalValues << ", xLocalValues = " << denseX.mLocalValues )

    HArrayUtils::unaryOp( mLocalValues, denseX.mLocalValues, op, getContextPtr() );
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseVector<ValueType>::binaryOp( const Vector<ValueType>& x, BinaryOp op, const Vector<ValueType>& y  )
{
    SCAI_ASSERT_EQ_ERROR( x.getDistribution(), y.getDistribution(), "distribution mismatch of vectors in binary op " << op )

    // Note: the logic of this routine is nearly the same as vectorPlusVector

    SCAI_LOG_INFO( logger, "binaryOp: this = x " << op << " y, with x = " << x << ", y = " << y );

    if ( x.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_UNSUPPORTED( "denseR = sparseX " << op << " y, tmpDenseVector = sparseX added" )
        DenseVector<ValueType> tmpX;
        tmpX.assign( x );
        binaryOp( tmpX, op, y );
        return;
    }

    if ( y.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_UNSUPPORTED( "denseR = x " << op << " sparseY, tmpDenseVector = sparseY added" )
        DenseVector<ValueType> tmpY;
        tmpY.assign( y );
        binaryOp( x, op, tmpY );
        return;
    }

    setDistributionPtr( x.getDistributionPtr() );  

    const DenseVector<ValueType>& denseX = static_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>& denseY = static_cast<const DenseVector<ValueType>&>( y );

    HArrayUtils::binaryOp( mLocalValues, denseX.getLocalValues(), denseY.getLocalValues(), op, getContextPtr() );
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseVector<ValueType>::binaryOpScalar( const Vector<ValueType>& x, const ValueType& alpha, const BinaryOp op, const bool swap )
{
    // Note: the logic of this routine is nearly the same as vectorPlusVector

    if ( swap )
    {
        SCAI_LOG_INFO( logger, "binaryOpScalar: this = " << alpha << " " << op << " x, with x = " << x );
    }
    else
    {
        SCAI_LOG_INFO( logger, "binaryOpScalar: this = x " << op << " " << alpha << ", with x = " << x );
    }
    
    if ( x.getVectorKind() != VectorKind::DENSE )
    {
        DenseVector<ValueType> tmpX;
        tmpX.assign( x );
        binaryOpScalar( tmpX, alpha, op, swap );
        return;
    }

    if ( &x != this )
    {
        allocate( x.getDistributionPtr() ); 
    }
    
    const DenseVector<ValueType>& denseX = static_cast<const DenseVector<ValueType>&>( x );
    
    HArrayUtils::binaryOpScalar( mLocalValues, denseX.getLocalValues(), alpha, op, swap, getContextPtr() );
}

/* ------------------------------------------------------------------------ */

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

/** 
 *  @brief Help function that computes the local part of a global array
 */

template<typename ValueType>
static void localize( HArray<ValueType>& out, const HArray<ValueType>& in, const Distribution& dist, ContextPtr ctx )
{
    SCAI_ASSERT_EQ_ERROR( in.size(), dist.getGlobalSize(), "input array must have global size of dist" );

    if ( dist.isReplicated() )
    {
        HArrayUtils::assign( out, in, ctx );
    }
    else
    {
        HArray<IndexType> myGlobalIndexes;
        dist.getOwnedIndexes( myGlobalIndexes );
        HArrayUtils::gather( out, in, myGlobalIndexes, common::BinaryOp::COPY, ctx );
    }

    SCAI_ASSERT_EQ_DEBUG( out.size(), dist.getLocalSize(), "serious mismatch: output array must have local size of dist" );
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
        localize( newLocalValues, mLocalValues, *distribution, getContextPtr() );
        mLocalValues.swap( newLocalValues );
        setDistributionPtr( distribution );
    }
    else if ( distribution->isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, *this << " will be replicated, i.e. combine all local parts" )
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

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseVector<ValueType>::redistribute( const Redistributor& redistributor )
{
    SCAI_ASSERT_EQ_ERROR( getDistribution(), *redistributor.getSourceDistributionPtr(), 
                          "redistributor does not match to actual distribution of this vector" );

    HArray<ValueType> newLocalValues( redistributor.getTargetLocalSize() );
    redistributor.redistribute( newLocalValues, mLocalValues );
    mLocalValues.swap( newLocalValues );
    setDistributionPtr( redistributor.getTargetDistributionPtr() );
}

/* -- IO ------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::writeLocalToFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::ScalarType dataType,
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

        std::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        if ( dataType != common::ScalarType::UNKNOWN )
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

    FileIO::read( mLocalValues, fileName, common::ScalarType::INTERNAL, first, n );

    return mLocalValues.size();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void DenseVector<ValueType>::clearValues()
{
    mLocalValues.clear();
}

/* ========================================================================= */
/*       Factory methods (must be provided to register)                      */
/* ========================================================================= */

template<typename ValueType>
_Vector* DenseVector<ValueType>::create()
{
    return new DenseVector<ValueType>();
}

template<typename ValueType>
VectorCreateKeyType DenseVector<ValueType>::createValue()
{
    return VectorCreateKeyType( VectorKind::DENSE, common::getScalarType<ValueType>() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const DenseVector<ValueType>& other ) :

    Vector<ValueType>( other )

{
    // implementation here can be simpler as DenseVector( const Vector& other )
    SCAI_LOG_INFO( logger,
                   "Copy of vector of global size " << size() << ", local size " << getDistribution().getLocalSize() )
    mLocalValues = other.getLocalValues();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DenseVector<ValueType>&& other ) noexcept :

    Vector<ValueType>( other )

{
    SCAI_LOG_INFO( logger,
                   "move of vector of global size " << size() << ", local size " << getDistribution().getLocalSize() )

    mLocalValues = std::move( other.mLocalValues );
    other.allocate( 0 );
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DenseVector, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
