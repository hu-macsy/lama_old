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
#include <scai/lama/VectorAssemblyAccess.hpp>

// local library
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
using utilskernel::LArray;

using namespace hmemo;
using namespace dmemo;

namespace lama
{

/* ------------------------------------------------------------------------- */
/*  Constructors of DenseVector<ValueType>                                   */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector() :
    Vector<ValueType>( 0 ),
    mLocalValues()
{
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( ContextPtr context ) :

    Vector<ValueType>( 0, context ),
    mLocalValues()
{
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const IndexType size ) :

    Vector<ValueType>( size ), 
    mLocalValues( size )

{
    SCAI_LOG_INFO( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">( size = " 
                            << size << " ), undefined values" )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution ) :

    Vector<ValueType>( distribution ), 
    mLocalValues( distribution->getLocalSize() )
{
    SCAI_LOG_INFO( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">( dist = " 
                   << *distribution << "), undefined values" )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const IndexType size, ContextPtr context ) : 

    Vector<ValueType>( size, context ), 
    mLocalValues( size )

{
    SCAI_LOG_INFO( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">( size = " 
                            << size << " ), undefined values, @ctx = " << getContext() )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, ContextPtr context ) : 

    Vector<ValueType>( distribution, context ), 
    mLocalValues( distribution->getLocalSize() )

{
    SCAI_LOG_INFO( logger, "DenseVector<" << common::TypeTraits<ValueType>::id() << ">( dist = " 
                   << *distribution << "), undefined values, @ctx = " << getContext() )
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
                   "Construct dense vector, size = " << distribution->getGlobalSize() << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", value = " << value )
}

template <typename ValueType>
void DenseVector<ValueType>::fillLinearValues( const ValueType startValue, const ValueType inc )
{
    const Distribution& dist = getDistribution();

    if (dist.isReplicated() )
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
    HArrayUtils::assignScalar( mLocalValues, inc, BinaryOp::MULT, context );
    HArrayUtils::assignScalar( mLocalValues, startValue, BinaryOp::ADD, context );
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const _Vector& other ) :

    Vector<ValueType>( other )

{
    allocate();
    assign( other );
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const _Vector& other, DistributionPtr distribution ) :

    Vector<ValueType>( other )

{
    allocate();       // make sure that local values array fits distribution, context
    assign( other );
    redistribute( distribution );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const std::string& filename ) : 

    Vector<ValueType>( 0 )

{
    SCAI_LOG_INFO( logger, "Construct dense vector from file " << filename )
    readFromFile( filename );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::fillRandom( const IndexType bound )
{
    mLocalValues.setRandom( bound, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::fillSparseRandom( const float fillRate, const IndexType bound )
{
    mLocalValues.setSparseRandom( fillRate, bound, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, const _HArray& localValues ) :

    Vector<ValueType>( distribution ),
    mLocalValues( localValues )
{
    SCAI_ASSERT_EQ_ERROR( localValues.size(), distribution->getLocalSize(), "size mismatch" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const _HArray& values ) :

    Vector<ValueType>( DistributionPtr( new NoDistribution( values.size() ) ) ),
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
DenseVector<ValueType>::DenseVector( const Expression_SV<ValueType>& expression ):

    Vector<ValueType>( expression.getArg2() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x )" )
    Vector<ValueType>::operator=( expression );
}

// linear algebra expression: a+x/x+a
template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SV_S<ValueType>& expression ) :

    Vector<ValueType>( expression.getArg1().getArg2() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x + beta)" )
    Vector<ValueType>::operator=( expression );
}

// linear algebra expression: x*y

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_VV<ValueType>& expression ) :

    Vector<ValueType>( expression.getArg1() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( x * y )" )
    Expression_SVV<ValueType> tmpExp( ValueType( 1 ), expression );
    Vector<ValueType>::operator=( tmpExp );
}

// linear algebra expression: s*x*y
template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SVV<ValueType>& expression ) :

    Vector<ValueType>( expression.getArg2().getArg1() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * y )" )
    Vector<ValueType>::operator=( expression );
}

// linear algebra expression: a*x+b*y, inherit distribution/context from vector x

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SV_SV<ValueType>& expression ) :

    Vector<ValueType>( expression.getArg1().getArg2() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x + beta * y )" )
    Vector<ValueType>::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SMV_SV<ValueType>& expression ) 
{
    const Matrix<ValueType>& m = expression.getArg1().getArg2().getArg1();
    this->setSpace( m.getRowDistributionPtr(), m.getContextPtr()  );
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x + b * y )" )
    Vector<ValueType>::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SVM_SV<ValueType>& expression ) 
{
    const Matrix<ValueType>& m = expression.getArg1().getArg2().getArg2();
    this->setSpace( m.getColDistributionPtr(), m.getContextPtr()  );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A + b * y )" )
    Vector<ValueType>::operator=( expression );
}

// linear algebra expression: a*A*x, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SMV<ValueType>& expression ) : 

    Vector<ValueType>( expression.getArg2().getArg1().getRowDistributionPtr(),
            expression.getArg2().getArg1().getContextPtr() )
{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x )" )
    Vector<ValueType>::operator=( expression );
}

// linear algebra expression: a*x*A, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression_SVM<ValueType>& expression ) :

    Vector<ValueType>( expression.getArg2().getArg2().getColDistributionPtr(),
            expression.getArg2().getArg2().getContextPtr() )

{
    allocate();
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A )" )
    Vector<ValueType>::operator=( expression );
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

    std::unique_ptr<ValueType[]> splitValues( new ValueType[ nPartitions + 1 ] );

    {
        ReadAccess<ValueType> rSortedValues( sortedValues );
        getSplitValues( splitValues.get(), comm, rSortedValues.get(), sortedValues.size(), ascending );
    }

    std::unique_ptr<IndexType[]> quantities( new IndexType[ nPartitions ] );

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

    SCAI_ASSERT_NE_ERROR( getDistribution().getBlockDistributionSize(), nIndex,
                          "scan only supported for block distribution" )

    const Communicator& comm = getDistribution().getCommunicator();
 
    HArray<ValueType> prefixValues;
    
    ValueType val = mLocalValues.sum();

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

    SCAI_ASSERT_NE_ERROR( other.getDistribution().getBlockDistributionSize(), nIndex,
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

    HArrayUtils::assign( mLocalValues, values, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::fillSparseData( 
    const HArray<IndexType>& nonZeroIndexes, 
    const _HArray& nonZeroValues,
    const BinaryOp op )
{
    // Note: scatter checks for legal indexes

    HArrayUtils::scatter( mLocalValues, nonZeroIndexes, false, nonZeroValues, op, getContextPtr() );
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
    std::unique_ptr<DenseVector<ValueType> > vector( new DenseVector<ValueType>() );
    vector->setContextPtr( this->getContextPtr() );
    return vector.release();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::concatenate( dmemo::DistributionPtr dist, const std::vector<const _Vector*>& vectors ) 
{
    DenseVector<ValueType> newVector( dist );
    {
        VectorAssemblyAccess<ValueType> assembly( newVector );

        IndexType offset = 0;

        for ( size_t k = 0; k < vectors.size(); ++k )
        {
            const _Vector& v = *vectors[k];

            if ( offset + v.size() > dist->getGlobalSize() )
            {
                COMMON_THROWEXCEPTION( "concatenate fails, exceeds global size of target vector" )
            }

            HArray<ValueType> localData;
  
            v.buildLocalValues( localData );

            ReadAccess<ValueType> rData( localData );

            for ( IndexType i = 0; i < rData.size(); ++i )
            {
                IndexType globalI = v.getDistribution().local2global( i );

                SCAI_LOG_ERROR( logger, "push " << offset + globalI << ", " << rData[i] )

                assembly.push( offset + globalI, rData[i] );
            }
 
            offset += v.size();
        }
    }

    swap( newVector );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseVector<ValueType>::getValue( IndexType globalIndex ) const
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

    if ( localIndex != nIndex )
    {
        mLocalValues[localIndex] = value;
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseVector<ValueType>::min() const
{
    // Note: min returns the maximal representation value on zero-sized vectors, TypeTraits<ValueType>::getMax()
    ValueType localMin = mLocalValues.min();
    return getDistribution().getCommunicator().min( localMin );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseVector<ValueType>::max() const
{
    // Note: max returns the minimal representation value on zero-sized vectors
    ValueType localMax = mLocalValues.max();
    return getDistribution().getCommunicator().max( localMax );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseVector<ValueType>::l1Norm() const
{
    NormType<ValueType> localL1Norm = mLocalValues.l1Norm();
    return getDistribution().getCommunicator().sum( localL1Norm );
}

/*---------------------------------------------------------------------------*/
template<typename ValueType>
ValueType DenseVector<ValueType>::sum() const
{
    ValueType localsum = mLocalValues.sum();
    return getDistribution().getCommunicator().sum( localsum );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseVector<ValueType>::l2Norm() const
{
    // Note: we do not call l2Norm here for mLocalValues to avoid sqrt
    NormType<ValueType> localDotProduct = mLocalValues.dotProduct( mLocalValues );
    NormType<ValueType> globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    return common::Math::sqrt( globalDotProduct );
}

template<>
IndexType DenseVector<IndexType>::l2Norm() const
{
    // Note: we do not call l2Norm here for mLocalValues to avoid sqrt

    double localDotProduct = mLocalValues.dotProduct( mLocalValues );
    double globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    return IndexType( common::Math::sqrt( globalDotProduct ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseVector<ValueType>::maxNorm() const
{
    NormType<ValueType> localMaxNorm = mLocalValues.maxNorm();
    const Communicator& comm = getDistribution().getCommunicator();
    NormType<ValueType> globalMaxNorm = comm.max( localMaxNorm );
    SCAI_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm: " << localMaxNorm
                   << ", max norm global = " << globalMaxNorm )
    return globalMaxNorm;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseVector<ValueType>::maxDiffNorm( const Vector<ValueType>& other ) const
{
    SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(), "maxDiffNorm: distribution mismatch" )

    if ( other.getVectorKind() != VectorKind::DENSE    )
    {
        // ToDo: this might be improved but requires some more kernel support

        DenseVector<ValueType> tmpOther( other, getDistributionPtr() );
        return maxDiffNorm( tmpOther );
    }

    SCAI_ASSERT_DEBUG( dynamic_cast<const DenseVector<ValueType>*>( &other ), "wrong cast" )

    const DenseVector<ValueType>& denseOther = reinterpret_cast<const DenseVector<ValueType>&>( other );

    NormType<ValueType> localMaxNorm = mLocalValues.maxDiffNorm( denseOther.getLocalValues() );

    const Communicator& comm = getDistribution().getCommunicator();

    NormType<ValueType> globalMaxNorm = comm.max( localMaxNorm );

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
        const DenseVector<ValueType>& denseOther = reinterpret_cast<const DenseVector<ValueType>&>( other );
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

        const SparseVector<ValueType>& sparseX = reinterpret_cast<const SparseVector<ValueType>&>( x );

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

        const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );

        const LArray<ValueType>& xValues = denseX.getLocalValues();

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
        DenseVector<ValueType> xTmp( x );
        vectorPlusVector( alpha, xTmp, beta, y );
        return;
    }

    if ( y.getVectorKind() != VectorKind::DENSE  )
    {
        SCAI_UNSUPPORTED( "vectorPlusVector: use temporary DenseVector<" << getValueType() << "> for y = " << y )
        DenseVector<ValueType> yTmp( y );
        vectorPlusVector( alpha, x, beta, yTmp );
        return;
    }

    SCAI_REGION( "Vector.Dense.VplusV" )

    if ( x.getDistribution() != getDistribution() )
    {
        allocate( x.getDistributionPtr() );   // result inherits distribution of operands
    }

    const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>& denseY = reinterpret_cast<const DenseVector<ValueType>&>( y );

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

    allocate( x.getDistributionPtr() );   // result inherits (same) space of operands, mostly its same

    if ( x.getVectorKind() != VectorKind::DENSE || x.getValueType() != getValueType() )
    {
        DenseVector<ValueType> xTmp( x );
        vectorTimesVector( alpha, xTmp, y );
        return;
    }

    if ( y.getVectorKind() != VectorKind::DENSE || y.getValueType() != getValueType() )
    {
        DenseVector<ValueType> yTmp( y );
        vectorTimesVector( alpha, x, yTmp );
        return;
    }

    const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>& denseY = reinterpret_cast<const DenseVector<ValueType>&>( y );

    SCAI_ASSERT_EQ_DEBUG( mLocalValues.size(), denseX.mLocalValues.size(), "serious space mismatch" )
    SCAI_ASSERT_EQ_DEBUG( mLocalValues.size(), denseY.mLocalValues.size(), "serious space mismatch" )

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

        DenseVector<ValueType> xTmp( x );
        vectorPlusScalar( alpha, xTmp, beta );
        return;
    }

    const DenseVector<ValueType>& denseX = dynamic_cast<const DenseVector<ValueType>&>( x );

    SCAI_LOG_INFO( logger, "z = alpha * x + beta, z = " << *this << ", alpha=  " << alpha
                   << " , x = " << x << " , beta = " << beta )
    SCAI_LOG_DEBUG( logger, "dist of x = " << x.getDistribution() )

    allocate( x.getDistributionPtr() );

    SCAI_ASSERT_EQ_DEBUG( mLocalValues.size(), denseX.mLocalValues.size(), "serious mismatch" )

    SCAI_LOG_DEBUG( logger, "call arrayPlusScalar" )

    utilskernel::HArrayUtils::arrayPlusScalar( mLocalValues, alpha, denseX.mLocalValues, beta, getContextPtr() );
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
ValueType DenseVector<ValueType>::dotProduct( const _Vector& other ) const
{
    SCAI_REGION( "Vector.Dense.dotP" )
    SCAI_LOG_INFO( logger, "Calculating dot product for " << *this << " * " << other )

    // add other->getVectorKind() == VectorKind::DENSE, if sparse is also supported

    SCAI_ASSERT_EQ_ERROR( getValueType(), other.getValueType(),
                          "dotProduct not supported for different value types. "
                          << *this << " x " << other )

    SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(),
                          "dotProduct not supported for vectors with different distributions. "
                          << *this  << " x " << other )

    ValueType localDotProduct;

    if ( other.getVectorKind() == VectorKind::DENSE )
    {
        SCAI_ASSERT_DEBUG( dynamic_cast<const DenseVector<ValueType>*>( &other ), "dynamic cast failed, other = " << other )

        const DenseVector<ValueType>& denseOther = reinterpret_cast<const DenseVector<ValueType>&>( other );

        localDotProduct = mLocalValues.dotProduct( denseOther.getLocalValues() );
    }
    else
    {
        SCAI_ASSERT_DEBUG( dynamic_cast<const SparseVector<ValueType>*>( &other ), "dynamic cast failed, other = " << other )

        const SparseVector<ValueType>& sparseOther = reinterpret_cast<const SparseVector<ValueType>&>( other );
    
        LArray<ValueType> myValues;   // build values at same position as sparse vector

        gatherLocalValues( myValues, sparseOther.getNonZeroIndexes(), BinaryOp::COPY, getContextPtr() );

        localDotProduct = myValues.dotProduct( sparseOther.getNonZeroValues() );
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

template<typename ValueType, typename TList> struct VectorWrapperT;

template<typename ValueType>
struct VectorWrapperT<ValueType, common::mepr::NullType>
{
    static void assign(
        DenseVector<ValueType>& target,
        const _Vector& source )
    {
        COMMON_THROWEXCEPTION( "vector assign = " << target << ", source = " << source  )
    }
};

template<typename ValueType, typename H, typename Tail>
struct VectorWrapperT< ValueType, common::mepr::TypeList<H, Tail> >
{
    static void assign(
        DenseVector<ValueType>& target,
        const _Vector& source )
    {
        if ( common::getScalarType<H>() ==  source.getValueType() )
        {
            if ( source.getVectorKind() == VectorKind::SPARSE )
            {
                const SparseVector<H>& typedSource = reinterpret_cast<const SparseVector<H>&>( source );
                target.assignImpl( typedSource );
            }
            else if ( source.getVectorKind() == VectorKind::DENSE )
            {
                const DenseVector<H>& typedSource = reinterpret_cast<const DenseVector<H>&>( source );
                target.assignImpl( typedSource );
            }
            else
            {
                COMMON_THROWEXCEPTION( "unsupported vector kind for assign to sparse vector" )
            }
        }
        else
        {
            VectorWrapperT< ValueType, Tail >::assign( target, source );
        }
    }
};

template<typename ValueType>
void DenseVector<ValueType>::assign( const _Vector& other )
{
    // use metaprogramming to call the correspoing assignImpl method

    VectorWrapperT<ValueType, SCAI_ARRAY_TYPES_HOST_LIST>::assign( *this, other );
}

template<typename ValueType>
template<typename OtherValueType>
void DenseVector<ValueType>::assignImpl( const SparseVector<OtherValueType>& other )
{
    allocate( other.getDistributionPtr() );
    setScalar( other.getZero() );
    fillSparseData( other.getNonZeroIndexes(), other.getNonZeroValues(), BinaryOp::COPY );
}

template<typename ValueType>
template<typename OtherValueType>
void DenseVector<ValueType>::assignImpl( const DenseVector<OtherValueType>& other )
{
    allocate( other.getDistributionPtr() );
    setDenseValues( other.getLocalValues() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::buildLocalValues(

     _HArray& localValues, 
     const BinaryOp op,
     ContextPtr prefLoc ) const

{
    if ( op == BinaryOp::COPY )
    {
        HArrayUtils::assign( localValues, mLocalValues, prefLoc );
    }
    else
    {
        SCAI_ASSERT_EQ_ERROR( localValues.size(), mLocalValues.size(), "size mismatch" )
        HArrayUtils::setArray( localValues, mLocalValues, op, prefLoc );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::gatherLocalValues(

     _HArray& localValues,
     const HArray<IndexType>& indexes,
     const BinaryOp op,
     ContextPtr prefLoc ) const
{
    HArrayUtils::gather( localValues, mLocalValues, indexes, op, prefLoc );
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
        DenseVector<ValueType> tmpX( x );
        unaryOp( x, op );
        return;
    }

    allocate( x.getDistributionPtr() );

    const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );

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
        DenseVector<ValueType> tmpX( x );
        binaryOp( tmpX, op, y );
        return;
    }

    if ( y.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_UNSUPPORTED( "denseR = x " << op << " sparseY, tmpDenseVector = sparseY added" )
        DenseVector<ValueType> tmpY( y );
        binaryOp( x, op, tmpY );
        return;
    }

    allocate( x.getDistributionPtr() );

    const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );
    const DenseVector<ValueType>& denseY = reinterpret_cast<const DenseVector<ValueType>&>( y );

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
        DenseVector<ValueType> tmpX( x );
        binaryOpScalar( tmpX, alpha, op, swap );
        return;
    }

    if ( &x != this )
    {
        allocate( x.getDistributionPtr() ); 
    }
    
    const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );
    
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

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DenseVector, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
