/**
 * @file SparseVector.cpp
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
 * @brief Implementations of constructors/methods for class SparseVector.
 * @author Thomas Brandes
 * @date 16.01.2017
 */

// hpp
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/DenseVector.hpp>

// local library
#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/expression/Expression.hpp>

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/PartitionIO.hpp>

#include <scai/lama/mepr/VectorWrapper.hpp>

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/freeFunction.hpp>
#include <scai/common/BinaryOp.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/dmemo/RedistributePlan.hpp>
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

using common::Math;
using common::TypeTraits;
using utilskernel::HArrayUtils;

using namespace hmemo;
using namespace dmemo;

namespace lama
{

/* ------------------------------------------------------------------------- */
/*  Implementation of constructors for SparseVector                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::SparseVector( ContextPtr context ) :

    Vector<ValueType>( 0, context ),
    mNonZeroIndexes( context ),
    mNonZeroValues( context ),
    mZeroValue( 0 )
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const IndexType n, ValueType zero, ContextPtr context ) :

    Vector<ValueType>( n, context ),
    mNonZeroIndexes( context ),
    mNonZeroValues( context ),
    mZeroValue( zero )
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( DistributionPtr distribution, ValueType zero, ContextPtr context ) : 

    Vector<ValueType>( distribution, context ), 
    mNonZeroIndexes( context ), 
    mNonZeroValues( context ),
    mZeroValue( zero )

{
    SCAI_LOG_INFO( logger, "Construct sparse vector on context = " << context << ", size = " << distribution->getGlobalSize()
                         << ", distribution = " << *distribution << ", all zero" )
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const DenseVector<ValueType>& other, const ValueType zeroValue ) :

    Vector<ValueType>( other ),
    mZeroValue( zeroValue )

{
    // Note: setDenseValues uses the zero value
    setDenseValuesImpl( other.getLocalValues() );
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( SparseVector<ValueType>&& other ) noexcept :

    Vector<ValueType>( other )

{
    mZeroValue = other.mZeroValue;
    mNonZeroIndexes = std::move( other.mNonZeroIndexes );
    mNonZeroValues  = std::move( other.mNonZeroValues );

    // Note: other vector remains as a constant 'zero' vector
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::assign( const _Vector& other )
{
    // translate virtual call to specific template call via wrapper

    mepr::VectorWrapper<SparseVector, SCAI_ARRAY_TYPES_HOST_LIST>::assignImpl( this, other );
}

template<typename ValueType>
template<typename OtherValueType>
void SparseVector<ValueType>::assignImpl( const Vector<OtherValueType>& other )
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
void SparseVector<ValueType>::assignSparse( const SparseVector<OtherValueType>& other )
{
    SCAI_REGION( "Vector.Sparse.assginSparse" )

    SCAI_LOG_INFO( logger, "sparseVector<" << common::TypeTraits<ValueType>::id() << "> = "
                        << "sparseVector<" << common::TypeTraits<OtherValueType>::id() )

    if ( static_cast<const _Vector*>( &other ) == this )
    {
        return;
    } 

    allocate( other.getDistributionPtr() );

    mZeroValue = static_cast<ValueType>( other.getZero() );

    HArrayUtils::setArray( mNonZeroIndexes, other.getNonZeroIndexes(), common::BinaryOp::COPY, getContextPtr() );
    HArrayUtils::_setArray( mNonZeroValues, other.getNonZeroValues(), common::BinaryOp::COPY, getContextPtr() );
}

template<typename ValueType>
template<typename OtherValueType>
void SparseVector<ValueType>::assignDense( const DenseVector<OtherValueType>& other )
{   
    SCAI_REGION( "Vector.Sparse.assignDense" )

    SCAI_LOG_INFO( logger, "sparseVector<" << common::TypeTraits<ValueType>::id() << "> = "
                        << "denseVector<" << common::TypeTraits<OtherValueType>::id() << ">" )

    allocate( other.getDistributionPtr() );
    setDenseValuesImpl( other.getLocalValues() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::SparseVector( 
    dmemo::DistributionPtr distribution,
    hmemo::HArray<IndexType> nonZeroIndexes, 
    hmemo::HArray<ValueType> nonZeroValues, 
    const ValueType zero,
    hmemo::ContextPtr context ) :

    Vector<ValueType>( distribution, context ),
    mNonZeroIndexes( std::move( nonZeroIndexes ) ),
    mNonZeroValues( std::move( nonZeroValues ) ),
    mZeroValue( zero )
{
    SCAI_ASSERT_EQ_ERROR( mNonZeroIndexes.size(), mNonZeroValues.size(), 
                          "#indexes and #values must be same for sparse vector" )

    const IndexType size = getDistribution().getLocalSize();

    bool isValid = HArrayUtils::validIndexes( mNonZeroIndexes, size, getContextPtr() );

    if ( !isValid )
    {
        COMMON_THROWEXCEPTION( "at least one illegal index, local size = " << size )
    }

    HArrayUtils::sortSparseEntries( mNonZeroIndexes, mNonZeroValues, true, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::SparseVector(
    const IndexType n,
    hmemo::HArray<IndexType> nonZeroIndexes, 
    hmemo::HArray<ValueType> nonZeroValues, 
    const ValueType zero,
    hmemo::ContextPtr context ) :

    Vector<ValueType>( std::make_shared<NoDistribution>( n ), context ),
    mNonZeroIndexes( std::move( nonZeroIndexes ) ),
    mNonZeroValues( std::move( nonZeroValues ) ),
    mZeroValue( zero )
{
    SCAI_ASSERT_EQ_ERROR( mNonZeroIndexes.size(), mNonZeroValues.size(), 
                          "#indexes and #values must be same for sparse vector" )

    bool isValid = HArrayUtils::validIndexes( nonZeroIndexes, n, getContextPtr() );
    
    if ( !isValid )
    {   
        COMMON_THROWEXCEPTION( "at least one illegal index, size = " << n )
    }

    HArrayUtils::sortSparseEntries( mNonZeroIndexes, mNonZeroValues, true, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::fillRandom( const IndexType bound )
{
    const IndexType localSize = getDistribution().getLocalSize();
    auto localValues = utilskernel::randomHArray<ValueType>( localSize, bound, getContextPtr() );
    setDenseValuesImpl( localValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::fillSparseRandom( const float fillRate, const IndexType bound )
{
    SCAI_ASSERT_EQ_ERROR( 0, mNonZeroIndexes.size(), "SparseRandom illegal, vector has already non-zero elements" )

    const IndexType localSize = getDistribution().getLocalSize();

    HArrayUtils::randomSparseIndexes( mNonZeroIndexes, localSize, fillRate );
    mNonZeroValues.resize( mNonZeroIndexes.size() );
    HArrayUtils::fillRandom( mNonZeroValues, bound, 1.0f );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::~SparseVector()
{
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::splitUp( DistributionPtr& dist, 
                                       hmemo::HArray<IndexType>& nonZeroIndexes,
                                       hmemo::HArray<ValueType>& nonZeroValues,
                                       ValueType& zero )
{
    dist  = getDistributionPtr();

    nonZeroIndexes = std::move( mNonZeroIndexes );
    nonZeroValues  = std::move( mNonZeroValues );

    zero  = mZeroValue;

    // this sparse vector remains consistent
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>& SparseVector<ValueType>::operator=( const SparseVector<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "SparseVector<" << TypeTraits<ValueType>::id() << "> = " <<
                   "SparseVector<" << TypeTraits<ValueType>::id() << ">" )

    assign( other );
    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>& SparseVector<ValueType>::operator=( SparseVector<ValueType>&& other ) noexcept
{
    // Note: we do not inherit the context from the other one

    setDistributionPtr( other.getDistributionPtr() );

    mZeroValue = other.mZeroValue;
    mNonZeroIndexes = std::move( other.mNonZeroIndexes );
    mNonZeroValues  = std::move( other.mNonZeroValues );

    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool SparseVector<ValueType>::isConsistent() const
{
    // for a spare vector we have to check that the nonZeroIndexes are 
    // all legal and well sorted (strongly increasing).

    const Distribution& dist = getDistribution();

    IndexType consistencyErrors = 0;

    const IndexType localSize = dist.getLocalSize();

    const HArray<IndexType>& nonZeroIndexes = getNonZeroIndexes();

    if ( getNonZeroValues().size() != nonZeroIndexes.size() )
    {
        SCAI_LOG_INFO( logger, "sizes of nonZeroValues/nonZeroIndexes do not match" )
        consistencyErrors++;
    }

    if ( !HArrayUtils::validIndexes( nonZeroIndexes, localSize ) )
    {
        SCAI_LOG_INFO( logger, "sparse indexes not valid for localSize = " << localSize )
        consistencyErrors++;
    }

    // indexes for the non-zero values must be sorted

    if ( !HArrayUtils::isSorted( nonZeroIndexes, common::CompareOp::LT ) )
    {
        SCAI_LOG_INFO( logger, "sparse indexes not strong increasing" )
        consistencyErrors++;
    }

    // not checked: there should be no double values in indexes

    // use communicator for global reduction to make sure that all processors return same value.

    consistencyErrors = dist.getCommunicator().sum( consistencyErrors );

    return 0 == consistencyErrors;
}

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
void SparseVector<ValueType>::buildLocalValues( 
    _HArray& values, 
    const common::BinaryOp op,
    ContextPtr loc ) const
{
    SCAI_LOG_INFO( logger, *this << ": build local values, op = " << op << " into " << values )

    const bool UNIQUE = true;

    // size of values will be local size of vector

    const IndexType size = getDistribution().getLocalSize();

    // convert the local sparse data to local dense data

    if ( op == common::BinaryOp::COPY )
    {
        // build values array from scratch 

        HArrayUtils::_setSameValue( values, size, mZeroValue, getContextPtr() );
        HArrayUtils::_scatter( values, mNonZeroIndexes, UNIQUE, mNonZeroValues, op, getContextPtr() );
    }
    else if ( mZeroValue == common::zeroBinary<ValueType>( op ) ) 
    {
        // ZERO element of this sparse vector is ZERO element for op, that is fine, we only apply non-zero values

        SCAI_ASSERT_EQ_ERROR( values.size(), size, "size mismatch" )
        HArrayUtils::_scatter( values, mNonZeroIndexes, UNIQUE, mNonZeroValues, op, loc );
    }
    else 
    {
        // temporary array needed for this operation

        SCAI_UNSUPPORTED( *this << ", mZero = " << mZeroValue << " is not ZERO element of " << op << ", temporary dense values are built" )

        auto denseArray = utilskernel::fillHArray<ValueType>( size, mZeroValue, loc );
        HArrayUtils::scatter( denseArray, mNonZeroIndexes, UNIQUE, mNonZeroValues, common::BinaryOp::COPY, loc );
        HArrayUtils::_setArray( values, denseArray, op, loc );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::gatherLocalValues(
    _HArray& values,
    const HArray<IndexType>& indexes,
    const common::BinaryOp op,
    ContextPtr loc ) const
{
    // for operations on sparse vectors with same patterns it is worth to check for same indexes

    bool samePattern = mNonZeroIndexes.size() == indexes.size();

    if ( samePattern )
    {
        samePattern = HArrayUtils::all( mNonZeroIndexes, common::CompareOp::EQ, indexes, loc );
    }

    if ( samePattern )
    {
        HArrayUtils::_assign( values, mNonZeroValues, loc );
    }
    else
    {
        // this operation is not very efficient

        HArrayUtils::_sparseGather( values, mZeroValue, mNonZeroValues, mNonZeroIndexes, indexes, op, loc );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setDenseValues( const _HArray& localValues )
{
    mepr::VectorWrapper<SparseVector, SCAI_ARRAY_TYPES_HOST_LIST>::setDenseValuesImpl( this, localValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void SparseVector<ValueType>::setDenseValuesImpl( const hmemo::HArray<OtherValueType>& values )
{
    SCAI_REGION( "Vector.Sparse.setDense" )

    const IndexType size = getDistribution().getLocalSize();

    SCAI_ASSERT_EQ_ERROR( size, values.size(), "size of local values does not match local size of vector" )

    OtherValueType zero = static_cast<OtherValueType>( mZeroValue );

    SCAI_LOG_INFO( logger, "setDenseValues<" << common::TypeTraits<ValueType>::id() << ", " << common::TypeTraits<OtherValueType>::id() << ">"
                            << ", local size = " << size << ", zero = " << mZeroValue )

    utilskernel::HArrayUtils::buildSparseArray( mNonZeroValues, mNonZeroIndexes, values, zero, getContextPtr() );

    SCAI_ASSERT_EQ_DEBUG( mNonZeroValues.size(), mNonZeroIndexes.size(), "serious mismatch for sparse arrays." )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::swapSparseValues( HArray<IndexType>& nonZeroIndexes, HArray<ValueType>& nonZeroValues )
{
    SCAI_ASSERT_EQ_ERROR( nonZeroIndexes.size(), nonZeroValues.size(), "size mismatch for arrays with non-zero indexes/values" )

    const IndexType size = getDistribution().getLocalSize();

    bool isValid = HArrayUtils::validIndexes( nonZeroIndexes, size, getContextPtr() );

    if ( !isValid )
    {
        COMMON_THROWEXCEPTION( "at least one illegal index, local size = " << size )
    }

    mNonZeroIndexes.swap( nonZeroIndexes );
    mNonZeroValues.swap( nonZeroValues );

    HArrayUtils::sortSparseEntries( mNonZeroIndexes, mNonZeroValues, true, getContextPtr() );

    SCAI_ASSERT_DEBUG( HArrayUtils::isSorted( mNonZeroIndexes, common::CompareOp::LT ), "sort sparse entries failed" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::fillSparseData( const HArray<IndexType>& nonZeroIndexes, const _HArray& nonZeroValues, common::BinaryOp op )
{
    SCAI_LOG_INFO( logger, "fillSparseData, nnz = " << nonZeroIndexes.size() )

    SCAI_ASSERT_EQ_ERROR( nonZeroIndexes.size(), nonZeroValues.size(), "arrays for sparse indexes and values must have same size" )

    if ( nonZeroIndexes.size() == 0 )
    {
        return;
    }
   
    const IndexType localSize = getDistribution().getLocalSize();

    bool isValid = HArrayUtils::validIndexes( nonZeroIndexes, localSize, getContextPtr() );

    if ( !isValid )
    {
        COMMON_THROWEXCEPTION( "at least one illegal index, local size = " << localSize )
    }

    if ( mNonZeroIndexes.size() !=  0 )
    {
        // sort the new sparse entries so merge is more efficient

        HArray<IndexType> newIndexes;
        HArray<ValueType> newValues;

        HArrayUtils::setArray( newIndexes, nonZeroIndexes, common::BinaryOp::COPY, getContextPtr() );
        HArrayUtils::_setArray( newValues, nonZeroValues, common::BinaryOp::COPY, getContextPtr() );

        // then we have to merge zero indexes

        HArrayUtils::sortSparseEntries( newIndexes, newValues, true, getContextPtr() );

        HArray<IndexType> resultIndexes;
        HArray<ValueType> resultValues;

        // Note: mergeSparse will also eleminate double values in one single input set

        HArrayUtils::mergeSparse( resultIndexes, resultValues,
                                  mNonZeroIndexes, mNonZeroValues,
                                  newIndexes, newValues, op );

        mNonZeroIndexes.swap( resultIndexes );
        mNonZeroValues.swap( resultValues );
    }
    else
    {
        HArrayUtils::setArray( mNonZeroIndexes, nonZeroIndexes, common::BinaryOp::COPY, getContextPtr() );
        HArrayUtils::_setArray( mNonZeroValues, nonZeroValues, common::BinaryOp::COPY, getContextPtr() );

        HArrayUtils::sortSparseEntries( mNonZeroIndexes, mNonZeroValues, true, getContextPtr() );
        HArrayUtils::elimDoubles( mNonZeroIndexes, mNonZeroValues, op );
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void SparseVector<ValueType>::readFromFile( FileIO& file )
{
    auto comm = file.getCommunicatorPtr();

    if ( file.getDistributedIOMode() == DistributedIOMode::MASTER )
    {
        const PartitionId MASTER = 0;
        const PartitionId myRank = comm->getRank();

        IndexType globalSize;

        if ( myRank == MASTER )
        {
            SCAI_LOG_INFO( logger, *comm << ": read sing sparse array from file " )
            file.readSparse( globalSize, &mZeroValue, mNonZeroIndexes, mNonZeroValues );
        }
        else
        {
            mNonZeroIndexes.clear();     
            mNonZeroValues.clear();     
        }

        comm->bcast( &mZeroValue, 1, MASTER );
        comm->bcast( &globalSize, 1, MASTER );

        setDistributionPtr( singleDistribution( globalSize, comm, MASTER ) );
    }
    else
    {
        // INDEPENDENT or COLLECTIVE 

        IndexType localSize;

        file.readSparse( localSize, &mZeroValue, mNonZeroIndexes, mNonZeroValues );
        setDistributionPtr( genBlockDistributionBySize( localSize, comm ) );
    }

    HArrayUtils::sortSparseEntries( mNonZeroIndexes, mNonZeroValues, true, getContextPtr() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void SparseVector<ValueType>::writeLocalToFile( FileIO& file ) const
{
    // in contrary to writeToFile( file ) here data is already redistributed correctly

    const IndexType localSize = getDistribution().getLocalSize();

    file.writeSparse( localSize, &mZeroValue, mNonZeroIndexes, mNonZeroValues );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void SparseVector<ValueType>::clearValues()
{
    mNonZeroValues.clear();
    mNonZeroIndexes.clear();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>* SparseVector<ValueType>::copy() const
{
    // create a new sparse vector with the copy constructor

    return new SparseVector<ValueType>( *this );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>* SparseVector<ValueType>::newVector() const
{
    std::unique_ptr<SparseVector<ValueType>> vector( new SparseVector<ValueType>( getContextPtr() ) );
    vector->setSameValue( this->getDistributionPtr(), 0 );
    return vector.release();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType SparseVector<ValueType>::getValue( IndexType globalIndex ) const
{
    ValueType myValue = 0;

    const IndexType localIndex = getDistribution().global2Local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": getValue( globalIndex = " << globalIndex << " ) -> local : " << localIndex )

    if ( localIndex != invalidIndex )
    {
        // we have here sparse data, so look for local index among the sparse indexes of non-zero values

        IndexType pos = HArrayUtils::findPosInSortedIndexes( mNonZeroIndexes, localIndex );

        if ( pos != invalidIndex )
        {
            myValue = mNonZeroValues[pos];
        }
        else
        {
            myValue = mZeroValue;
        }
    }

    ValueType allValue = getDistribution().getCommunicator().sum( myValue );

    // works also fine for replicated distributions with NoCommunicator

    SCAI_LOG_TRACE( logger, "myValue = " << myValue << ", allValue = " << allValue )

    return allValue;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setValue( const IndexType globalIndex, const ValueType value )
{
    SCAI_ASSERT_VALID_INDEX_ERROR( globalIndex, size(), "illegal index" )

    SCAI_LOG_TRACE( logger, *this << ": setValue( globalIndex = " << globalIndex << " ) = " <<  value )

    const IndexType localIndex = getDistribution().global2Local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": set @g " << globalIndex << " is @l " << localIndex << " : " << value )

    if ( localIndex != invalidIndex )
    {
        // This partition is the owner, add it locally

        IndexType pos = HArrayUtils::findPosInSortedIndexes( mNonZeroIndexes, localIndex );

        if ( pos != invalidIndex )
        {
            mNonZeroValues[pos] = value;
        }
        else
        {
            // add a new entry in mNonZeroIndexes, mNonZeroValues

            pos = HArrayUtils::insertSorted( mNonZeroIndexes, localIndex );
            SCAI_LOG_TRACE( logger, "setValue, local index = " << localIndex << " at pos = " << pos )
            HArrayUtils::insertAtPos( mNonZeroValues, pos, value );
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType SparseVector<ValueType>::min() const
{
    // Note: min returns the maximal representation value on zero-sized vectors, TypeTraits<ValueType>::getMax()

    ValueType localMin = HArrayUtils::min( mNonZeroValues );

    // if there are implicit zero values they must be used for min computation

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        localMin = Math::min( localMin, mZeroValue );
    }

    return getDistribution().getCommunicator().min( localMin );
}

#ifdef SCAI_COMPLEX_SUPPORTED
template<>
ComplexFloat SparseVector<ComplexFloat>::min() const
{
    COMMON_THROWEXCEPTION( "min unsupported on complex (float) (sparse) vector." )
}

template<>
ComplexDouble SparseVector<ComplexDouble>::min() const
{
    COMMON_THROWEXCEPTION( "min unsupported on complex (double) (sparse) vector." )
}

template<>
ComplexLongDouble SparseVector<ComplexLongDouble>::min() const
{
    COMMON_THROWEXCEPTION( "min unsupported on complex (long double) (sparse) vector." )
}
#endif

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType SparseVector<ValueType>::max() const
{
    // Note: max returns the minimal representation value on zero-sized vectors

    ValueType localMax = HArrayUtils::max( mNonZeroValues );

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        localMax = Math::max( localMax, mZeroValue );
    }

    return getDistribution().getCommunicator().max( localMax );
}

#ifdef SCAI_COMPLEX_SUPPORTED
template<>
ComplexFloat SparseVector<ComplexFloat>::max() const
{
    COMMON_THROWEXCEPTION( "max unsupported for complex vectors." )
}

template<>
ComplexDouble SparseVector<ComplexDouble>::max() const
{
    COMMON_THROWEXCEPTION( "max unsupported for complex vectors." )
}

template<>
ComplexLongDouble SparseVector<ComplexLongDouble>::max() const
{
    COMMON_THROWEXCEPTION( "max unsupported for complex vectors." )
}
#endif

/* ------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> SparseVector<ValueType>::l1Norm() const
{
    SCAI_REGION( "Vector.sparse.l1Norm" )

    auto localL1Norm = HArrayUtils::l1Norm( mNonZeroValues );

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        // ToDo: replace ABS with ASUM, is different for complex numbers

        RealType<ValueType> zeroNorm = common::applyUnary( common::UnaryOp::ABS, mZeroValue );
        localL1Norm += zeroNorm * RealType<ValueType>( nZero );
    }

    return getDistribution().getCommunicator().sum( localL1Norm );
}

/*---------------------------------------------------------------------------*/

template<typename ValueType>
ValueType SparseVector<ValueType>::sum() const
{
    ValueType localSum = HArrayUtils::sum( mNonZeroValues );

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        localSum += mZeroValue * ValueType( nZero );
    }

    return getDistribution().getCommunicator().sum( localSum );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> SparseVector<ValueType>::l2Norm() const
{
    SCAI_REGION( "Vector.sparse.l2Norm" )

    // Note: we do not call l2Norm here for mNonZeroValues to avoid sqrt

    RealType<ValueType> localDotProduct = HArrayUtils::dotProduct( mNonZeroValues, mNonZeroValues );

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        RealType<ValueType> zeroNorm = mZeroValue * Math::conj( mZeroValue ) * ValueType( nZero );
        localDotProduct += zeroNorm;
    }
 
    RealType<ValueType> globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );

    return Math::sqrt( globalDotProduct );
}

/* ------------------------------------------------------------------------- */

template<>
IndexType SparseVector<IndexType>::l2Norm() const
{
    SCAI_REGION( "Vector.sparse.l2Norm" )

    // Note: we do not call l2Norm here for mNonZeroValues to avoid sqrt

    double localDotProduct = static_cast<double>( HArrayUtils::dotProduct( mNonZeroValues, mNonZeroValues ) );
    double globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    return static_cast<IndexType>( Math::sqrt( globalDotProduct ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> SparseVector<ValueType>::maxNorm() const
{
    SCAI_REGION( "Vector.sparse.maxNorm" )

    auto localMaxNorm = HArrayUtils::maxNorm( mNonZeroValues );

    // the ZERO element must also be considered if at least one element is zero

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        SCAI_LOG_DEBUG( logger, "maxNorm, zero = " << mZeroValue << ", non-zero = " << localMaxNorm )
        localMaxNorm = Math::max( Math::abs( localMaxNorm ), Math::abs( mZeroValue ) );
    }

    const Communicator& comm = getDistribution().getCommunicator();

    RealType<ValueType> globalMaxNorm = comm.max( localMaxNorm );

    SCAI_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm: " << localMaxNorm
                   << ", max norm global = " << globalMaxNorm )
    return globalMaxNorm;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> SparseVector<ValueType>::maxDiffNorm( const Vector<ValueType>& other ) const
{
    SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(), "distribution mismatch for maxDiffNorm" )

    // ToDo: find some more efficient solutions wherever possible

    SCAI_LOG_INFO( logger, "maxDiffNorm, this = " << *this << ", other = " << other )

    SparseVector<ValueType> tmp;
    tmp.binaryOp( *this, common::BinaryOp::SUB, other );

    return tmp.maxNorm();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool SparseVector<ValueType>::all( const common::CompareOp op, const ValueType value ) const
{
    // all non-zero values must fulfill the condition

    bool localAll = HArrayUtils::allScalar( mNonZeroValues, op, value );

    if ( mNonZeroValues.size() != getDistribution().getLocalSize() )
    {
        // at least one entry has the ZERO value, so we compare it

        localAll = localAll && common::compare( mZeroValue, op, value );
    }

    bool globalAll = getDistribution().getCommunicator().all( localAll );

    return globalAll;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool SparseVector<ValueType>::all( const common::CompareOp op, const Vector<ValueType>& other ) const
{
    SCAI_ASSERT_EQ_ERROR( other.getDistribution(), getDistribution(), "distribution mismatch for all compare, op = " << op )

    if ( other.getVectorKind() == VectorKind::DENSE )
    {
        // dense vector can deal with sparse vector

        return other.all( op, *this );
    }

    if ( other.getValueType() != getValueType() )
    {
        return all( op, convert<SparseVector<ValueType>>( other ) );
    }

    // both vectors are sparse and have same value type

    bool localAll;

    const SparseVector<ValueType>& otherSparse = static_cast<const SparseVector<ValueType>& >( other );

    // ValueType otherZero = otherSparse.getZero().getValue<ValueType>();

    ValueType otherZero = otherSparse.getZero();

    IndexType n = HArrayUtils::allSparse( localAll,
                                          mNonZeroIndexes, mNonZeroValues, mZeroValue,
                                          otherSparse.getNonZeroIndexes(), otherSparse.getNonZeroValues(), otherZero, op );

    if ( n != getDistribution().getLocalSize() )
    {
        // at least at one position we use the zero values
    
        localAll = localAll && common::compare( mZeroValue, op, otherZero );
    }

    bool globalAll = getDistribution().getCommunicator().all( localAll );

    return globalAll;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::swap( _Vector& other )
{
    SCAI_LOG_DEBUG( logger, "swap:" << *this << " with " << other )

    SCAI_ASSERT_EQ_ERROR( getVectorKind(), other.getVectorKind(), "Swap only for same kind of vector allowed" )
    SCAI_ASSERT_EQ_ERROR( getValueType(), other.getValueType(), "Swap only for same value type of vector allowed" )

    SparseVector<ValueType>& typedOther = static_cast<SparseVector<ValueType>&>( other );

    _Vector::swapVector( other );   // swap sizes, distributions

    mNonZeroValues.swap( typedOther.mNonZeroValues );
    mNonZeroIndexes.swap( typedOther.mNonZeroIndexes );
    std::swap( mZeroValue, typedOther.mZeroValue );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::writeAt( std::ostream& stream ) const
{
    const Distribution& dist = getDistribution();

    stream << "SparseVector<" << getValueType() << ">" << "( size = " << size() << ", zero = " << mZeroValue 
           <<", local nnz = " << mNonZeroIndexes.size() << ", dist = " << dist << ", loc  = " << *getContextPtr() << " )";
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::unaryOp( const Vector<ValueType>& x, common::UnaryOp op )
{
    if ( x.getVectorKind() != VectorKind::SPARSE )
    {
        SCAI_UNSUPPORTED( "sparseVector = unaryOp( denseVector ), uses temporary dense vector" )
        DenseVector<ValueType> tmpResult;
        tmpResult.unaryOp( x, op );
        assign( tmpResult );
        return;
    }

    const SparseVector<ValueType>& sparseX = static_cast<const SparseVector<ValueType>&>( x );

    if ( &x != this )
    {
        // allocation and copy of non-zero indexes only required if there is no alias

        allocate( x.getDistributionPtr() );
        mNonZeroValues.resize( sparseX.mNonZeroValues.size() );
        mNonZeroIndexes = sparseX.mNonZeroIndexes;
    }

    mZeroValue = common::applyUnary( op, sparseX.mZeroValue );
    HArrayUtils::unaryOp( mNonZeroValues, sparseX.mNonZeroValues, op, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::binaryOp( const Vector<ValueType>& x, const common::BinaryOp op, const Vector<ValueType>& y )
{
    SCAI_ASSERT_EQ_ERROR( x.getDistribution(), y.getDistribution(), "serious space mismatch" )

    SCAI_LOG_INFO( logger, "binaryOp: this = x " << op << " y, with x = " << x << ", y = " << y );

    if ( x.getVectorKind() == VectorKind::SPARSE && y.getVectorKind() == VectorKind::SPARSE )
    {
        SCAI_REGION( "Vector.Sparse.binOpSS" )
        const SparseVector<ValueType>& sparseX = static_cast<const SparseVector<ValueType>&>( x );
        const SparseVector<ValueType>& sparseY = static_cast<const SparseVector<ValueType>&>( y );

        binaryOpSparse( sparseX, op, sparseY );
    }
    else
    {
        SCAI_REGION( "Vector.Sparse.binOpSD" )
        SCAI_UNSUPPORTED( "SparseVector<" << common::TypeTraits<ValueType>::id() << ">::binaryOp x " << op << " y" )
        DenseVector<ValueType> tmp;
        tmp.binaryOp( x, op, y );
        assign( tmp );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::binaryOpScalar( const Vector<ValueType>& x, const ValueType& alpha, const common::BinaryOp op, const bool swap )
{
    if ( swap )
    {
        SCAI_LOG_INFO( logger, "binaryOp: this = " << alpha << " " << op << " x, with x = " << x );
    }
    else
    {
        SCAI_LOG_INFO( logger, "binaryOp: this = x " << op << " " << alpha << ", with x = " << x );
    }

    if ( x.getVectorKind() != VectorKind::SPARSE )
    {
        // thisSparse = xDense op alpha -> tmpResult = xDense op alpha; thisSparse = tmpsDense 

        SCAI_UNSUPPORTED( "sparseVector = denseVector op scalar, uses temporaray dense vector" )
        DenseVector<ValueType> tmpResult;
        tmpResult.binaryOpScalar( x, alpha, op, swap );
        assign( tmpResult );
        return;
    }

    const SparseVector<ValueType> sparseX = static_cast<const SparseVector<ValueType>&>( x );

    if ( &x != this )
    {
        // allocation and copy of non-zero indexes only required if there is no alias

        allocate( x.getDistributionPtr() );
        mNonZeroValues.resize( sparseX.mNonZeroValues.size() );
        mNonZeroIndexes = sparseX.mNonZeroIndexes;
    }

    // the following code works for any kind of alias

    HArrayUtils::binaryOpScalar( mNonZeroValues, sparseX.mNonZeroValues, alpha, op, swap, getContextPtr() );

    if ( !swap )
    {
        mZeroValue = common::applyBinary( sparseX.mZeroValue, op, alpha );
    }
    else
    {
        mZeroValue = common::applyBinary( alpha, op, sparseX.mZeroValue );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::binaryOpSparseSamePattern( const SparseVector<ValueType>& x, const common::BinaryOp op, const SparseVector<ValueType>& y )
{
    // here we know that x.getNonZeroIndexes() == y.getNonZeroIndexes() )

    SCAI_REGION( "Vector.Sparse.binOpSparseSamePat" )

    SCAI_ASSERT_EQ_DEBUG( x.getDistribution(), y.getDistribution(), "serious space mismatch" )

    if ( getDistribution() != x.getDistribution() )
    {
        // there is no alias of this vector, neither with x nor with y

        allocate( x.getDistributionPtr() );
    }

    SCAI_LOG_INFO( logger, "binaryOpSparseSamePattern: this = x " << op << " y, with x = " << x << ", y = " << y );

    const HArray<ValueType>& xValues  = x.getNonZeroValues();
    const HArray<ValueType>& yValues  = y.getNonZeroValues();

    // copy the non-zero indexes

    mNonZeroIndexes = x.getNonZeroIndexes();

    // binary operation on values array

    HArrayUtils::binaryOp( mNonZeroValues, xValues, yValues, op, getContextPtr() );
    
    mZeroValue = common::applyBinary( x.getZero(), op, y.getZero() );

    SCAI_LOG_INFO( logger, "result = " << *this )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::binaryOpSparse( const SparseVector<ValueType>& x, const common::BinaryOp op, const SparseVector<ValueType>& y )
{
    SCAI_REGION( "Vector.Sparse.binOpSparse" )

    SCAI_ASSERT_EQ_DEBUG( x.getDistribution(), y.getDistribution(), "serious space mismatch" )

    if ( getDistribution() != x.getDistribution() )
    {
        // there is no alias of this vector, neither with x nor with y

        allocate( x.getDistributionPtr() );
    }

    SCAI_LOG_INFO( logger, "binaryOpSparse, x = " << x << ", y = " << y )

    bool samePattern = x.getNonZeroIndexes().size() == y.getNonZeroIndexes().size();

    if ( samePattern )
    {
        samePattern = HArrayUtils::all( x.getNonZeroIndexes(), common::CompareOp::EQ, y.getNonZeroIndexes() );
    }

    if ( samePattern )
    {
        binaryOpSparseSamePattern( x, op, y );
    }
    else
    {
        binaryOpSparseNewPattern( x, op, y );
    }
}

template<typename ValueType>
void SparseVector<ValueType>::binaryOpSparseNewPattern( const SparseVector<ValueType>& x, const common::BinaryOp op, const SparseVector<ValueType>& y )
{
    SCAI_REGION( "Vector.Sparse.binOpSparseNewPat" )

    SCAI_LOG_INFO( logger, "binaryOpSparseNewPattern: this = x " << op << " y, with x = " << x << ", y = " << y );

    const HArray<IndexType>& xIndexes = x.getNonZeroIndexes();
    const HArray<ValueType>& xValues  = x.getNonZeroValues();
    ValueType xZero = x.getZero();

    const HArray<IndexType>& yIndexes = y.getNonZeroIndexes();
    const HArray<ValueType>& yValues  = y.getNonZeroValues();
    ValueType yZero = y.getZero();

    // binary operation on sparse vectors

    HArray<ValueType> resultValues;
    HArray<IndexType> resultIndexes;

    HArrayUtils::binaryOpSparse( resultIndexes, resultValues,
                                 xIndexes, xValues, xZero,
                                 yIndexes, yValues, yZero, op, getContextPtr() );
    
    // apply binary for zero values, but be careful:
    //   if operation is illegal here, it might be still valid if there are no zero elements

    if ( op == common::BinaryOp::DIVIDE && yZero == common::Constants::ZERO )
    {
        SCAI_ASSERT_EQ_ERROR( resultIndexes.size(), getDistribution().getLocalSize(), "Sparse vector in division has zero elements: " << y )
        mZeroValue = xZero;
    }
    else
    {
        mZeroValue = common::applyBinary( xZero, op, yZero);
    }

    mNonZeroIndexes.swap( resultIndexes );
    mNonZeroValues.swap( resultValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::vectorPlusVector( const ValueType& alpha, const Vector<ValueType>& x, 
                                                const ValueType& beta, const Vector<ValueType>& y )
{
    if ( x.getVectorKind() == VectorKind::SPARSE  )
    {
        const SparseVector<ValueType>& spX = static_cast<const SparseVector<ValueType>&>( x );

        if ( y.getVectorKind() == VectorKind::SPARSE )
        {
            const SparseVector<ValueType>& spY = static_cast<const SparseVector<ValueType>&>( y );
         
            vectorPlusVectorImpl( alpha, spX, beta, spY );
            return;
        }
    }

    // just get it running: use DenseVector as temporary

    SCAI_UNSUPPORTED( "SparseVector<" << common::TypeTraits<ValueType>::id() << ">::vectorPlusVector( " 
                       << alpha << " * x + " << beta << " * y ) uses temporary dense vector" )

    DenseVector<ValueType> tmp;
    tmp.vectorPlusVector( alpha, x, beta, y );
    assign( tmp );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::vectorPlusVectorImpl( 
    const ValueType alpha, const SparseVector<ValueType>& x, 
    const ValueType beta, const SparseVector<ValueType>& y )
{
    SCAI_REGION( "Vector.Sparse.VplusV" )

    SCAI_ASSERT_EQ_ERROR( x.getDistribution(), y.getDistribution(), "mismatch distribution" );

    if ( &x == this || &y == this )
    {
         SCAI_LOG_INFO( logger, "result = " << alpha << " * x + " << beta << " * y, alias result = x or result = y" )

         // alias of input and output array, needs a temporay for sparse vectors

         SparseVector<ValueType> tmp( this->getContextPtr() );
         tmp.vectorPlusVectorImpl( alpha, x, beta, y );
         swap( tmp );
         return;
    }

    // Now we can just call addSparse for the local vectors

    setDistributionPtr( x.getDistributionPtr() );

    SCAI_LOG_INFO( logger, "addSparse: " << alpha << " * x + " << beta << " * y, x = " << x << ", y = " << y )

    HArrayUtils::addSparse( mNonZeroIndexes, mNonZeroValues,
                            x.mNonZeroIndexes, x.mNonZeroValues, x.mZeroValue, alpha, 
                            y.mNonZeroIndexes, y.mNonZeroValues, y.mZeroValue, beta, this->getContextPtr() );

    mZeroValue = alpha * x.mZeroValue + beta * y.mZeroValue;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::vectorTimesVector( 
    const ValueType& alpha, 
    const Vector<ValueType>& x, 
    const Vector<ValueType>& y )
{
    // just get it running: use DenseVector as temporary

    SCAI_UNSUPPORTED( "SparseVector<" << common::TypeTraits<ValueType>::id() << ">::vectorTimesVector( " 
                       << alpha << " * x * y ) uses temporary dense vector" )

    DenseVector<ValueType> tmp;
    tmp.vectorTimesVector( alpha, x, y );
    this->assign( tmp );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::vectorPlusScalar( const ValueType& alpha, const Vector<ValueType>& x, const ValueType& beta )
{
    if ( beta == common::Constants::ZERO )
    {
        binaryOpScalar( x, alpha, common::BinaryOp::MULT, true );
        return;
    }

    if ( alpha == common::Constants::ZERO )
    {
        // be careful, we inherit the space of x
       
        allocate( x.getDistributionPtr() );
        setScalar( beta );
        return;
    }

    if ( alpha == common::Constants::ONE )
    {
        // be careful, we inherit the space of x
       
        binaryOpScalar( x, beta, common::BinaryOp::ADD, false );
        return;
    }

    if ( x.getVectorKind() != VectorKind::SPARSE )
    {
        SCAI_UNSUPPORTED( "sparseVector = alpha * denseVector + beta not supported, use denseResult as temporary" )
        DenseVector<ValueType> tmpResult;
        tmpResult.vectorPlusScalar( alpha, x, beta );
        assign( tmpResult );
        return;
    }

    const SparseVector<ValueType>& sparseX = static_cast<const SparseVector<ValueType>&>( x );

    if ( &x != this )
    {
        // allocation and copy of non-zero indexes only required if there is no alias

        allocate( x.getDistributionPtr() );
        mNonZeroValues.resize( sparseX.mNonZeroValues.size() );
        mNonZeroIndexes = sparseX.mNonZeroIndexes;
    }

    mZeroValue = alpha * sparseX.mZeroValue + beta;

    utilskernel::HArrayUtils::arrayPlusScalar( mNonZeroValues, alpha, sparseX.mNonZeroValues, beta, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType SparseVector<ValueType>::dotProduct( const Vector<ValueType>& other ) const
{
    SCAI_REGION( "Vector.Sparse.dotP" )

    SCAI_LOG_INFO( logger, "Calculating dot product: " << *this << " * " << other )

    SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(),
                          "dotProduct not supported for vectors with different distributions. "
                          << *this  << " x " << other )

    ValueType localDotProduct;

    if ( &other == this )
    {
        // dot product with this sparse vector

        localDotProduct = HArrayUtils::dotProduct( mNonZeroValues, mNonZeroValues );
    }
    else if ( mZeroValue == common::Constants::ZERO )
    {
        HArray<ValueType> otherNonZeroValues;  //  the values form other at my non-zero indexes

        other.gatherLocalValues( otherNonZeroValues, mNonZeroIndexes );

        // now build dotproduct( mNonZeroValues, otherNonZeroValues )

        localDotProduct = HArrayUtils::dotProduct( mNonZeroValues, otherNonZeroValues );
    }
    else 
    {
         HArray<ValueType> multValues;

         buildLocalValues( multValues, common::BinaryOp::COPY, getContextPtr() );
         other.buildLocalValues( multValues, common::BinaryOp::MULT, getContextPtr() );
         localDotProduct = HArrayUtils::sum( multValues );
    }

    SCAI_LOG_DEBUG( logger, "Calculating global dot product form local dot product = " << localDotProduct )

    ValueType dotProduct = getDistribution().getCommunicator().sum( localDotProduct );

    SCAI_LOG_DEBUG( logger, "Global dot product = " << dotProduct )

    return dotProduct;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::selectComplexPart( Vector<RealType<ValueType> >& x, const common::ComplexPart part ) const
{
    HArray<RealType<ValueType>> localX;

    HArrayUtils::selectComplexPart( localX, getNonZeroValues(), part );

    RealType<ValueType> zero = part == common::ComplexPart::REAL ? common::Math::real( mZeroValue ) : common::Math::imag( mZeroValue );

    // now call virtual methods to set the sparse data to any vector

    x.setSameValue( getDistributionPtr(), zero );
    x.fillSparseData( getNonZeroIndexes(), localX, common::BinaryOp::COPY );
}

// template specializaton needed for IndexType, as Math::real and Math::imag are not supported for it

template<>
void SparseVector<IndexType>::selectComplexPart( Vector<IndexType>& x, const common::ComplexPart part ) const
{
    if ( part == common::ComplexPart::REAL )
    {
        x = *this;
    }
    else
    {
        x.setSameValue( getDistributionPtr(), IndexType( 0 ) );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::buildComplex( const Vector<RealType<ValueType> >& x, const Vector<RealType<ValueType> >& y )
{
    // ToDo: provide more efficient solutions as this one with two temporary complex vectors

    DenseVector<ValueType> x1;
    x1.assign( x );
    DenseVector<ValueType> y1;
    y1.assign( y );
    ValueType i = common::TypeTraits<ValueType>::imaginaryUnit();
    vectorPlusVector( 1, x1, i, y1 );
}

template<>
void SparseVector<IndexType>::buildComplex( const Vector<IndexType>& x, const Vector<IndexType>& )
{
    *this = x;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setVector( const _Vector& other, common::BinaryOp op, const bool swapArgs )
{
    SCAI_REGION( "Vector.Sparse.setVector" )

    SCAI_LOG_INFO( logger, "set " << *this << " with " << other << ", op = " << op )

    // op == COPY is special case 

    if ( op == common::BinaryOp::COPY )
    {
        SCAI_REGION( "Vector.Sparse.setVectorCopy" )

        SCAI_ASSERT_ERROR( !swapArgs, "swapping for binary COPY operator not allowed" )

        if ( &other == this )
        {
            return;
        }

        allocate( other.getDistributionPtr() );
        assign( other );

        return;
    }

    SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(), "setVector only with same distributions supported" );

    SCAI_ASSERT_ERROR( !swapArgs, "swapping arguments not supported yet" )

    if ( mZeroValue == ValueType( 0 ) && op == common::BinaryOp::MULT )
    {
        SCAI_REGION( "Vector.Sparse.setVectorMult0" )

        SCAI_LOG_INFO( logger, "setVector: thisSparse( zero = 0 ) *= otherSparse, other = " << other )

        // gather the values from other vector at the non-zero positions 

        HArray<ValueType> otherValues;  // = other[ nonZeroIndexes ]

        other.gatherLocalValues( otherValues, mNonZeroIndexes, common::BinaryOp::COPY, getContextPtr() );

        HArrayUtils::binaryOp( mNonZeroValues, mNonZeroValues, otherValues, op, getContextPtr() );
    }
    else if ( other.getValueType() == getValueType() )
    {
        SCAI_REGION( "Vector.Sparse.setVectorBinary" )

        const Vector<ValueType>& otherTyped = static_cast<const Vector<ValueType>&>( other );

        if ( !swapArgs )
        {
            binaryOp( *this, op, otherTyped );
        }
        else
        {
            binaryOp( otherTyped, op, *this );
        }
    }
    else
    {
        SCAI_REGION( "Vector.Sparse.setVectorConvert" )

        // Maybe not very efficient if other vector is dense

        SparseVector<ValueType> tmpOther;
        tmpOther.assign( other );

        SCAI_LOG_DEBUG( logger, "binary op, converted other = " << other << " -> " << tmpOther )

        if ( !swapArgs )
        {
            binaryOpSparse( *this, op, tmpOther );
        }
        else
        {
            binaryOpSparse( tmpOther, op, *this );
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::allocate( DistributionPtr distribution )
{
    setDistributionPtr( distribution );
    mNonZeroValues.clear();
    mNonZeroIndexes.clear();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::allocate( const IndexType n )
{
    setDistributionPtr( DistributionPtr( new NoDistribution( n ) ) );

    mNonZeroValues.clear();
    mNonZeroIndexes.clear();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setScalar( const ValueType& alpha )
{
    SCAI_LOG_INFO( logger, *this << ": set scalar " << alpha )

    // just set this values as the ZERO value for the sparse vector

    mNonZeroIndexes.clear();
    mNonZeroValues.clear();
    mZeroValue = alpha;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::prefetch( const ContextPtr location ) const
{
    mNonZeroValues.prefetch( location );
    mNonZeroIndexes.prefetch( location );
}

template<typename ValueType>
void SparseVector<ValueType>::wait() const
{
    mNonZeroValues.wait();
    mNonZeroIndexes.wait();
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
size_t SparseVector<ValueType>::getMemoryUsage() const
{
    // Note: memory of mHaloValues is not counted, is just a temporary

    IndexType localSize = mNonZeroValues.size();

    // Note: do sum with IndexType, as size_t is not yet handled by TypeTraits

    IndexType globalSize = getDistribution().getCommunicator().sum( localSize );

    // for each non zero value we have one index and one value

    return ( sizeof( ValueType ) + sizeof( IndexType ) ) * globalSize;
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void SparseVector<ValueType>::redistribute( DistributionPtr distribution )
{
    SCAI_LOG_INFO( logger, *this << ", redistribute to dist = " << *distribution )

    SCAI_ASSERT_EQ_ERROR( size(), distribution->getGlobalSize(), "global size mismatch between old/new distribution" )

    if ( getDistribution() == *distribution )
    {
        SCAI_LOG_INFO( logger, *this << " redistribute to same distribution " << *distribution )
        // we can keep local/global values, but just set dist pointer
        setDistributionPtr( distribution );
    }
    else if ( getDistribution().isReplicated() )
    {
        // each processor has all values, so just pick up the local values

        SCAI_LOG_DEBUG( logger, "localize replicated vector" )

        // we just compress the non-zero indexes/values owned by this process

        IndexType oldSize = mNonZeroIndexes.size();

        ContextPtr hostContext = Context::getHostPtr();

        {
            IndexType newSize = 0;

            WriteAccess<IndexType> wNonZeroIndexes( mNonZeroIndexes, hostContext );
            WriteAccess<ValueType> wNonZeroValues( mNonZeroValues, hostContext );

            for ( IndexType i = 0; i < oldSize; ++i )
            {
                IndexType globalIndex = wNonZeroIndexes[i];

                if ( distribution->isLocal( globalIndex ) )
                {
                    const IndexType localIndex = distribution->global2Local( globalIndex );
                    wNonZeroIndexes[newSize] = localIndex;
                    wNonZeroValues[newSize] = wNonZeroValues[i];
                    newSize++;
                }
            }
            
            wNonZeroIndexes.resize( newSize );
            wNonZeroValues.resize( newSize );
        
        }

        SCAI_LOG_DEBUG( logger, "Kept locally " << mNonZeroIndexes.size() << " of " << oldSize << " non-zero values" )

        setDistributionPtr( distribution );

        SCAI_ASSERT( HArrayUtils::validIndexes( mNonZeroIndexes, distribution->getLocalSize(), getContextPtr() ), "serious" )
    }
    else if ( distribution->isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, "replicate distributed sparse vector" )

        // translate the local non-zero indexes to 'global' non-zero indexes

        IndexType nLocalIndexes = mNonZeroIndexes.size();

        const Distribution& currentDist = getDistribution();

        currentDist.local2GlobalV( mNonZeroIndexes, mNonZeroIndexes );

        SCAI_LOG_DEBUG( logger, "translated " << nLocalIndexes << " local indexes to global indexes" )

        // ToDo: actually this is merge of sorted arrays

        HArray<IndexType> allNonZeroIndexes;
        HArray<ValueType> allNonZeroValues;

        getDistribution().getCommunicator().joinArray( allNonZeroIndexes, mNonZeroIndexes );
        getDistribution().getCommunicator().joinArray( allNonZeroValues, mNonZeroValues );

        // sort the non-zero indexes ascending

        HArrayUtils::sortSparseEntries( allNonZeroIndexes, allNonZeroValues, true, getContextPtr() );

        mNonZeroIndexes.swap( allNonZeroIndexes );
        mNonZeroValues.swap( allNonZeroValues );

        setDistributionPtr( distribution );

        SCAI_LOG_DEBUG( logger, "Here is the replicated sparse vector: " << *this )
    }
    else
    {
        SCAI_ASSERT_EQ_ERROR( distribution->getCommunicator(), getDistribution().getCommunicator(),
                              "redistribute only supported within same communicator" )

        SCAI_LOG_INFO( logger, *this << " will be redistributed to " << *distribution << " in two steps: replicate/localize" )

        auto repDist = std::make_shared<NoDistribution>( getDistribution().getGlobalSize() );

        redistribute( repDist );
        redistribute( distribution );

        // optimized pattern : shift all parts between all processors and pick up the new local ones
    }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void SparseVector<ValueType>::redistribute( const RedistributePlan& redistributor )
{
    // use a temporary dense vector for redistribution

    SCAI_LOG_WARN( logger, "use dense vector for redistribution" )

    auto tmp = convert<DenseVector<ValueType>>( *this );
    tmp.redistribute( redistributor );
    assign( tmp );
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void SparseVector<ValueType>::resize( DistributionPtr distribution )
{
    // disassemble this vector and fill it up again

    VectorAssembly<ValueType> assembly;

    this->disassemble( assembly );

    IndexType newSize = distribution->getGlobalSize();

    // truncate elements if necessary to avoid ERROR messages when we fil up

    if ( newSize < size() )
    {
        assembly.truncate( newSize );
    }

    // and now fill the assembly back

    allocate( distribution );

    this->fillFromAssembly( assembly );
}

/* ========================================================================= */
/*       Factory methods (must be provided to register)                      */
/* ========================================================================= */

template<typename ValueType>
_Vector* SparseVector<ValueType>::create()
{
    return new SparseVector<ValueType>();
}

template<typename ValueType>
VectorCreateKeyType SparseVector<ValueType>::createValue()
{
    return VectorCreateKeyType( VectorKind::SPARSE, common::getScalarType<ValueType>() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const SparseVector<ValueType>& other )

    : Vector<ValueType>( other ),
      mNonZeroIndexes( other.mNonZeroIndexes ),
      mNonZeroValues( other.mNonZeroValues ),
      mZeroValue( other.mZeroValue )

{
    // implementation here can be simpler as SparseVector( const _Vector& other )

    SCAI_LOG_INFO( logger,
                   "CopyConstructor of SparseVector " << size() << ", local size " << getDistribution().getLocalSize() )
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( SparseVector, SCAI_ARRAY_TYPES_HOST )

#define SPARSE_VECTOR_INST_LVL2( ValueType, OtherValueType )                                    \
    template void SparseVector<ValueType>::assignImpl( const Vector<OtherValueType>& );         \
    template void SparseVector<ValueType>::setDenseValuesImpl( const HArray<OtherValueType>& );   

#define SPARSE_VECTOR_INST_LVL1( ValueType )                                                    \
    SCAI_COMMON_LOOP_LVL2( ValueType, SPARSE_VECTOR_INST_LVL2, SCAI_ARRAY_TYPES_HOST )

SCAI_COMMON_LOOP( SPARSE_VECTOR_INST_LVL1, SCAI_ARRAY_TYPES_HOST )

#undef SPARSE_VECTOR_INST_LVL2
#undef SPARSE_VECTOR_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
