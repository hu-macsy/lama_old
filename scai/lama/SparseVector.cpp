/**
 * @file SparseVector.cpp
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

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/common/BinaryOp.hpp>

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

/* ------------------------------------------------------------------------- */
/*  Implementation of constructors for SparseVector                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::SparseVector() :

    _SparseVector( 0 ),
    mNonZeroIndexes(),
    mNonZeroValues(),
    mZeroValue( 0 )
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const IndexType n ) :

    _SparseVector( n ),
    mNonZeroIndexes(),
    mNonZeroValues(),
    mZeroValue( 0 )
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( ContextPtr context ) :

    _SparseVector( 0, context ),
    mNonZeroIndexes( context ),
    mNonZeroValues( context ),
    mZeroValue( 0 )
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const IndexType n, ContextPtr context ) :

    _SparseVector( n, context ),
    mNonZeroIndexes( context ),
    mNonZeroValues( context ),
    mZeroValue( 0 )
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( DistributionPtr distribution ) :
    
    _SparseVector( distribution ), 
    mNonZeroIndexes(),
    mNonZeroValues(),
    mZeroValue( 0 )
{
    SCAI_LOG_INFO( logger, "Construct sparse vector, size = " << distribution->getGlobalSize()
                         << ", distribution = " << *distribution << ", all zero" )
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( DistributionPtr distribution, ContextPtr context ) : 

    _SparseVector( distribution, context ), 
    mNonZeroIndexes( context ), 
    mNonZeroValues( context ),
    mZeroValue( 0 )

{
    SCAI_LOG_INFO( logger, "Construct sparse vector on context = " << context << ", size = " << distribution->getGlobalSize()
                         << ", distribution = " << *distribution << ", all zero" )
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const IndexType size, const ValueType value, ContextPtr context ) :

    _SparseVector( size, context ),
    mZeroValue( value )
{
    SCAI_LOG_INFO( logger, "Construct sparse vector, size = " << size << ", ZERO =" << value )
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Vector& other ) : 

    _SparseVector( other )

{
    allocate( getDistributionPtr() );
    assign( other );
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Vector& other, DistributionPtr distribution )

    : _SparseVector( other )

{
    assign( other );
    redistribute( distribution );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::SparseVector( dmemo::DistributionPtr distribution, const hmemo::_HArray& localValues ) :

    _SparseVector( distribution )

{
    setDenseValues( localValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const hmemo::_HArray& localValues ) :

    _SparseVector( DistributionPtr( new NoDistribution( localValues.size() ) ) )

{
    setDenseValues( localValues );   // builds the sparse version
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::SparseVector( 
    dmemo::DistributionPtr distribution,
    const hmemo::HArray<IndexType>& indexes, 
    const hmemo::_HArray& values, 
    const Scalar zero ) :

    _SparseVector( distribution )

{
    setSparseValues( indexes, values, zero );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const std::string& filename ) : 

    _SparseVector( 0 )

{
    SCAI_LOG_INFO( logger, "Construct sparse vector from file " << filename )
    readFromFile( filename );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setRandom( dmemo::DistributionPtr distribution, const float fillRate )
{
    allocate( distribution );

    // build dense random array and compress it

    LArray<ValueType> localValues;

    const IndexType localSize = getDistribution().getLocalSize();

    localValues.setRandom( localSize, fillRate, getContextPtr() );

    setDenseValues( localValues );
}

/* ------------------------------------------------------------------------- */

/*
 * Constructors with Expressions as arguments
 */

// linear algebra expression: a*x
template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_SV& expression )

    : _SparseVector( expression.getArg2() )

{
    SCAI_LOG_INFO( logger, "Constructor( alpha * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a+x/x+a
template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_SV_S& expression )

    : _SparseVector( expression.getArg1().getArg2() )
{
    SCAI_LOG_INFO( logger, "Constructor( alpha * x + beta)" )
    Vector::operator=( expression );
}

// linear algebra expression: x*y
template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_VV& expression )

    : _SparseVector( expression.getArg1() )

{
    Vector::operator=( expression );
}

// linear algebra expression: s*x*y
template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_SVV& expression ) :

    _SparseVector( expression.getArg2().getArg1() )

{
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x+b*y, inherit distribution/context from vector x

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_SV_SV& expression ) :

    _SparseVector( expression.getArg1().getArg2() )

{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x + beta * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_SMV_SV& expression )

    : _SparseVector( expression.getArg1().getArg2().getArg1().getRowDistributionPtr(),
              expression.getArg1().getArg2().getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_SVM_SV& expression )
    : _SparseVector( expression.getArg1().getArg2().getArg2().getColDistributionPtr(),
              expression.getArg1().getArg2().getArg2().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x, inherit distribution/context from matrix A

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_SMV& expression )

    : _SparseVector( expression.getArg2().getArg1().getRowDistributionPtr(),
              expression.getArg2().getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x*A, inherit distribution/context from matrix A

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Expression_SVM& expression )
    : _SparseVector( expression.getArg2().getArg2().getColDistributionPtr(),
              expression.getArg2().getArg2().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A )" )
    Vector::operator=( expression );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseVector<ValueType>::~SparseVector()
{
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

bool _SparseVector::isConsistent() const
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

    if ( !HArrayUtils::isSorted( nonZeroIndexes, common::binary::LT ) )
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
common::scalar::ScalarType SparseVector<ValueType>::getValueType() const
{
    return TypeTraits<ValueType>::stype;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::buildLocalValues( 
    _HArray& values, 
    const common::binary::BinaryOp op,
    ContextPtr loc ) const
{
    SCAI_LOG_INFO( logger, *this << ": build local values, op = " << op << " into " << values )

    const bool UNIQUE = true;

    // size of values will be local size of vector

    const IndexType size = getDistribution().getLocalSize();

    // convert the local sparse data to local dense data

    if ( op == common::binary::COPY )
    {
        // build values array from scratch 

        values.clear();
        values.resize( size );
        HArrayUtils::assignScalar( values, mZeroValue, op, getContextPtr() );
        HArrayUtils::scatter( values, mNonZeroIndexes, UNIQUE, mNonZeroValues, op, getContextPtr() );
    }
    else if ( mZeroValue == common::zeroBinary<ValueType>( op ) ) 
    {
        // ZERO element of this sparse vector is ZERO element for op, that is fine, we only apply non-zero values

        SCAI_ASSERT_EQ_ERROR( values.size(), size, "size mismatch" )
        HArrayUtils::scatter( values, mNonZeroIndexes, UNIQUE, mNonZeroValues, op, loc );
    }
    else 
    {
        // temporary array needed for this operation

        SCAI_LOG_WARN( logger, *this << ", is not ZERO element of " << op << ", temporary dense values are built" )

        utilskernel::LArray<ValueType> myDenseValues( size, mZeroValue );
        HArrayUtils::scatterImpl( myDenseValues, mNonZeroIndexes, UNIQUE, mNonZeroValues, common::binary::COPY, loc );
        HArrayUtils::setArray( values, myDenseValues, op, loc );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::gatherLocalValues(
    _HArray& values,
    const HArray<IndexType>& indexes,
    const common::binary::BinaryOp op,
    ContextPtr loc ) const
{
    HArrayUtils::sparseGather( values, mNonZeroValues, mNonZeroIndexes, indexes, op, loc );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setDenseValues( const _HArray& values )
{
    const IndexType size = getDistribution().getLocalSize();

    SCAI_ASSERT_EQ_ERROR( size, values.size(), "size of local values does not match local size of vector" )

    mZeroValue = ValueType( 0 );

    // ToDo: use the current ZERO value for building the sparse data structures

    HArrayUtils::buildSparseArray( mNonZeroValues, mNonZeroIndexes, values, getContextPtr() );
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

    SCAI_ASSERT_DEBUG( HArrayUtils::isSorted( mNonZeroIndexes, common::binary::LT ), "sort sparse entries failed" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setSparseValues( const HArray<IndexType>& nonZeroIndexes, const _HArray& nonZeroValues, const Scalar zeroValue )
{
    SCAI_LOG_INFO( logger, "setSparseValues, nnz = " << nonZeroIndexes.size() << ", ZERO = " << zeroValue )

    const IndexType size = getDistribution().getLocalSize();

    bool isValid = HArrayUtils::validIndexes( nonZeroIndexes, size, getContextPtr() );

    if ( !isValid )
    {
        COMMON_THROWEXCEPTION( "at least one illegal index, local size = " << size )
    }

    HArrayUtils::setArray( mNonZeroIndexes, nonZeroIndexes, common::binary::COPY, getContextPtr() );
    HArrayUtils::setArray( mNonZeroValues, nonZeroValues, common::binary::COPY, getContextPtr() );

    HArrayUtils::sortSparseEntries( mNonZeroIndexes, mNonZeroValues, true, getContextPtr() );

    mZeroValue = zeroValue.getValue<ValueType>();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType SparseVector<ValueType>::readLocalFromFile( const std::string& fileName, const IndexType first, const IndexType n )
{
    SCAI_REGION( "Vector.sparse.readLocal" )

    SCAI_LOG_INFO( logger, "read local array from file " << fileName )

    IndexType localN;   // for local size of the array data

    FileIO::read( localN, mNonZeroIndexes, mNonZeroValues, fileName );

    HArrayUtils::sortSparseEntries( mNonZeroIndexes, mNonZeroValues, true, getContextPtr() );

    // ToDo: read block from sparse array

    SCAI_ASSERT_EQ_ERROR( 0, first, "block read not supported for sparse data" )
    SCAI_ASSERT_EQ_ERROR( nIndex, n, "block read not supported for sparse data" )

    return localN;
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
    common::unique_ptr<SparseVector<ValueType> > vector( new SparseVector<ValueType>() );
    vector->setContextPtr( this->getContextPtr() );
    return vector.release();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::getValue( IndexType globalIndex ) const
{
    ValueType myValue = 0;

    const IndexType localIndex = getDistribution().global2local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": getValue( globalIndex = " << globalIndex << " ) -> local : " << localIndex )

    if ( localIndex != nIndex )
    {
        // we have here sparse data, so look for local index among the sparse indexes of non-zero values

        IndexType pos = HArrayUtils::findPosInSortedIndexes( mNonZeroIndexes, localIndex );

        if ( pos != nIndex )
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

    return Scalar( allValue );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setValue( const IndexType globalIndex, const Scalar value )
{
    SCAI_ASSERT_VALID_INDEX_ERROR( globalIndex, size(), "illegal index" )

    SCAI_LOG_TRACE( logger, *this << ": setValue( globalIndex = " << globalIndex << " ) = " <<  value )

    const IndexType localIndex = getDistribution().global2local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": set @g " << globalIndex << " is @l " << localIndex << " : " << value )

    if ( localIndex != nIndex )
    {
        // This partition is the owner, add it locally

        IndexType pos = HArrayUtils::findPosInSortedIndexes( mNonZeroIndexes, localIndex );

        ValueType typedValue = value.getValue<ValueType>();

        if ( pos != nIndex )
        {
            mNonZeroValues[pos] = typedValue;
        }
        else
        {
            // add a new entry in mNonZeroIndexes, mNonZeroValues

            pos = HArrayUtils::insertSorted( mNonZeroIndexes, localIndex );
            SCAI_LOG_TRACE( logger, "setValue, local index = " << localIndex << " at pos = " << pos )
            HArrayUtils::insertAtPos( mNonZeroValues, pos, typedValue );
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::min() const
{
    // Note: min returns the maximal representation value on zero-sized vectors, TypeTraits<ValueType>::getMax()

    ValueType localMin = mNonZeroValues.min();

    // if there are implicit zero values they must be used for min computation

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
         localMin = common::Math::min( localMin, mZeroValue );
    }

    return Scalar( getDistribution().getCommunicator().min( localMin ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::max() const
{
    // Note: max returns the minimal representation value on zero-sized vectors

    ValueType localMax = mNonZeroValues.max();

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        // ToDo: replace ABS with ASUM, is different for complex numbers

        localMax = common::Math::max( localMax, mZeroValue );
    }

    return Scalar( getDistribution().getCommunicator().max( localMax ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::l1Norm() const
{
    SCAI_REGION( "Vector.sparse.l1Norm" )

    ValueType localL1Norm = mNonZeroValues.l1Norm();

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        // ToDo: replace ABS with ASUM, is different for complex numbers

        localL1Norm += common::applyUnary( common::unary::ABS, mZeroValue ) * nZero;
    }

    return Scalar( getDistribution().getCommunicator().sum( localL1Norm ) );
}

/*---------------------------------------------------------------------------*/
template<typename ValueType>
Scalar SparseVector<ValueType>::sum() const
{
    ValueType localSum = mNonZeroValues.sum();

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        localSum += mZeroValue * nZero;
    }

    return Scalar( getDistribution().getCommunicator().sum( localSum ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::l2Norm() const
{
    SCAI_REGION( "Vector.sparse.l2Norm" )

    // Note: we do not call l2Norm here for mNonZeroValues to avoid sqrt

    ValueType localDotProduct = mNonZeroValues.dotProduct( mNonZeroValues );

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        localDotProduct += mZeroValue * mZeroValue * nZero;
    }
 
    ValueType globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );

    return Scalar( common::Math::sqrt( globalDotProduct ) );
}

/* ------------------------------------------------------------------------- */

template<>
Scalar SparseVector<IndexType>::l2Norm() const
{
    SCAI_REGION( "Vector.sparse.l2Norm" )

    // Note: we do not call l2Norm here for mNonZeroValues to avoid sqrt

    ScalarRepType localDotProduct = mNonZeroValues.dotProduct( mNonZeroValues );
    ScalarRepType globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );
    return Scalar( common::Math::sqrt( globalDotProduct ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::maxNorm() const
{
    SCAI_REGION( "Vector.sparse.maxNorm" )

    ValueType localMaxNorm = mNonZeroValues.maxNorm();

    // the ZERO element must also be considered if at least one element is zero

    IndexType nZero = getDistribution().getLocalSize() - mNonZeroValues.size();

    if ( nZero > 0 )
    {
        SCAI_LOG_DEBUG( logger, "maxNorm, zero = " << mZeroValue << ", non-zero = " << localMaxNorm )
        localMaxNorm = common::Math::max( common::Math::abs( localMaxNorm ), common::Math::abs( mZeroValue ) );
    }

    const Communicator& comm = getDistribution().getCommunicator();
    ValueType globalMaxNorm = comm.max( localMaxNorm );
    SCAI_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm: " << localMaxNorm
                   << ", max norm global = " << globalMaxNorm )
    return Scalar( globalMaxNorm );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::maxDiffNorm( const Vector& other ) const
{
    // ToDo: find some more efficient solutions wherever possible

    SparseVector<ValueType> tmp( *this );
    tmp.setVector( other, common::binary::SUB, false );
    return tmp.maxNorm();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::swap( Vector& other )
{
    SCAI_LOG_DEBUG( logger, "swap:" << *this << " with " << other )

    SCAI_ASSERT_EQ_ERROR( getVectorKind(), other.getVectorKind(), "Swap only for same kind of vector allowed" )
    SCAI_ASSERT_EQ_ERROR( getValueType(), other.getValueType(), "Swap only for same value type of vector allowed" )

    SparseVector<ValueType>& typedOther = reinterpret_cast<SparseVector<ValueType>&>( other );

    Vector::swapVector( other );   // swap sizes, distributions

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
void SparseVector<ValueType>::vectorPlusVector( const Scalar& alpha, const Vector& x, const Scalar& beta, const Vector& y )
{
    if ( x.getValueType() == getValueType() && x.getVectorKind() == Vector::SPARSE  )
    {
        const SparseVector<ValueType>& spX = reinterpret_cast<const SparseVector<ValueType>&>( x );

        if ( y.getValueType() == getValueType() && y.getVectorKind() == Vector::SPARSE )
        {
            const SparseVector<ValueType>& spY = reinterpret_cast<const SparseVector<ValueType>&>( y );
         
            const ValueType alphaV = alpha.getValue<ValueType>();
            const ValueType betaV = beta.getValue<ValueType>();
 
            vectorPlusVectorImpl( alphaV, spX, betaV, spY );
            return;
        }
    }

    // just get it running: use DenseVector as temporary

    SCAI_LOG_WARN( logger, "SparseVector<" << common::TypeTraits<ValueType>::id() << ">::vectorPlusVector( " 
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
void SparseVector<ValueType>::vectorTimesVector( const Scalar& alpha, const Vector& x, const Vector& y )
{
    // just get it running: use DenseVector as temporary

    SCAI_LOG_WARN( logger, "SparseVector<" << common::TypeTraits<ValueType>::id() << ">::vectorTimesVector( " 
                          << alpha << " * x * y ) uses temporary dense vector" )

    DenseVector<ValueType> tmp;
    tmp.vectorTimesVector( alpha, x, y );
    assign( tmp );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::vectorPlusScalar( const Scalar& alpha, const Vector& x, const Scalar& beta )
{
    // just get it running: use DenseVector as temporary

    SCAI_LOG_WARN( logger, "SparseVector<" << common::TypeTraits<ValueType>::id() << ">::vectorAddScalar( " 
                          << alpha << " * x + " << beta << " ) uses temporary dense vector" )

    DenseVector<ValueType> tmp;
    tmp.vectorPlusScalar( alpha, x, beta );
    assign( tmp );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::dotProduct( const Vector& other ) const
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

        localDotProduct = mNonZeroValues.dotProduct( mNonZeroValues );
    }
    else if ( mZeroValue == common::constants::ZERO )
    {
        HArray<ValueType> otherNonZeroValues;  //  the values form other at my non-zero indexes

        other.gatherLocalValues( otherNonZeroValues, mNonZeroIndexes );

        // now build dotproduct( mNonZeroValues, otherNonZeroValues )

        localDotProduct = mNonZeroValues.dotProduct( otherNonZeroValues );
    }
    else 
    {
         utilskernel::LArray<ValueType> multValues;

         buildLocalValues( multValues, common::binary::COPY, mContext );
         other.buildLocalValues( multValues, common::binary::MULT, mContext );
         localDotProduct = multValues.sum();
    }

    SCAI_LOG_DEBUG( logger, "Calculating global dot product form local dot product = " << localDotProduct )

    ValueType dotProduct = getDistribution().getCommunicator().sum( localDotProduct );

    SCAI_LOG_DEBUG( logger, "Global dot product = " << dotProduct )

    return Scalar( dotProduct );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::binOpSparse( const _SparseVector& other,
                                           const common::binary::BinaryOp op,
                                           bool swapArgs )
{
    SCAI_LOG_INFO( logger, *this << " " << op << " = " << other )

    SCAI_ASSERT_EQ_ERROR( other.getValueType(), getValueType(), "type mismatch not handled yet" );

    const HArray<IndexType>& otherIndexes = other.getNonZeroIndexes();
    const HArray<ValueType>& otherValues = reinterpret_cast<const HArray<ValueType>&>( other.getNonZeroValues() );

    ValueType otherZero = other.getZero().getValue<ValueType>();

    // binary operation on sparse vectors

    HArray<ValueType> resultValues;
    HArray<IndexType> resultIndexes;

    if ( !swapArgs )
    {
        HArrayUtils::binaryOpSparse( resultIndexes, resultValues,
                                     mNonZeroIndexes, mNonZeroValues, mZeroValue,
                                     otherIndexes, otherValues, otherZero, op, mContext );
    
        mZeroValue = common::applyBinary( mZeroValue, op, otherZero );
    }
    else
    {
        HArrayUtils::binaryOpSparse( resultIndexes, resultValues,
                                     otherIndexes, otherValues, otherZero, 
                                     mNonZeroIndexes, mNonZeroValues, mZeroValue,
                                     op, mContext );
    
        mZeroValue = common::applyBinary( otherZero, op, mZeroValue );
    }

    // ToDo: remove entries in non-zero values that are now ZERO 

    swapSparseValues( resultIndexes, resultValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setVector( const Vector& other, common::binary::BinaryOp op, const bool swapArgs )
{
    SCAI_REGION( "Vector.Sparse.setVector" )

    SCAI_LOG_INFO( logger, "set " << *this << " with " << other << ", op = " << op )

    SCAI_ASSERT_EQ_ERROR( getDistribution(), other.getDistribution(), "setVector only with same distributions supported" );

    SCAI_ASSERT_ERROR( !swapArgs, "swapping arguments not supported yet" )

    if ( mZeroValue == ValueType( 0 ) && op == common::binary::MULT )
    {
        // gather the values from other vector at the non-zero positions 

        HArray<ValueType> otherValues;  // = other[ nonZeroIndexes ]
        other.gatherLocalValues( otherValues, mNonZeroIndexes, common::binary::COPY, getContextPtr() );
        HArrayUtils::binaryOp( mNonZeroValues, mNonZeroValues, otherValues, op, getContextPtr() );
    }
    else if ( other.getVectorKind() == Vector::SPARSE )
    {
        const _SparseVector& otherSparse = reinterpret_cast<const _SparseVector&>( other );

        binOpSparse( otherSparse, op, swapArgs );
    }
    else
    {
        COMMON_THROWEXCEPTION( "setVector with dense vector makes this sparse vector dense, unsupported" )
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
void SparseVector<ValueType>::assign( const Scalar value )
{
    SCAI_LOG_INFO( logger, *this << ": assign " << value )

    // just set this values as the ZERO value for the sparse vector

    mNonZeroIndexes.clear();
    mNonZeroValues.clear();
    mZeroValue = value.getValue<ValueType>();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setScalar( const Scalar value, common::binary::BinaryOp op, const bool swapScalar )
{
    SCAI_LOG_DEBUG( logger, *this << ": set " << value << ", op = " << op )

    ValueType val = value.getValue<ValueType>();

    if ( op == common::binary::COPY )
    {
        mNonZeroIndexes.clear();
        mNonZeroValues.clear();
        mZeroValue = val;
        return;
    }

    if ( swapScalar )
    {
        mZeroValue = common::applyBinary( val, op, mZeroValue );
    }
    else
    {
        mZeroValue = common::applyBinary( mZeroValue, op, val );
    }

    HArrayUtils::binaryOpScalar( mNonZeroValues, mNonZeroValues, val, op, swapScalar, mContext );
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
void SparseVector<ValueType>::applyUnary( common::unary::UnaryOp op )
{
    mZeroValue = common::applyUnary( op, mZeroValue );
    HArrayUtils::unaryOp( mNonZeroValues, mNonZeroValues, op, mContext );
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
                    const IndexType localIndex = distribution->global2local( globalIndex );
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

        ContextPtr hostContext = Context::getHostPtr();
        {
            WriteAccess<IndexType> wNonZeroIndexes( mNonZeroIndexes, hostContext );

            for ( IndexType i = 0; i < nLocalIndexes; ++i )
            {
                wNonZeroIndexes[i] = currentDist.local2global( wNonZeroIndexes[i] );
            }
        }

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
        SCAI_LOG_INFO( logger, *this << " will be redistributed to " << *distribution << " in two steps: replicate/localize" )

        DistributionPtr repDist ( new NoDistribution( getDistribution().getGlobalSize() ) );

        redistribute( repDist );
        redistribute( distribution );

        // optimized pattern : shift all parts between all processors and pick up the new local ones
    }
}

/* -- IO ------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::writeLocalToFile(
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

        // write the sparse data

        const IndexType size = getDistribution().getLocalSize();

        if ( mZeroValue == common::constants::ZERO )
        {
            fileIO->writeSparse( size, mNonZeroIndexes, mNonZeroValues, fileName );
        }
        else
        {
            // build a dense array on the host where it is used for the output

            hmemo::ContextPtr ctx = hmemo::Context::getHostPtr();
            LArray<ValueType> denseArray( size, mZeroValue, ctx );
            HArrayUtils::scatterImpl( denseArray, mNonZeroIndexes, true, mNonZeroValues, common::binary::COPY, ctx );
            fileIO->writeArray( denseArray, fileName );
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "File : " << fileName << ", unknown suffix" )
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Vector* SparseVector<ValueType>::create()
{
    return new SparseVector<ValueType>();
}

template<typename ValueType>
VectorCreateKeyType SparseVector<ValueType>::createValue()
{
    return VectorCreateKeyType( Vector::SPARSE, common::getScalarType<ValueType>() );
}

template<typename ValueType>
VectorCreateKeyType SparseVector<ValueType>::getCreateValue() const
{
    return createValue();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const SparseVector<ValueType>& other )

    : _SparseVector( other ),
      mNonZeroIndexes( other.mNonZeroIndexes ),
      mNonZeroValues( other.mNonZeroValues ),
      mZeroValue( other.mZeroValue )

{
    // implementation here can be simpler as SparseVector( const Vector& other )

    SCAI_LOG_INFO( logger,
                   "CopyConstructor of SparseVector " << size() << ", local size " << getDistribution().getLocalSize() )
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( SparseVector, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
