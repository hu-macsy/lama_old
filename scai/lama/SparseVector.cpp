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
 * @brief Implementations and instantiations for class SparseVector.
 * @author Jiri Kraus
 * @date 22.02.2011
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

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
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
/*  Implementation of methods/constructors for _SparseVector                 */
/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( _SparseVector::logger, "Vector.SparseVector" )

_SparseVector::_SparseVector( const IndexType n ) :

    Vector( n )
{
}

_SparseVector::_SparseVector( const IndexType n, ContextPtr context ) :

    Vector( n, context )
{
}

_SparseVector::_SparseVector( DistributionPtr dist ) :

    Vector( dist )
{
}

_SparseVector::_SparseVector( DistributionPtr dist, ContextPtr context ) :

    Vector( dist, context )
{
}

_SparseVector::_SparseVector( const _SparseVector& other ) :

    Vector( other )
{
}

_SparseVector::_SparseVector( const Vector& other ) :

    Vector( other )
{
}

_SparseVector* _SparseVector::create( common::scalar::ScalarType type )
{
    // There is only one factor for all vectors

    Vector* v = Vector::create( VectorCreateKeyType( Vector::SPARSE, type ) );

    // reinterpret cast is safe

    return reinterpret_cast<_SparseVector*>( v );
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector() :

    _SparseVector( 0 ),
    mNonZeroIndexes(),
    mNonZeroValues()
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const IndexType n ) :

    _SparseVector( n ),
    mNonZeroIndexes(),
    mNonZeroValues()
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( ContextPtr context ) :

    _SparseVector( 0, context ),
    mNonZeroIndexes( context ),
    mNonZeroValues( context )
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const IndexType n, ContextPtr context ) :

    _SparseVector( n, context ),
    mNonZeroIndexes( context ),
    mNonZeroValues( context )
{
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( DistributionPtr distribution ) :
    
    _SparseVector( distribution ), 
    mNonZeroIndexes(),
    mNonZeroValues()
{
    SCAI_LOG_INFO( logger, "Construct sparse vector, size = " << distribution->getGlobalSize()
                   // << ", type = " << typename(ValueType)
                   << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", no initialization" )
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( DistributionPtr distribution, ContextPtr context ) : 

    _SparseVector( distribution, context ), 
    mNonZeroIndexes( context ), 
    mNonZeroValues( context )

{
    SCAI_LOG_INFO( logger, "Construct sparse vector on context = " << context << ", size = " << distribution->getGlobalSize()
                   // << ", type = " << typename(ValueType)
                   << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", no initialization" )
}

template<typename ValueType>
SparseVector<ValueType>::SparseVector( const Vector& other )

    : _SparseVector( other )

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
SparseVector<ValueType>::SparseVector( const hmemo::_HArray& localValues, dmemo::DistributionPtr distribution ) :

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
SparseVector<ValueType>::SparseVector( const std::string& filename ) 

    : _SparseVector( 0 )
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
    const utilskernel::binary::BinaryOp op,
    ContextPtr loc ) const
{
    // size of values will be local size of vector

    const IndexType size = getDistribution().getLocalSize();

    // convert the local sparse data to local dense data

    if ( op == utilskernel::binary::COPY )
    {
        // values might be uninitialized

        HArrayUtils::buildDenseArray( values, size, mNonZeroValues, mNonZeroIndexes, loc );
    }
    else
    {
        SCAI_ASSERT_EQ_ERROR( values.size(), size, "size mismatch" )
        bool unique = true;   // no double entry in mNonZeroIndexes
        HArrayUtils::scatter( values, mNonZeroIndexes, unique, mNonZeroValues, op, loc );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setDenseValues( const _HArray& values )
{
    const IndexType size = getDistribution().getLocalSize();

    SCAI_ASSERT_EQ_ERROR( size, values.size(), "size of local values does not match local size of vector" )

    HArrayUtils::buildSparseArray( mNonZeroValues, mNonZeroIndexes, values, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::setSparseValues( const HArray<IndexType>& nonZeroIndexes, const _HArray& nonZeroValues )
{
    const IndexType size = getDistribution().getLocalSize();

    bool isValid = HArrayUtils::validIndexes( nonZeroIndexes, size, getContextPtr() );

    if ( !isValid )
    {
        COMMON_THROWEXCEPTION( "at least one illegal index, local size = " << size )
    }

    // we cannot check yet for double entries of one index, but for ascending order

    bool isSorted = HArrayUtils::isSorted( nonZeroIndexes, true );

    if ( !isSorted )
    {
        COMMON_THROWEXCEPTION( "indexes of non-zero values must be sorted" )
    }

    HArrayUtils::setArray( mNonZeroIndexes, nonZeroIndexes, utilskernel::binary::COPY, getContextPtr() );
    HArrayUtils::setArray( mNonZeroValues, nonZeroValues, utilskernel::binary::COPY, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType SparseVector<ValueType>::readLocalFromFile( const std::string& fileName, const IndexType first, const IndexType n )
{
    SCAI_LOG_INFO( logger, "read local array from file " << fileName )

    // sparse vector read not supported yet, so read a dense array

    LArray<ValueType> denseValues;

    FileIO::read( denseValues, fileName, common::scalar::INTERNAL, first, n );

    setDenseValues( denseValues );

    return denseValues.size();
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
    SCAI_LOG_TRACE( logger, *this << ": setValue( globalIndex = " << globalIndex << " ) = " <<  value )

    const IndexType localIndex = getDistribution().global2local( globalIndex );

    SCAI_LOG_TRACE( logger, *this << ": set @g " << globalIndex << " is @l " << localIndex << " : " << value )

    if ( localIndex != nIndex )
    {
        IndexType pos = HArrayUtils::findPosInSortedIndexes( mNonZeroIndexes, localIndex );

        if ( pos != nIndex )
        {
            mNonZeroValues[pos] = value.getValue<ValueType>();
        }
        else
        {
            COMMON_THROWEXCEPTION( "cannot set non-existing element in sparse vector, index = " << globalIndex )
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

    if ( mNonZeroValues.size() < getDistribution().getLocalSize() ) 
    {
         localMin = common::Math::min( localMin, ValueType( 0 ) );
    }

    return Scalar( getDistribution().getCommunicator().min( localMin ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::max() const
{
    // Note: max returns the minimal representation value on zero-sized vectors
    ValueType localMax = mNonZeroValues.max();
    return Scalar( getDistribution().getCommunicator().max( localMax ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::l1Norm() const
{
    SCAI_REGION( "Vector.sparse.l1Norm" )

    ValueType localL1Norm = mNonZeroValues.l1Norm();
    return Scalar( getDistribution().getCommunicator().sum( localL1Norm ) );
}

/*---------------------------------------------------------------------------*/
template<typename ValueType>
Scalar SparseVector<ValueType>::sum() const
{
    ValueType localsum = mNonZeroValues.sum();
    return Scalar( getDistribution().getCommunicator().sum( localsum ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseVector<ValueType>::l2Norm() const
{
    SCAI_REGION( "Vector.sparse.l2Norm" )

    // Note: we do not call l2Norm here for mNonZeroValues to avoid sqrt
    ValueType localDotProduct = mNonZeroValues.dotProduct( mNonZeroValues );
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
    const Communicator& comm = getDistribution().getCommunicator();
    ValueType globalMaxNorm = comm.max( localMaxNorm );
    SCAI_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm: " << localMaxNorm
                   << ", max norm global = " << globalMaxNorm )
    return Scalar( globalMaxNorm );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::swap( Vector& other )
{
    SCAI_LOG_DEBUG( logger, "swap:" << *this << " with " << other )
    SparseVector<ValueType>* otherPtr = dynamic_cast<SparseVector<ValueType>*>( &other );

    if ( !otherPtr )
    {
        COMMON_THROWEXCEPTION( "Tried to swap with a Vector of a different type." )
    }

    Vector::swapVector( other );
    mNonZeroValues.swap( otherPtr->mNonZeroValues );
    mNonZeroIndexes.swap( otherPtr->mNonZeroIndexes );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::writeAt( std::ostream& stream ) const
{
    const Distribution& dist = getDistribution();

    stream << "SparseVector<" << getValueType() << ">" << "( size = " << size() << ", local nnz = " << mNonZeroIndexes.size()
           << ", dist = " << dist << ", loc  = " << *getContextPtr() << " )";
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

    SCAI_ASSERT_NE_ERROR( &x, this, "alias result, x" )
    SCAI_ASSERT_NE_ERROR( &y, this, "alias result, y" )

    // Now we can just call addSparse for the local vectors

    setDistributionPtr( x.getDistributionPtr() );

    HArrayUtils::addSparse( mNonZeroIndexes, mNonZeroValues,
                            x.getNonZeroIndexes(), x.getNonZeroValues(), alpha, 
                            y.getNonZeroIndexes(), y.getNonZeroValues(), beta, this->getContextPtr() );
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

    HArray<ValueType> otherLocalValues;
    other.buildLocalValues( otherLocalValues );

    HArray<ValueType> otherNonZeroValues;   // get the values form other at my non-zero indexes

    utilskernel::HArrayUtils::gather( otherNonZeroValues, otherLocalValues, mNonZeroIndexes, utilskernel::binary::COPY );
 
    // now build dotproduct( mNonZeroValues, otherNonZeroValues )

    const ValueType localDotProduct = mNonZeroValues.dotProduct( otherNonZeroValues );

    SCAI_LOG_DEBUG( logger, "Calculating global dot product form local dot product = " << localDotProduct )

    ValueType dotProduct = getDistribution().getCommunicator().sum( localDotProduct );

    SCAI_LOG_DEBUG( logger, "Global dot product = " << dotProduct )

    return Scalar( dotProduct );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::scale( const Vector& other )
{
    SCAI_LOG_WARN( logger, "SparseVector<" << common::TypeTraits<ValueType>::id() << ">::scale( x )" 
                           << " uses temporary dense vector" )

    DenseVector<ValueType> tmp( *this );
    tmp.scale( other );
    assign( tmp );
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

    SCAI_UNSUPPORTED( "Assignment of scalar to a sparse vector not very useful" )
 
    const IndexType localSize = getDistribution().getLocalSize();

    HArrayUtils::setOrder( mNonZeroIndexes, localSize, getContextPtr() );

    mNonZeroValues.clear();
    mNonZeroValues.resize( localSize );

    HArrayUtils::assignScalar( mNonZeroValues, value.getValue<ValueType>(), utilskernel::binary::COPY, getContextPtr() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseVector<ValueType>::add( const Scalar value )
{
    SCAI_LOG_DEBUG( logger, *this << ": add " << value )

    COMMON_THROWEXCEPTION( "add scalar ( = " << value << " ) unsupported for sparse vector" )
}

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

template<typename ValueType>
void SparseVector<ValueType>::invert()
{
    // invert not very useful on a sparse vector that usually has a lot of zero values

    const IndexType localSize = getDistribution().getLocalSize();

    SCAI_ASSERT_EQ_ERROR( localSize, mNonZeroValues.size(), "invert on sparse vector with zero values" )
    SCAI_ASSERT_EQ_ERROR( localSize, mNonZeroIndexes.size(), "invert on sparse vector with zero values" )

    mNonZeroValues.invert();
}

template<typename ValueType>
void SparseVector<ValueType>::conj()
{
    mNonZeroValues.conj();
}

template<typename ValueType>
void SparseVector<ValueType>::exp()
{
    COMMON_THROWEXCEPTION( "log not supported for sparse vectors" )
}

template<typename ValueType>
void SparseVector<ValueType>::log()
{
    COMMON_THROWEXCEPTION( "log not supported for sparse vectors" )
}

template<typename ValueType>
void SparseVector<ValueType>::floor()
{
    mNonZeroValues.floor();
}

template<typename ValueType>
void SparseVector<ValueType>::ceil()
{
    mNonZeroValues.ceil();
}

template<typename ValueType>
void SparseVector<ValueType>::sqrt()
{
    mNonZeroValues.sqrt();
}

template<typename ValueType>
void SparseVector<ValueType>::sin()
{
    // supported for sparse vectors, as sin( 0 ) = 0
    mNonZeroValues.sin();
}

template<typename ValueType>
void SparseVector<ValueType>::cos()
{
    // unsupported for sparse vectors, as cos( 0 ) = 1

    COMMON_THROWEXCEPTION( "cos not supported for sparse vectors" )
}

template<typename ValueType>
void SparseVector<ValueType>::tan()
{
    COMMON_THROWEXCEPTION( "tan not supported for sparse vectors" )
}

template<typename ValueType>
void SparseVector<ValueType>::atan()
{
    COMMON_THROWEXCEPTION( "atan not supported for sparse vectors" )
}

template<typename ValueType>
void SparseVector<ValueType>::powExp( const Vector& other )
{
    COMMON_THROWEXCEPTION( "powExp not supported for sparse vectors, other = " << other )
}

template<typename ValueType>
void SparseVector<ValueType>::powBase( const Vector& other )
{
    COMMON_THROWEXCEPTION( "powBase not supported for sparse vectors, other = " << other )
}

template<typename ValueType>
void SparseVector<ValueType>::powBase( Scalar base )
{
    COMMON_THROWEXCEPTION( "powBase not supported for sparse vectors, base = " << base )
}

template<typename ValueType>
void SparseVector<ValueType>::powExp( Scalar exp )
{
    COMMON_THROWEXCEPTION( "powExp not supported for sparse vectors, exp = " << exp )
}

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

        IndexType globalSize = currentDist.getGlobalSize();

        HArray<IndexType> allNonZeroIndexes;
        HArray<ValueType> allNonZeroValues;

        {
            ReadAccess<IndexType> rNonZeroIndexes( mNonZeroIndexes, hostContext );
            ReadAccess<ValueType> rNonZeroValues( mNonZeroValues, hostContext );

            WriteOnlyAccess<IndexType> wAllNonZeroIndexes( allNonZeroIndexes, hostContext, globalSize );
            WriteOnlyAccess<ValueType> wAllNonZeroValues( allNonZeroValues, hostContext, globalSize );

            // Now replicate 

            currentDist.replicate( wAllNonZeroIndexes.get(), rNonZeroIndexes.get() );
            currentDist.replicate( wAllNonZeroValues.get(), rNonZeroValues.get() );
        }

        // sort the non-zero indexes ascending

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

        // build dense local vector data for output

        HArray<ValueType> denseLocalValues;

        buildLocalValues( denseLocalValues );

        fileIO->writeArray( denseLocalValues, fileName );
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
      mNonZeroValues( other.mNonZeroValues )

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
