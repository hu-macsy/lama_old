/**
 * @file DenseVector.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
#include <scai/lama/StorageIO.hpp>

#include <scai/lama/io/FileType.hpp>
#include <scai/lama/io/IOUtils.hpp>
#include <scai/lama/io/FileIO.hpp>

#include <scai/lama/mepr/IOWrapper.hpp>

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
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

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DenseVector<ValueType>::logger, "Vector.DenseVector" )

template<typename ValueType>
DenseVector<ValueType>::DenseVector() :
    Vector( 0 ),
    mLocalValues()
{
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( ContextPtr context ) :
    Vector( 0, context ),
    mLocalValues()
{
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution )
    : Vector( distribution ), mLocalValues( distribution->getLocalSize() )
{
    SCAI_LOG_INFO( logger, "Construct dense vector, size = " << distribution->getGlobalSize()
                   // << ", type = " << typename(ValueType)
                   << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", no initialization" )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, ContextPtr context )
    : Vector( distribution, context ), mLocalValues( distribution->getLocalSize() )
{
    SCAI_LOG_INFO( logger, "Construct dense vector on context = " << context << ", size = " << distribution->getGlobalSize()
                   // << ", type = " << typename(ValueType)
                   << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", no initialization" )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const IndexType size, const ValueType value, ContextPtr context )
    : Vector( size, context ), mLocalValues( size, value, context )
{
    SCAI_LOG_INFO( logger, "Construct dense vector, size = " << size << ", init =" << value )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( DistributionPtr distribution, const ValueType value, ContextPtr context )
    : Vector( distribution, context ), mLocalValues( distribution->getLocalSize(), value )
{
    SCAI_LOG_INFO( logger,
                   "Construct dense vector, size = " << distribution->getGlobalSize() << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", init = " << value )
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Vector& other )
    : Vector( other )
{
    allocate( getDistributionPtr() );
    assign( other );
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Vector& other, DistributionPtr distribution )

    : Vector( other )
{
    assign( other );
    redistribute( distribution );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const std::string& filename )
{
    SCAI_LOG_INFO( logger, "Construct dense vector from file " << filename )
    readFromFile( filename );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::readFromFile( const std::string& filename )
{
    SCAI_LOG_INFO( logger, "read dense vector from file " << filename )
    // Take the current default communicator
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    IndexType myRank = comm->getRank();
    IndexType host = 0; // reading processor

    if ( myRank == host )
    {
        // Only host reads the values

        std::string suffix = FileIO::getSuffix( filename );

        if ( FileIO::canCreate( suffix ) )
        {
            // okay, we can use FileIO class from factory

            common::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );
    
            fileIO->readArray( mLocalValues, filename );
        }
        else
        {
            // ToDo: readFromFile( filename + ".<suffix>" ) for all known suffixes

            COMMON_THROWEXCEPTION( "File : " << filename << ", unknown file type " << suffix )
        }
    }
    else
    {
        // other processors have to clear their local values
        mLocalValues.clear();
    }

    IndexType numElements = mLocalValues.size();
    comm->bcast( &numElements, 1, host );
    DistributionPtr dist( new CyclicDistribution( numElements, numElements, comm ) );
    SCAI_ASSERT_EQ_DEBUG( dist->getLocalSize(), mLocalValues.size(), "wrong distribution" );
    SCAI_ASSERT_EQ_DEBUG( dist->getGlobalSize(), numElements, "wrong distribution" );
    // this is safe, we have allocated it correctly
    setDistributionPtr( dist );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const _HArray& localValues, DistributionPtr distribution )
    : Vector( distribution )
{
    SCAI_ASSERT_EQ_ERROR( localValues.size(), distribution->getLocalSize(), "size mismatch" )
    HArrayUtils::assign( mLocalValues, localValues ); // can deal with type conversions
}

/* ------------------------------------------------------------------------- */

/*
 * Constructors with Expressions as arguments
 */

// linear algebra expression: a*x
template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Scalar, Vector, Times>& expression )

    : Vector( expression.getArg2() )
{
    SCAI_LOG_INFO( logger, "Constructor( alpha * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x+b*y, inherit distribution/context from vector x

template<typename ValueType>
DenseVector<ValueType>::DenseVector(
    const Expression<Expression<Scalar, Vector, Times>, Expression<Scalar, Vector, Times>, Plus>& expression ) //Expression_SV_SV

    : Vector( expression.getArg1().getArg2() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x + beta * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector(
    const Expression<Expression<Scalar, Expression<Matrix, Vector, Times>, Times>, Expression<Scalar, Vector, Times>, Plus>& expression ) //Expression_SMV_SV

    : Vector( expression.getArg1().getArg2().getArg1().getRowDistributionPtr(),
              expression.getArg1().getArg2().getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector(
    const Expression<Expression<Scalar, Expression<Vector, Matrix, Times>, Times>, Expression<Scalar, Vector, Times>, Plus>& expression ) //Expression_SVM_SV
    : Vector( expression.getArg1().getArg2().getArg2().getColDistributionPtr(),
              expression.getArg1().getArg2().getArg2().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Scalar, Expression<Matrix, Vector, Times>, Times>& expression )

    : Vector( expression.getArg2().getArg1().getRowDistributionPtr(),
              expression.getArg2().getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x*A, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Scalar, Expression<Vector, Matrix, Times>, Times>& expression )
    : Vector( expression.getArg2().getArg2().getColDistributionPtr(),
              expression.getArg2().getArg2().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A )" )
    Vector::operator=( expression );
}

// linear algebra expression: A*x, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Matrix, Vector, Times>& expression )
    : Vector( expression.getArg1().getRowDistributionPtr(), expression.getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( A * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: x*A, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Vector, Matrix, Times>& expression )
    : Vector( expression.getArg2().getColDistributionPtr(), expression.getArg2().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( x * A )" )
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
    assign( other );
    return *this;
}

template<typename ValueType>
DenseVector<ValueType>& DenseVector<ValueType>::operator=( const Scalar value )
{
    assign( value );
    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::scalar::ScalarType DenseVector<ValueType>::getValueType() const
{
    return TypeTraits<ValueType>::stype;
}

template<typename ValueType>
void DenseVector<ValueType>::buildValues( _HArray& values ) const
{
    // size of values will be local size of vecotr
    HArrayUtils::assign( values, mLocalValues );
}

template<typename ValueType>
void DenseVector<ValueType>::setValues( const _HArray& values )
{
    SCAI_ASSERT_ERROR(
        values.size() == mLocalValues.size(),
        "Size of values = " << values.size() << ", does not match local size of vector = " << mLocalValues.size() );
    HArrayUtils::assign( mLocalValues, values );
}

template<typename ValueType>
DenseVector<ValueType>* DenseVector<ValueType>::copy() const
{
    // create a new dense vector with the copy constructor
    return new DenseVector<ValueType>( *this );
}

template<typename ValueType>
DenseVector<ValueType>* DenseVector<ValueType>::newVector() const
{
    common::unique_ptr<DenseVector<ValueType> > vector( new DenseVector<ValueType>() );
    vector->setContextPtr( this->getContextPtr() );
    return vector.release();
}

template<typename ValueType>
void DenseVector<ValueType>::updateHalo( const dmemo::Halo& halo ) const
{
    const IndexType haloSize = halo.getHaloSize();
    SCAI_LOG_DEBUG( logger, "Acquiring halo write access on " << *mContext )
    mHaloValues.clear();
    WriteAccess<ValueType> haloAccess( mHaloValues, mContext );
    haloAccess.reserve( haloSize );
    haloAccess.release();
    getDistribution().getCommunicator().updateHalo( mHaloValues, mLocalValues, halo );
}

template<typename ValueType>
tasking::SyncToken* DenseVector<ValueType>::updateHaloAsync( const dmemo::Halo& halo ) const
{
    const IndexType haloSize = halo.getHaloSize();
    // create correct size of Halo
    {
        WriteOnlyAccess<ValueType> haloAccess( mHaloValues, mContext, haloSize );
    }
    return getDistribution().getCommunicator().updateHaloAsync( mHaloValues, mLocalValues, halo );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::getValue( IndexType globalIndex ) const
{
    SCAI_LOG_TRACE( logger, *this << ": getValue( globalIndex = " << globalIndex << " )" )
    ValueType myValue = static_cast<ValueType>( 0.0 );
    const IndexType localIndex = getDistribution().global2local( globalIndex );

    if ( localIndex != nIndex )
    {
        myValue = mLocalValues[localIndex];
    }

    ValueType allValue = getDistribution().getCommunicator().sum( myValue );
    SCAI_LOG_TRACE( logger, "myValue = " << myValue << ", allValue = " << allValue )
    return Scalar( allValue );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::setValue( const IndexType globalIndex, const Scalar value )
{
    SCAI_LOG_TRACE( logger, *this << ": setValue( globalIndex = " << globalIndex << " ) = " <<  value )

    const IndexType localIndex = getDistribution().global2local( globalIndex );

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
void DenseVector<ValueType>::conj()
{
    HArrayUtils::conj( mLocalValues, mContext );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::l1Norm() const
{
    ValueType localL1Norm = mLocalValues.l1Norm();
    return Scalar( getDistribution().getCommunicator().sum( localL1Norm ) );
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
    mHaloValues.swap( otherPtr->mHaloValues );
}

template<typename ValueType>
void DenseVector<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "DenseVector<" << getValueType() << ">" << "( size = " << size() << ", local = " << mLocalValues.size()
           << ", dist = " << getDistribution() << ", loc  = " << *getContextPtr() << " )";
}

template<typename ValueType>
void DenseVector<ValueType>::assign( const Expression_SV_SV& expression )
{
    const Expression_SV& exp1 = expression.getArg1();
    const Expression_SV& exp2 = expression.getArg2();
    const ValueType alpha = exp1.getArg1().getValue<ValueType>();
    const Vector& x = exp1.getArg2();
    const ValueType beta = exp2.getArg1().getValue<ValueType>();
    const Vector& y = exp2.getArg2();
    SCAI_LOG_INFO( logger, "z = " << alpha << " * x + " << beta << " * y, with  x = " << x << ", y = " << y << ", z = " << *this )
    SCAI_LOG_DEBUG( logger, "dist of x = " << x.getDistribution() )
    SCAI_LOG_DEBUG( logger, "dist of y = " << y.getDistribution() )

    if ( x.getDistribution() != y.getDistribution() )
    {
        COMMON_THROWEXCEPTION(
            "distribution do not match for z = alpha * x + beta * y, z = " << *this << " , x = " << x << " , y = " << y )
    }

    if ( x.getDistribution() != getDistribution() || x.size() != size() )
    {
        allocate( x.getDistributionPtr() );
    }

    if ( typeid( *this ) == typeid( x ) && typeid( *this ) == typeid( y ) )
    {
        const DenseVector<ValueType>& denseX = dynamic_cast<const DenseVector<ValueType>&>( x );
        const DenseVector<ValueType>& denseY = dynamic_cast<const DenseVector<ValueType>&>( y );

        if ( mLocalValues.size() != denseX.mLocalValues.size() )
        {
            SCAI_LOG_DEBUG( logger, "resize local values of z = this" )
            mLocalValues.clear();
            WriteAccess<ValueType> localAccess( mLocalValues, mContext );
            localAccess.resize( denseX.mLocalValues.size() );
        }

#ifdef NOT_SWITCHED_ON
        {
            // useful output to identify aliases between arguments, write should be the last one
            ReadAccess<ValueType> rX( denseX.mLocalValues, mContext );
            ReadAccess<ValueType> rY( denseY.mLocalValues, mContext );
            WriteAccess<ValueType> rZ( mLocalValues, mContext );
            SCAI_LOG_DEBUG( logger, " z = " << rZ.get() << ", x = " << rX.get() << ", y = " << rY.get() )
        }
#endif
        SCAI_LOG_DEBUG( logger, "call arrayPlusArray" )
        utilskernel::HArrayUtils::arrayPlusArray( mLocalValues, alpha, denseX.mLocalValues, beta, denseY.mLocalValues, mContext );
    }
    else
    {
        COMMON_THROWEXCEPTION(
            "Can not calculate z = alpha * x + beta * y, z = " << *this << ", x = " << x << ", y = " << y << " because of type mismatch." );
    }
}

template<typename ValueType>
Scalar DenseVector<ValueType>::dotProduct( const Vector& other ) const
{
    SCAI_REGION( "Vector.Dense.dotP" )
    SCAI_LOG_INFO( logger, "Calculating dot product for " << *this << " * " << other )

    // add other->getVectorKind() == DENSE, if sparse is also supported

    if ( this->getValueType() == other.getValueType() )
    {
        if ( getDistribution() != other.getDistribution() )
        {
            COMMON_THROWEXCEPTION( "distribution do not match for this * other, this = " << *this << " , other = " << other )
        }

        const DenseVector<ValueType>* denseOther = dynamic_cast<const DenseVector<ValueType>*>( &other );
        SCAI_ASSERT_DEBUG( denseOther, "dynamic_cast failed for other = " << other )
        SCAI_LOG_DEBUG( logger, "Calculating local dot product at " << *mContext )
        const IndexType localSize = mLocalValues.size();
        SCAI_ASSERT_EQ_DEBUG( localSize, getDistribution().getLocalSize(), "size mismatch" )
        const ValueType localDotProduct = mLocalValues.dotProduct( denseOther->mLocalValues );
        SCAI_LOG_DEBUG( logger, "Calculating global dot product form local dot product = " << localDotProduct )
        ValueType dotProduct = getDistribution().getCommunicator().sum( localDotProduct );
        SCAI_LOG_DEBUG( logger, "Global dot product = " << dotProduct )
        return Scalar( dotProduct );
    }

    COMMON_THROWEXCEPTION(
        "Can not calculate a dot product of " << typeid( *this ).name() << " and " << typeid( other ).name() )
}

template<typename ValueType>
void DenseVector<ValueType>::allocate( DistributionPtr distribution )
{
    setDistributionPtr( distribution );
    // resize the local values at its context
    WriteOnlyAccess<ValueType> dummyWAccess( mLocalValues, mContext, getDistribution().getLocalSize() );
    // local values are likely to be uninitialized
}

template<typename ValueType>
void DenseVector<ValueType>::assign( const Vector& other )
{
    setDistributionPtr( other.getDistributionPtr() );
    // Note: we cannot use other.getLocalValues() as other might be a sparse vector
    other.buildLocalValues( mLocalValues ); // but this works fine also for format conversion
}

template<typename ValueType>
void DenseVector<ValueType>::assign( const Scalar value )
{
    SCAI_LOG_DEBUG( logger, *this << ": assign " << value )
    // assign the scalar value on the home of this dense vector.
    HArrayUtils::setScalar( mLocalValues, value.getValue<ValueType>(), utilskernel::reduction::COPY, mContext );
}

template<typename ValueType>
void DenseVector<ValueType>::assign( const _HArray& localValues, DistributionPtr dist )
{
    SCAI_LOG_INFO( logger, "assign vector with localValues = " << localValues << ", dist = " << *dist )
    SCAI_ASSERT_EQ_ERROR( localValues.size(), dist->getLocalSize(), "size mismatch" )
    setDistributionPtr( dist );
    HArrayUtils::assign( mLocalValues, localValues );
}

template<typename ValueType>
void DenseVector<ValueType>::buildLocalValues( _HArray& localValues ) const
{
    HArrayUtils::assign( localValues, mLocalValues );
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
    utilskernel::HArrayUtils::invert( mLocalValues, this->getContextPtr() );
}

template<typename ValueType>
size_t DenseVector<ValueType>::getMemoryUsage() const
{
    // Note: memory of mHaloValues is not counted, is just a temporary
    size_t memoryUsage = sizeof( ValueType ) * mLocalValues.size();
    return getDistribution().getCommunicator().sum( memoryUsage );
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseVector<ValueType>::redistribute( DistributionPtr distribution )
{
    SCAI_ASSERT_EQ_ERROR( size(), distribution->getGlobalSize(), "global size mismatch between old/new distribution" )

    if ( getDistribution() == *distribution )
    {
        SCAI_LOG_INFO( logger, *this << " redistribute to same distribution " << *distribution )
        // we can keep local/global values, but just set dist pointer
        setDistributionPtr( distribution );
    }
    else if ( getDistribution().isReplicated() )
    {
        SCAI_LOG_INFO( logger, *this << ": replicated vector" << " will be localized to " << *distribution )
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
        SCAI_LOG_INFO( logger, *this << " will be replicated" )
        // replicate a distributed vector
        HArray<ValueType> globalValues;
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
        SCAI_LOG_INFO( logger, *this << " will be redistributed to " << *distribution )
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
void DenseVector<ValueType>::writeToFile(
    const std::string& fileName,
    const std::string& fileType,               /* = "", take IO type by suffix   */
    const common::scalar::ScalarType dataType, /* = UNKNOWN, take defaults of IO type */
    const FileIO::FileMode fileMode            /* = DEFAULT_MODE */ 
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
DenseVector<ValueType>::DenseVector( const DenseVector<ValueType>& other )
    : Vector( other )
{
    // implementation here can be simpler as DenseVector( const Vector& other )
    SCAI_LOG_INFO( logger,
                   "Copy of vector of global size " << size() << ", local size " << getDistribution().getLocalSize() )
    mLocalValues = other.getLocalValues();
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DenseVector, SCAI_ARITHMETIC_HOST )

} /* end namespace lama */

} /* end namespace scai */
