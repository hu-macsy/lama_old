/**
 * @file DenseVector.hpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief DenseVector.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * $Id$
 */

// hpp
#include <lama/DenseVector.hpp>

// others
#include <lama/LAMAArrayUtils.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/ContextAccess.hpp>
//#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/CommunicatorFactory.hpp>

#include <lama/distribution/NoDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>
#include <lama/distribution/Redistributor.hpp>

#include <lama/matrix/Matrix.hpp>
#include <lama/expression/Expression.hpp>

#include <boost/scoped_array.hpp>

namespace lama
{

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, DenseVector<T>::logger, "Vector.DenseVector" )

template<typename T>
DenseVector<T>::DenseVector( ) : Vector( 0 ), mLocalValues()
{
}

template<typename T>
DenseVector<T>::DenseVector( DistributionPtr distribution )
    : Vector( distribution ), mLocalValues( distribution->getLocalSize() )
{
    LAMA_LOG_INFO( logger, "Construct dense vector, size = " << distribution->getGlobalSize()
                   // << ", type = " << typename(T)
                   << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", no initialization" )
}

template<typename T>
DenseVector<T>::DenseVector( const IndexType size, const ValueType value )
    : Vector( size ), mLocalValues( size, value )
{
    LAMA_LOG_INFO( logger, "Construct dense vector, size = " << size << ", init =" << value )
}

template<typename T>
DenseVector<T>::DenseVector( DistributionPtr distribution, const ValueType value )
    : Vector( distribution ), mLocalValues( distribution->getLocalSize(), value )
{
    LAMA_LOG_INFO( logger,
                   "Construct dense vector, size = " << distribution->getGlobalSize() << ", distribution = " << *distribution << ", local size = " << distribution->getLocalSize() << ", init = " << value )
}

template<typename T>
DenseVector<T>::DenseVector( const Vector& other )
    : Vector( other )
{
    allocate( getDistributionPtr() );
    assign( other );
}

template<typename T>
DenseVector<T>::DenseVector( const Vector& other, DistributionPtr distribution )

    : Vector( other )
{
    assign( other );
    redistribute( distribution );
}

template<typename T>
DenseVector<T>::DenseVector( const std::string& filename )
{
    LAMA_LOG_INFO( logger, "Construct dense vector from file " << filename )

    // Take the current default communicator

    CommunicatorPtr comm = CommunicatorFactory::get();

    IndexType myRank = comm->getRank();
    IndexType host = 0; // reading processor

    IndexType numElements = 0; // will be the size of the vector

    File::FileType fileType = File::XDR;
    long dataTypeSize = -1;
    std::string vecFileName( filename + ".vec" );

    if ( myRank == host )
    {
        readVectorHeader( filename + ".frv", fileType, dataTypeSize );

        // broadcast the size to all other processors

        numElements = size();
    }

    comm->bcast( &numElements, 1, host );

    // host gets all elements

    DistributionPtr dist( new CyclicDistribution( numElements, numElements, comm ) );

    allocate( dist );

    if ( myRank == host )
    {
        // host owns all elements, all others have no elements

        switch ( fileType )
        {
        case File::FORMATTED: //ASCII
            readVectorFromFormattedFile( vecFileName );
            break;
        case File::BINARY: //Binary without number of elements in the vector file
            readVectorFromBinaryFile( vecFileName, dataTypeSize );
            break;
        case File::XDR: //XDR following the IEEE standard
            readVectorFromXDRFile( vecFileName, dataTypeSize );
            break;
        default:
            throw Exception( "Unknown File Type." );
        }

        // @TODO do not throw exception as application will hang
    }
}

/* ------------------------------------------------------------------------- */

template<typename T>
DenseVector<T>::DenseVector( const _LAMAArray& localValues, DistributionPtr distribution )
    : Vector( distribution )
{
    LAMA_ASSERT_EQUAL_ERROR( localValues.size(), distribution->getLocalSize() )

    LAMAArrayUtils::assign( mLocalValues, localValues ); // can deal with type conversions
}

/* ------------------------------------------------------------------------- */

/*
 * Constructors with Expressions as arguments
 */

// linear algebra expression: a*x
template<typename T>
DenseVector<T>::DenseVector( const Expression<Scalar,Vector,Times>& expression )

    : Vector( expression.getArg2() )
{
    LAMA_LOG_INFO( logger, "Constructor( alpha * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x+b*y, inherit distribution/context from vector x

template<typename T>
DenseVector<T>::DenseVector(
    const Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>& expression )

    : Vector( expression.getArg1().getArg2() )
{
    allocate( getDistributionPtr() );
    LAMA_LOG_INFO( logger, "Constructor( alpha * x + beta * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename T>
DenseVector<T>::DenseVector(
    const Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>& expression )

    : Vector( expression.getArg1().getArg2().getArg1().getDistributionPtr(),
              expression.getArg1().getArg2().getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    LAMA_LOG_INFO( logger, "Constructor( alhpa * A * x + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x, inherit distribution/context from matrix A

template<typename T>
DenseVector<T>::DenseVector( const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& expression )

    : Vector( expression.getArg2().getArg1().getDistributionPtr(),
              expression.getArg2().getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    LAMA_LOG_INFO( logger, "Constructor( alpha * A * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: A*x, inherit distribution/context from matrix A

template<typename T>
DenseVector<T>::DenseVector( const Expression<Matrix,Vector,Times>& expression )

    : Vector( expression.getArg1().getDistributionPtr(), expression.getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    LAMA_LOG_INFO( logger, "Constructor( A * x )" )
    Vector::operator=( expression );
}

/* ------------------------------------------------------------------------- */

template<typename T>
DenseVector<T>::~DenseVector()
{
}

/* ------------------------------------------------------------------------- */

template<typename T>
DenseVector<T>& DenseVector<T>::operator=( const DenseVector<T>& other )
{
    assign( other );
    return *this;
}

template<typename T>
DenseVector<T>& DenseVector<T>::operator=( const Scalar value )
{
    assign( value );
    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename T>
Scalar::ScalarType DenseVector<T>::getValueType() const
{
    return Scalar::getType<ValueType>();
}

template<typename T>
void DenseVector<T>::buildValues( _LAMAArray& values ) const
{
    // size of values will be local size of vecotr

    LAMAArrayUtils::assign( values, mLocalValues );
}

template<typename T>
void DenseVector<T>::setValues( const _LAMAArray& values )
{
    LAMA_ASSERT_ERROR(
        values.size() == mLocalValues.size(),
        "Size of values = " << values.size() << ", does not match local size of vector = " << mLocalValues.size() );

    LAMAArrayUtils::assign( mLocalValues, values );
}

template<typename T>
Vector* DenseVector<T>::create() const
{
    return create( getDistributionPtr() );
}

template<typename T>
Vector* DenseVector<T>::create( DistributionPtr distribution ) const
{
    LAMA_LOG_INFO( logger, "DenseVector<T>::create" )

    std::auto_ptr<Vector> newDenseVector( new DenseVector<T>( distribution ) );

    newDenseVector->setContext( mContext );

    // give back the new vector and its ownership

    return newDenseVector.release();
}

template<typename T>
void DenseVector<T>::updateHalo( const Halo& halo ) const
{
    const IndexType haloSize = halo.getHaloSize();

    LAMA_LOG_DEBUG( logger, "Acquiring halo write access on " << *mContext )

    mHaloValues.clear();
    WriteAccess<T> haloAccess( mHaloValues, mContext );

    haloAccess.reserve( haloSize );
    haloAccess.release();

    getDistribution().getCommunicator().updateHalo( mHaloValues, mLocalValues, halo );
}

template<typename T>
SyncToken* DenseVector<T>::updateHaloAsync( const Halo& halo ) const
{
    const IndexType haloSize = halo.getHaloSize();

    // create correct size of Halo

    {
        WriteOnlyAccess<T> haloAccess( mHaloValues, mContext, haloSize );
    }

    return getDistribution().getCommunicator().updateHaloAsync( mHaloValues, mLocalValues, halo );
}

template<typename T>
Scalar DenseVector<T>::getValue( IndexType globalIndex ) const
{
    LAMA_LOG_TRACE( logger, *this << ": getValue( globalIndex = " << globalIndex << " )" )
    ValueType myValue = 0.0;
    const IndexType localIndex = getDistribution().global2local( globalIndex );

    if ( localIndex != nIndex )
    {
        HostReadAccess<ValueType> localAccess( mLocalValues );

        LAMA_LOG_TRACE( logger, "index "<< globalIndex << " is local " << localIndex )
        myValue = localAccess[localIndex];
    }
    ValueType allValue = getDistribution().getCommunicator().sum( myValue );
    LAMA_LOG_TRACE( logger, "myValue = " << myValue << ", allValue = " << allValue )
    return Scalar( allValue );
}

template<typename T>
Scalar DenseVector<T>::min() const
{
    //TODO: need a interface function for this
    HostReadAccess<ValueType> localValues( mLocalValues );
    ValueType localMin = localValues[0];
    #pragma omp parallel
    {
        ValueType myLocalMin = localMin;
        #pragma omp for
        for ( IndexType i = 0; i < localValues.size(); ++i )
        {
            myLocalMin = std::min( localValues[i], myLocalMin );
        }

        #pragma omp critical
        {
            localMin = std::min( localMin, myLocalMin );
        }
    }
    return getDistribution().getCommunicator().min( localMin );
}

template<typename T>
Scalar DenseVector<T>::max() const
{
    LAMA_ASSERT_ERROR( mLocalValues.size() > 0 , "no local values for max" )

    //TODO: need a interface function for this
    HostReadAccess<ValueType> localValues( mLocalValues );
    ValueType localMax = localValues[0];
    #pragma omp parallel
    {
        ValueType myLocalMax = localMax;
        #pragma omp for
        for ( IndexType i = 0; i < localValues.size(); ++i )
        {
            myLocalMax = std::max( localValues[i], myLocalMax );
        }

        #pragma omp critical
        {
            localMax = std::max( localMax, myLocalMax );
        }
    }
    return getDistribution().getCommunicator().max( localMax );
}

/* ------------------------------------------------------------------------- */

template<typename T>
Scalar DenseVector<T>::l1Norm() const
{
    IndexType nnu = mLocalValues.size();

    ValueType localL1Norm = static_cast<ValueType>( 0 );

    if ( nnu > 0 )
    {
        //choose preferred context

        ContextPtr loc = mContext;

        // get function pointer BLAS::BLAS1<T>::asum in appropriate context

        LAMA_INTERFACE_FN_DEFAULT_T( asum, loc, BLAS, BLAS1, T )

        ReadAccess<T> read( mLocalValues, loc );

        LAMA_CONTEXT_ACCESS( loc )

        localL1Norm = asum( nnu, read.get(), 1, NULL );
    }

    return getDistribution().getCommunicator().sum( localL1Norm );
}

/* ------------------------------------------------------------------------- */

template<typename T>
Scalar DenseVector<T>::l2Norm() const
{
    IndexType nnu = mLocalValues.size();
    
    ValueType localDotProduct = static_cast<ValueType>( 0 );

    if ( nnu > 0 )
    {
        // choose preferred context as context of vector, might be changed by availability

        ContextPtr loc = mContext;

        // get function pointer BLAS::BLAS1<T>::dot in appropriate context

        LAMA_INTERFACE_FN_DEFAULT_T( dot, loc, BLAS, BLAS1, T )

        ReadAccess<T> read( mLocalValues, loc );

        LAMA_CONTEXT_ACCESS( mContext )

        localDotProduct = dot( nnu, read.get(), 1, read.get(), 1, NULL );
    }

    ValueType globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );

    return sqrt( globalDotProduct );
}

/* ------------------------------------------------------------------------- */

template<typename T>
Scalar DenseVector<T>::maxNorm() const
{
    IndexType nnu = mLocalValues.size(); // number of local rows

    ValueType localMaxNorm = static_cast<ValueType>( 0 );

    if ( nnu > 0 )
    {
        ContextPtr loc = mContext; // loc might be set to Host

        LAMA_INTERFACE_FN_DEFAULT_T( absMaxVal, loc, Utils, Reductions, ValueType )

        ReadAccess<T> read( mLocalValues, loc );

        LAMA_CONTEXT_ACCESS( loc )

        localMaxNorm = absMaxVal( read.get(), nnu );
    }

    const Communicator& comm = getDistribution().getCommunicator();

    ValueType globalMaxNorm = comm.max( localMaxNorm );

    LAMA_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm of " << nnu << " elements: "
                   << localMaxNorm << ", max norm global = " << globalMaxNorm )

    return globalMaxNorm;
}

/* ------------------------------------------------------------------------- */

template<typename T>
void DenseVector<T>::swap( Vector& other )
{
    LAMA_LOG_DEBUG( logger, "swap:" << *this << " with " << other )

    DenseVector* otherPtr = dynamic_cast<DenseVector*>( &other );

    if ( !otherPtr )
    {
        LAMA_THROWEXCEPTION( "Tried to swap with a Vector of a different type." )
    }
    Vector::swapVector( other );
    mLocalValues.swap( otherPtr->mLocalValues );
    mHaloValues.swap( otherPtr->mHaloValues );
}

template<typename T>
void DenseVector<T>::writeAt( std::ostream& stream ) const
{
    stream << "DenseVector<" << getValueType() << ">(" << size() << ", local = " << mLocalValues.size() << ", dist = "
           << getDistribution() << ")";
}

template<typename T>
void DenseVector<T>::vectorPlusVector(
    ContextPtr context,
    LAMAArrayView<T> result,
    const T alpha,
    const LAMAArrayConstView<T> x,
    const T beta,
    const LAMAArrayConstView<T> y )
{
    LAMA_LOG_DEBUG( logger,
                    "vectorPlusVector: result:" << result << " = " << alpha << " * x:" << x << " + " << beta << " * y:" << y )

    // get function pointers, do not use fallbacks here

    LAMA_INTERFACE_FN_T( scal, context, BLAS, BLAS1, T )
    LAMA_INTERFACE_FN_T( axpy, context, BLAS, BLAS1, T )
    LAMA_INTERFACE_FN_T( sum, context, BLAS, BLAS1, T )

    const IndexType nnu = result.size();

    if ( result == x && result == y ) //result = alpha * result + beta * result
    {
        //result = alpha * result + beta * result
        //=>
        //result = ( alpha + beta ) * result
        //=>
        //result *= ( alpha + beta )

        LAMA_LOG_DEBUG( logger,
                        "vectorPlusVector: x = y = result, result *= " << "alpha(" << alpha << ") + beta(" << beta << ")" )

        WriteAccess<T> resultAccess( result, context, true );

        LAMA_CONTEXT_ACCESS( context )
        scal( nnu, alpha + beta, resultAccess.get(), 1, NULL );
    }
    else if ( result == x ) //result = alpha * result + beta * y
    {
        ReadAccess<T> yAccess( y, context );
        WriteAccess<T> resultAccess( result, context, true );

        if ( beta == 0.0 )
        {
            LAMA_LOG_DEBUG( logger, "vectorPlusVector: result *= alpha" )

            if ( alpha != 1.0 ) // result *= alpha
            {
                LAMA_CONTEXT_ACCESS( context )
                scal( nnu, alpha, resultAccess.get(), 1, NULL );
            }
            else
            {
                // do nothing: result = 1 * result
            }
        }
        else if ( beta == 1.0 ) // result = alpha * result + y
        {
            LAMA_LOG_DEBUG( logger, "vectorPlusVector: result = alpha * result + y" )

            if ( alpha != 1.0 ) // result = alpha * result + y
            {
                // result *= alpha
                LAMA_CONTEXT_ACCESS( context )
                scal( nnu, alpha, resultAccess.get(), 1, NULL );
            }
            // result += y
            LAMA_CONTEXT_ACCESS( context )
            axpy( nnu, 1/*alpha*/, yAccess.get(), 1, resultAccess.get(), 1, NULL );
        }
        else // beta != 1.0 && beta != 0.0 --> result = alpha * result + beta * y
        {
            LAMA_LOG_DEBUG( logger,
                            "vectorPlusVector: result = alpha(" << alpha << ")" << " * result + beta(" << beta << ") * y" )

            if ( alpha != 1.0 )
            {
                LAMA_CONTEXT_ACCESS( context )
                scal( nnu, alpha, resultAccess.get(), 1, NULL );
            }
            LAMA_CONTEXT_ACCESS( context )
            axpy( nnu, beta, yAccess.get(), 1, resultAccess.get(), 1, NULL );
        }
    }
    else if ( result == y ) // result = alpha * x + beta * result
    {
        LAMA_LOG_DEBUG( logger,
                        "vectorPlusVector: result = alpha(" << alpha << ")" << " * x + beta(" << beta << ") * result" )

        // so we do here:  result = beta * result, result += alpha * x

        ReadAccess<T> xAccess( x, context );
        WriteAccess<T> resultAccess( result, context, true );

        if ( beta != 1.0 ) // result = [alpha * x + ] beta * result
        {
            // result *= beta
            LAMA_CONTEXT_ACCESS( context )
            scal( nnu, beta, resultAccess.get(), 1, NULL );
        }

        if ( alpha != 0.0 )
        {
            // result = alpha * x + result
            LAMA_CONTEXT_ACCESS( context )
            axpy( nnu, alpha, xAccess.get(), 1, resultAccess.get(), 1, NULL );
        }
    }
    else // result = alpha * x + beta * y
    {
        LAMA_LOG_DEBUG( logger,
                        "vectorPlusVector: result = alpha(" << alpha << ")" << " * x + beta(" << beta << ") * y" )

        ReadAccess<T> xAccess( x, context );
        ReadAccess<T> yAccess( y, context );
        // no need to keep old values of result
        WriteAccess<T> resultAccess( result, context, false );

        LAMA_CONTEXT_ACCESS( context )
        sum( nnu, alpha, xAccess.get(), beta, yAccess.get(), resultAccess.get(), NULL );
    }

    LAMA_LOG_INFO( logger, "vectorPlusVector done" )
}

template<typename T>
void DenseVector<T>::assign(
    const Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>& expression )
{
    const ValueType alpha = expression.getArg1().getArg1().getValue<ValueType>();
    const Vector& x = expression.getArg1().getArg2();
    const ValueType beta = expression.getArg2().getArg1().getValue<ValueType>();
    const Vector& y = expression.getArg2().getArg2();

    LAMA_LOG_DEBUG( logger, *this << ": assign" << alpha << " * x:" << x << " + " << beta << " * y:" << y )

    LAMA_LOG_DEBUG( logger, "dist of x = " << x.getDistribution() )
    LAMA_LOG_DEBUG( logger, "dist of y = " << y.getDistribution() )

    if ( x.getDistribution() != y.getDistribution() )
    {
        LAMA_THROWEXCEPTION(
            "distribution do not match for z = alpha * x + beta * y, z = "<< *this <<" , x = "<< x <<" , y = "<< y )
    }

    if ( x.getDistribution() != getDistribution() || x.size() != size() )
    {
        resize( x.getDistributionPtr() );
    }

    if ( typeid( *this ) == typeid( x ) && typeid( *this ) == typeid( y ) )
    {
        const DenseVector<ValueType>& denseX = dynamic_cast<const DenseVector<ValueType>&>( x );
        const DenseVector<ValueType>& denseY = dynamic_cast<const DenseVector<ValueType>&>( y );

        if ( mLocalValues.size() != denseX.mLocalValues.size() )
        {
            mLocalValues.clear();
            WriteAccess<ValueType> localAccess( mLocalValues, mContext );
            localAccess.resize( denseX.mLocalValues.size() );
        }

#ifdef LAMA_LOG_DEBUG_ENABLED
        {
            // useful output to identify aliases between arguments

            ReadAccess<ValueType> rZ ( mLocalValues, mContext );
            ReadAccess<ValueType> rX ( denseX.mLocalValues, mContext );
            ReadAccess<ValueType> rY ( denseY.mLocalValues, mContext );

            LAMA_LOG_DEBUG( logger, " z = " << rZ.get() << ", x = " << rX.get() << ", y = " << rY.get() )
        }
#endif

        vectorPlusVector( mContext, mLocalValues, alpha, denseX.mLocalValues, beta, denseY.mLocalValues );
    }
    else
    {
        LAMA_THROWEXCEPTION(
            "Can not calculate z = alpha * x + beta * y, z = " << *this << ", x = " << x << ", y = " << y << " because of type mismatch." );
    }
}

template<typename T>
Scalar DenseVector<T>::dotProduct( const Vector& other ) const
{
    LAMA_LOG_INFO( logger, "Calculating dot product for " << *this << " * " << other )
    if ( typeid( *this ) == typeid( other ) )
    {
        if ( getDistribution() != other.getDistribution() )
        {
            LAMA_THROWEXCEPTION(
                "distribution do not match for this * other, this = "<< *this <<" , other = "<< other )
        }
        const DenseVector<ValueType>& denseOther = dynamic_cast<const DenseVector<ValueType>&>( other );

        LAMA_LOG_DEBUG( logger, "Calculating local dot product at " << *mContext )

        ContextPtr loc = mContext;   // prefered location is context of this vector

        LAMA_INTERFACE_FN_DEFAULT_T( dot, loc, BLAS, BLAS1, T );

        // Now do the dot production at location loc ( might have been changed to other location  )

        ReadAccess<T> localRead( mLocalValues, loc );
        ReadAccess<T> otherRead( denseOther.mLocalValues, loc );

        LAMA_CONTEXT_ACCESS( loc )

        const IndexType localSize = mLocalValues.size();

        LAMA_ASSERT_EQUAL_DEBUG( localSize, getDistribution().getLocalSize() )

        const ValueType localDotProduct = dot( localSize, localRead.get(), 1, otherRead.get(), 1, NULL );

        LAMA_LOG_DEBUG( logger, "Calculating global dot product form local dot product = " << localDotProduct )

        ValueType dotProduct = getDistribution().getCommunicator().sum( localDotProduct );

        LAMA_LOG_DEBUG( logger, "Global dot product = " << dotProduct )

        return dotProduct;
    }

    LAMA_THROWEXCEPTION(
        "Can not calculate a dot product of "<< typeid( *this ).name() <<" and "<< typeid( other ).name() )
}

template<typename T>
void DenseVector<T>::allocate( DistributionPtr distribution )
{
    setDistributionPtr( distribution );

    // resize the local values at its context

    WriteOnlyAccess<T> dummyWAccess( mLocalValues, mContext, getDistribution().getLocalSize() );

    // local values are likely to be uninitialized
}

template<typename T>
void DenseVector<T>::assign( const Vector& other )
{
    setDistributionPtr( other.getDistributionPtr() );

    // Note: we cannot use other.getLocalValues() as other might be a sparse vector

    other.buildLocalValues( mLocalValues ); // but this works fine also for format conversion
}

template<typename T>
void DenseVector<T>::assign( const Scalar value )
{
    LAMA_LOG_DEBUG( logger, *this << ": assign " << value )

    ContextPtr ctx = mLocalValues.getValidContext( mContext->getType() );
    LAMAArrayUtils::assign( mLocalValues, value, ctx );
}

template<typename T>
void DenseVector<T>::assign( const _LAMAArray& localValues, DistributionPtr dist )
{
    LAMA_LOG_INFO( logger, "assign vector with localValues = " << localValues << ", dist = " << *dist )

    LAMA_ASSERT_EQUAL_ERROR( localValues.size(), dist->getLocalSize() )

    setDistributionPtr( dist );
    LAMAArrayUtils::assign( mLocalValues, localValues );
}

template<typename T>
void DenseVector<T>::buildLocalValues( _LAMAArray& localValues ) const
{
    LAMAArrayUtils::assign( localValues, mLocalValues );
}

template<typename T>
void DenseVector<T>::prefetch( const ContextPtr location ) const
{
    mLocalValues.prefetch( location );
}

template<typename T>
void DenseVector<T>::wait() const
{
    mLocalValues.wait();
}

template<typename T>
void DenseVector<T>::invert()
{
    const IndexType size = mLocalValues.size();

    const ContextPtr loc = getContext();

    LAMA_INTERFACE_FN_T( invert, loc, Utils, Math, T )

    WriteAccess<ValueType> wValues( mLocalValues, loc );

    LAMA_CONTEXT_ACCESS( loc );

    invert( wValues.get(), size );
}

template<typename T>
size_t DenseVector<T>::getMemoryUsage() const
{
    // Note: memory of mHaloValues is not counted, is just a temporary

    size_t memoryUsage = sizeof(T) * mLocalValues.size();
    return getDistribution().getCommunicator().sum( memoryUsage );
}

/* ------------------------------------------------------------------------ */

template<typename T>
void DenseVector<T>::redistribute( DistributionPtr distribution )
{
    LAMA_ASSERT_EQUAL_ERROR( size(), distribution->getGlobalSize() )

    if ( getDistribution() == *distribution )
    {
        LAMA_LOG_INFO( logger, *this << " redistribute to same distribution " << *distribution )

        // we can keep local/global values, but just set dist pointer

        setDistributionPtr( distribution );
    }

    else if ( getDistribution().isReplicated() )

    {
        LAMA_LOG_INFO( logger, *this << ": replicated vector" << " will be localized to " << *distribution )

        LAMAArray<T> newLocalValues;

        {
            const IndexType newSize = distribution->getLocalSize();

            HostReadAccess<T> rLocalValues( mLocalValues );
            HostWriteOnlyAccess<T> wNewLocalValues( newLocalValues, newSize );

            #pragma omp parallel for
            for ( IndexType i = 0; i < size(); ++i )
            {
                if ( distribution->isLocal( i ) )
                {
                    const IndexType iLocal = distribution->global2local( i );
                    LAMA_ASSERT_DEBUG( iLocal < newSize, "illegal index " << iLocal )
                    wNewLocalValues[iLocal] = rLocalValues[i];
                }
            }
        }

        mLocalValues.swap( newLocalValues );
        setDistributionPtr( distribution );
    }

    else if ( distribution->isReplicated() )
    {
        LAMA_LOG_INFO( logger, *this << " will be replicated" )

        // replicate a distributed vector

        LAMAArray<ValueType> globalValues;

        {
            HostReadAccess<ValueType> localData( mLocalValues );
            HostWriteOnlyAccess<ValueType> globalData( globalValues, size() );
            getDistribution().replicate( globalData.get(), localData.get() );
        }

        mLocalValues.swap( globalValues );
        setDistributionPtr( distribution );
    }
    else
    {
        LAMA_LOG_INFO( logger, *this << " will be redistributed to " << *distribution )

        // so we have now really a redistibution, build a Redistributor

        LAMAArray<ValueType> newLocalValues( distribution->getLocalSize() );

        Redistributor redistributor( distribution, getDistributionPtr() ); // target, source distributions

        redistributor.redistribute( newLocalValues, mLocalValues );

        mLocalValues.swap( newLocalValues );

        setDistributionPtr( distribution );
    }
}

/* ------------------------------------------------------------------------- */

template<typename T>
void DenseVector<T>::resizeImpl()
{
    WriteAccess<T> wLocalValues( mLocalValues );
    wLocalValues.resize( getDistribution().getLocalSize() );
}

/* -- IO ------------------------------------------------------------------- */

template<typename T>
void DenseVector<T>::readVectorHeader( const std::string& filename, File::FileType& fileType, long& dataTypeSize )
{
    LAMA_LOG_INFO( logger, "Read Vector Header" )

    char charFileType;
    std::ifstream inFile( filename.c_str(), std::ios::in );

    if ( !inFile.is_open() )
    {
        LAMA_THROWEXCEPTION( "Unable to open vector header file " + filename + "." )
    }

    IndexType numrows;

    inFile >> charFileType;
    inFile >> numrows;
    inFile >> dataTypeSize;
    inFile.close();

    // Cyclic distribution where first processor gets all

    DistributionPtr distribution( new NoDistribution( numrows ) );
    setDistributionPtr( distribution );

    switch ( charFileType )
    {
    case 'b':
        fileType = File::BINARY;
        break;
    case 'f':
        fileType = File::FORMATTED;
        break;
    case 'x':
        fileType = File::XDR;
        break;
    default:
        LAMA_THROWEXCEPTION( "Invalid header file." )
    }
    LAMA_LOG_TRACE( logger, "Read Vector Header, size = " << size() )
}

template<typename T>
void DenseVector<T>::writeToFile(
    const std::string& fileBaseName,
    const File::FileType fileType/*=XDR*/,
    const File::DataType dataType/*=DOUBLE*/) const
{
    long dataTypeSize;
    std::string file = fileBaseName.c_str();
    file += ".vec";

    dataTypeSize = getDataTypeSize( dataType );

    switch ( fileType )
    {
    case File::FORMATTED:
    {
        writeVectorToFormattedFile( file );
        break;
    }
    case File::BINARY:
    {
        writeVectorToBinaryFile( file, dataTypeSize );
        break;
    }
    case File::XDR:
    {
        writeVectorToXDRFile( file, dataTypeSize );
        break;
    }
    case File::MATRIX_MARKET:
    {
        writeVectorToMMFile( file, dataType );
        break;
    }
    default:
        throw Exception( "Unknown file type definition." );
    } //switch(fileType)

    writeVectorHeader( fileBaseName + ".frv", fileType, dataTypeSize );
}

template<typename T>
void DenseVector<T>::writeVectorHeader(
    const std::string& fileName,
    const File::FileType& fileType,
    const long dataTypeSize ) const
{
    char charFileType;

    switch ( fileType )
    {
    case File::BINARY:
        charFileType = 'b';
        break;
    case File::FORMATTED:
        charFileType = 'f';
        break;
    case File::XDR:
        charFileType = 'x';
        break;
    case File::MATRIX_MARKET:
        // nothing to do for MATRIX_MARKET
        return;
    default:
        LAMA_THROWEXCEPTION( "Invalid header file." )
    }

    std::ofstream outFile( fileName.c_str(), std::ios::out );

    if ( !outFile.is_open() )
    {
        LAMA_THROWEXCEPTION( "Unable to open vector header file " + fileName + "." )
    }

    outFile << charFileType << std::endl;
    outFile << size() << std::endl;
    outFile << dataTypeSize;
    outFile.close();
}

template<typename T>
void DenseVector<T>::writeVectorToMMFile( const std::string& filename, const File::DataType& dataType ) const
{
    MM_typecode veccode;
    mm_initialize_typecode( &veccode );
    mm_set_array( &veccode );
    mm_set_dense( &veccode );
    if ( dataType == File::DOUBLE || dataType == File::FLOAT )
    {
        mm_set_real( &veccode );
    }
    else if ( dataType == File::COMPLEX )
    {
        mm_set_complex( &veccode );
    }
    else if ( dataType == File::INTEGER )
    {
        mm_set_integer( &veccode );
    }
    else if ( dataType == File::PATTERN )
    {
        mm_set_pattern( &veccode );
    }
    else
    {
        LAMA_THROWEXCEPTION( "DenseVector<T>::writeVectorToMMFile: "
                             "Unknown datatype." )
    }

    std::FILE* file;
    if ( !( file = std::fopen( filename.c_str(), "w+" ) ) )
    {
        LAMA_THROWEXCEPTION( "DenseVector<T>::writeVectorToMMFile: '" + filename + "' could not be opened." )
    }

    IndexType numRows = size();

    mm_write_banner( file, veccode );
    mm_write_mtx_array_size( file, numRows, numRows );

    if ( std::fclose( file ) != 0 )
    {
        LAMA_THROWEXCEPTION( "DenseVector<T>::writeVectorToMMFile: '" + filename + "' could not be closed." )
    }
    file = 0;

    std::ofstream ofile;
    ofile.open( filename.c_str(), std::ios::out | std::ios::app );

    if ( ofile.fail() )
    {
        LAMA_THROWEXCEPTION( "DenseVector<T>::writeVectorToMMFile: '" + filename + "' could not be reopened." )
    }

    HostReadAccess<ValueType> dataRead( mLocalValues );

    for ( IndexType ii = 0; ii < numRows; ++ii )
    {
        ofile << ii + 1;
        if ( dataType != File::PATTERN )
        {
            ofile << " " << dataRead[ii];
        }
        ofile << std::endl;
    }
    ofile.close();
}

template<typename T>
long DenseVector<T>::getDataTypeSize( const File::DataType dataType ) const
{
    switch ( dataType )
    {
    case File::DOUBLE:
        return TypeTraits<double>::size;
    case File::FLOAT:
        return TypeTraits<float>::size;
    case File::INTERNAL:
        if ( sizeof(ValueType) == TypeTraits<float>::size )
        {
            return TypeTraits<float>::size;
        }
        else if ( sizeof(ValueType) == TypeTraits<double>::size )
        {
            return TypeTraits<double>::size;
        }
        else
        {
            LAMA_THROWEXCEPTION( "Unknown vector value type size." )
        }
    default:
        LAMA_THROWEXCEPTION( "Unknown vector data type for writing the XDR file." )
    }
}

template<typename T>
void DenseVector<T>::writeVectorToBinaryFile( const std::string& file, const long dataTypeSize ) const
{
    std::fstream outFile( file.c_str(), std::ios::out | std::ios::binary );

    writeVectorDataToBinaryFile( outFile, dataTypeSize );

    outFile.close();
}

template<typename T>
void DenseVector<T>::writeVectorToXDRFile( const std::string& file, const long dataTypeSize ) const
{
    XDRFileStream outFile( file.c_str(), std::ios::out );

    IndexType numRows = size();

    long nnu = static_cast<long>( numRows );
    outFile.write( &nnu );
    outFile.write( &dataTypeSize );

    HostReadAccess<ValueType> dataRead( mLocalValues );

    if ( dataTypeSize == sizeof(ValueType) )
    {
        outFile.write( dataRead.get(), numRows );
    }
    else if ( dataTypeSize == TypeTraits<double>::size )
    {
        double* buffer = new double[numRows];

        for ( IndexType i = 0; i < numRows; i++ )
        {
            buffer[i] = dataRead[i];
        }
        outFile.write( buffer, numRows );
        delete[] buffer;
    }
    else if ( dataTypeSize == TypeTraits<float>::size )
    {
        float* buffer = new float[numRows];

        for ( IndexType i = 0; i < numRows; i++ )
        {
            buffer[i] = static_cast<float>( dataRead[i] );
        }

        outFile.write( buffer, numRows );
        delete[] buffer;
    }

    outFile.write( &nnu );
    outFile.write( &dataTypeSize );
    outFile.close();
}

template<typename FileType,typename DataType>
static void writeBinaryData( std::fstream& outFile, const DataType data[], const IndexType n )
{
    if ( typeid(FileType) == typeid(DataType) )
    {
        // no type conversion needed

        outFile.write( reinterpret_cast<const char*>( data ), sizeof(DataType) * n );
        outFile.flush();
        return;
    }

    // allocate buffer for type conversion

    boost::scoped_array<FileType> buffer( new FileType[n] );

    for ( IndexType i = 0; i < n; i++ )
    {
        buffer[i] = static_cast<FileType>( data[i] );
    }

    outFile.write( reinterpret_cast<const char*>( buffer.get() ), sizeof(FileType) * n );
    outFile.flush();
    return;
}

template<typename T>
void DenseVector<T>::writeVectorDataToBinaryFile( std::fstream& outFile, const long dataTypeSize ) const
{
    IndexType numRows = size();
    HostReadAccess<ValueType> dataRead( mLocalValues );

    if ( sizeof(ValueType) == dataTypeSize )
    {
        writeBinaryData<ValueType,ValueType>( outFile, dataRead.get(), numRows );
    }
    else if ( dataTypeSize == TypeTraits<double>::size )
    {
        writeBinaryData<double,ValueType>( outFile, dataRead.get(), numRows );
    }
    else if ( dataTypeSize == TypeTraits<float>::size )
    {
        writeBinaryData<float,ValueType>( outFile, dataRead.get(), numRows );
    }
    else
    {
        LAMA_THROWEXCEPTION( "unsupported dataTypeSize = " << dataTypeSize )
    }
}

template<typename T>
void DenseVector<T>::writeVectorToFormattedFile( const std::string& file ) const
{
    std::fstream outFile( file.c_str(), std::ios::out );
    HostReadAccess<ValueType> dataRead( mLocalValues );

    for ( IndexType i = 0; i < size(); ++i )
    {
        outFile << dataRead[i] << std::endl;
    }
    outFile.close();
}

template<typename T>
void DenseVector<T>::readVectorFromFormattedFile( const std::string& fileName )
{
    std::ifstream inFile( fileName.c_str(), std::ios::in );

    if ( !inFile.is_open() )
    {
        LAMA_THROWEXCEPTION( "Could not open formatted ascii vector file." )
    }

    const IndexType n = size();

    HostWriteOnlyAccess<T> dataWrite( mLocalValues, n );

    for ( IndexType i = 0; i < n; ++i )
    {
        inFile >> dataWrite[i];
    }
    inFile.close();
    LAMA_LOG_TRACE( logger, "read Vector From Formatted File, " << size() << " values" )
}

template<typename T>
void DenseVector<T>::readVectorFromBinaryFile( const std::string& fileName, const long dataTypeSize )
{
    std::fstream inFile( fileName.c_str(), std::ios::in | std::ios::binary );

    if ( !inFile.is_open() )
    {
        LAMA_THROWEXCEPTION( "Could not open binary vector file." )
    }

//    m_data = m_allocator.allocate(m_capacity);

    readVectorDataFromBinaryFile( inFile, dataTypeSize );

    inFile.close();
}

template<typename T>
void DenseVector<T>::readVectorFromXDRFile( const std::string& fileName, const long dataTypeSizeHeader )
{
    XDRFileStream inFile( fileName.c_str(), std::ios::in | std::ios::binary );

    if ( !inFile.is_open() )
    {
        LAMA_THROWEXCEPTION( "Could not open XDR vector file." )
    }

    // Number of elements
    long nnuLong = 0;
    inFile.read( &nnuLong );

    IndexType nnu = static_cast<IndexType>( nnuLong );
    if ( size() != nnu )
    {
        LAMA_THROWEXCEPTION( "Header file doesn't fit to vector data file. Unequal nnu value." )
    }

    // double or flaot vector data
    long dataTypeSize = 0;
    inFile.read( &dataTypeSize );

    if ( dataTypeSizeHeader != dataTypeSize )
    {
        LAMA_THROWEXCEPTION( "Header file doesn't fit to vector data file. Unequal data type size." )
    }

//    m_data = m_allocator.allocate( m_capacity )

    // read IEEE float
    if ( dataTypeSize == TypeTraits<float>::size )
    {
        if ( sizeof(ValueType) == TypeTraits<float>::size )
        {   //read float

            HostWriteAccess<ValueType> writeData( mLocalValues );

            float *buffer = new float[nnu];
            inFile.read( buffer, nnu );
            for ( IndexType i = 0; i < nnu; i++ )
            {
                writeData[i] = buffer[i];
            }
            delete[] buffer;
        }
        else
        {   //read float and cast to ValueType
            HostWriteAccess<ValueType> writeData( mLocalValues );

            float *buffer = new float[nnu];
            inFile.read( buffer, nnu );
            for ( IndexType i = 0; i < nnu; i++ )
            {
                writeData[i] = static_cast<ValueType>( buffer[i] );
            }
            delete[] buffer;
        }
    }
    // read IEEE double

    else if ( dataTypeSize == TypeTraits<double>::size )
    {
        if ( sizeof(ValueType) == TypeTraits<double>::size )
        {   //read double
            HostWriteAccess<ValueType> writeData( mLocalValues );
            writeData.resize( size() );

            double *buffer = new double[nnu];
            inFile.read( buffer, nnu );
            for ( IndexType i = 0; i < nnu; i++ )
            {
                writeData[i] = static_cast<ValueType>( buffer[i] );
            }
            delete[] buffer;
        }
        else
        {   //read double and cast to ValueType
            HostWriteAccess<ValueType> writeData( mLocalValues );

            double *buffer = new double[nnu];
            inFile.read( buffer, nnu );
            for ( IndexType i = 0; i < nnu; i++ )
            {
                writeData[i] = static_cast<ValueType>( buffer[i] );
            }
            delete[] buffer;
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "Invalid data type size of vector data." )
    }

    // Validate Header
    nnuLong = 0;
    inFile.read( &nnuLong );
    if ( size() != static_cast<IndexType>( nnuLong ) )
    {
        LAMA_THROWEXCEPTION( "Invalid header of the vector file. Unequal nnu." )
    }

    long checkDataType = 0;
    inFile.read( &checkDataType );
    if ( checkDataType != dataTypeSize )
    {
        LAMA_THROWEXCEPTION( "Invalid header of the vector file. Unequal data type size." )
    }
}

/* -------------------------------------------------------------------------- */

template<typename FileType,typename DataType>
static void readBinaryData( std::fstream& inFile, DataType data[], const IndexType n )
{
    if ( typeid(FileType) == typeid(DataType) )
    {
        // no type conversion needed

        inFile.read( reinterpret_cast<char*>( data ), sizeof(DataType) * n );
        return;
    }

    // allocate buffer for type conversion

    boost::scoped_array<FileType> buffer( new FileType[n] );

    inFile.read( reinterpret_cast<char*>( buffer.get() ), sizeof(FileType) * n );

    for ( IndexType i = 0; i < n; i++ )
    {
        data[i] = static_cast<DataType>( buffer[i] );
    }
}

/* -------------------------------------------------------------------------- */

template<typename T>
void DenseVector<T>::readVectorDataFromBinaryFile( std::fstream &inFile, const long dataTypeSize )
{
    LAMA_LOG_INFO( logger, "read vector from binary file, dataTypeSize = " << dataTypeSize )

    IndexType n = size();

    HostWriteOnlyAccess<ValueType> writeData( mLocalValues, n );

    if ( dataTypeSize == sizeof(ValueType) )
    {
        readBinaryData<ValueType,ValueType>( inFile, writeData.get(), n );
    }
    else if ( dataTypeSize == TypeTraits<float>::size )
    {
        readBinaryData<float,ValueType>( inFile, writeData.get(), n );
    }
    else if ( dataTypeSize == TypeTraits<double>::size )
    {
        readBinaryData<double,ValueType>( inFile, writeData.get(), n );
    }
    else
    {
        LAMA_THROWEXCEPTION( "Invalid data type size of vector data: " << dataTypeSize )
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename T>
DenseVector<T>::DenseVector( const DenseVector<T>& other )

    : Vector( other )

{
    // implementation here can be simpler as DenseVector( const Vector& other )

    LAMA_LOG_INFO( logger,
                   "Copy of vector of global size " << size() << ", local size " << getDistribution().getLocalSize() )

    mLocalValues = other.getLocalValues();
}

/* ---------------------------------------------------------------------------------*/


// Template instantiation for all relevant types
template class LAMA_DLL_IMPORTEXPORT DenseVector<float> ;
template class LAMA_DLL_IMPORTEXPORT DenseVector<double> ;

} //namespace
