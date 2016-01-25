/**
 * @file DenseVector.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @since 1.0.0
 */

// hpp
#include <scai/lama/DenseVector.hpp>

// local library
#include <scai/lama/HArrayUtils.hpp>
#include <scai/lama/UtilKernelTrait.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/lama/LAMAKernel.hpp>

#include <scai/lama/distribution/NoDistribution.hpp>
#include <scai/lama/distribution/CyclicDistribution.hpp>
#include <scai/lama/distribution/Redistributor.hpp>

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/expression/Expression.hpp>
#include <scai/lama/StorageIO.hpp>

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/FileType.hpp>

// internal scai libraries
#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/unique_ptr.hpp>
#include <scai/common/exception/UnsupportedException.hpp>
#include <scai/common/Constants.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <ostream>

namespace scai
{

using common::scoped_array;
using common::TypeTraits;

namespace context = scai::common::context;

using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DenseVector<ValueType>::logger, "Vector.DenseVector" )

template<typename ValueType>
DenseVector<ValueType>::DenseVector()
                : Vector( 0 ), mLocalValues()
{
}

template<typename ValueType>
DenseVector<ValueType>::DenseVector( ContextPtr context )
                : Vector( 0, context ), mLocalValues()
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
    CommunicatorPtr comm = Communicator::get();

    IndexType myRank = comm->getRank();
    IndexType host = 0; // reading processor

    IndexType numElements = 0; // will be the size of the vector

    File::FileType fileType; // = File::XDR;
    std::string suffix;
	std::string baseFileName = filename;
	std::string vecFileName;
    long dataTypeSize = -1;

    if( filename.size() >= 4 )
    {
        suffix = filename.substr( filename.size() - 4, 4 );
    }

    if( suffix == ".frv" )
    {
        baseFileName = filename.substr( 0, filename.size() - 4 );
    }

    if( suffix == ".mtx" )
    {
    	fileType = File::MATRIX_MARKET;
    } else
    {
    	std::string frvFileName = baseFileName + ".frv";

		if( myRank == host )
		{
			readVectorHeader( frvFileName, fileType, dataTypeSize );

			// broadcast the size to all other processors

			numElements = size();
		}

		vecFileName = baseFileName + ".vec";

		comm->bcast( &numElements, 1, host );

		// host gets all elements

		DistributionPtr dist( new CyclicDistribution( numElements, numElements, comm ) );

		allocate( dist );
    }

    if( myRank == host )
    {
        // host owns all elements, all others have no elements

        switch( fileType )
        {
            case File::FORMATTED: //ASCII
                readVectorFromFormattedFile( vecFileName );
                break;

            case File::BINARY: //Binary without number of elements in the vector file
                readVectorFromBinaryFile( vecFileName, getDataType<ValueType>( dataTypeSize ) );
                break;

            case File::XDR: //XDR following the IEEE standard
                readVectorFromXDRFile( vecFileName, dataTypeSize );
                break;
            case File::MATRIX_MARKET:
            	readVectorFromMMFile( filename );
            	break;

            default:
                COMMON_THROWEXCEPTION( "Unknown File Type." );
        }

        // @TODO do not throw exception as application will hang
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const _HArray& localValues, DistributionPtr distribution )
                : Vector( distribution )
{
    SCAI_ASSERT_EQUAL_ERROR( localValues.size(), distribution->getLocalSize() )

    HArrayUtils::assign( mLocalValues, localValues ); // can deal with type conversions
}

/* ------------------------------------------------------------------------- */

/*
 * Constructors with Expressions as arguments
 */

// linear algebra expression: a*x
template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Scalar,Vector,Times>& expression )

                : Vector( expression.getArg2() )
{
    SCAI_LOG_INFO( logger, "Constructor( alpha * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x+b*y, inherit distribution/context from vector x

template<typename ValueType>
DenseVector<ValueType>::DenseVector(
    const Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>& expression ) //Expression_SV_SV

                : Vector( expression.getArg1().getArg2() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x + beta * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector(
    const Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>& expression ) //Expression_SMV_SV

                : Vector( expression.getArg1().getArg2().getArg1().getDistributionPtr(),
                          expression.getArg1().getArg2().getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x+b*y, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector(
    const Expression<Expression<Scalar,Expression<Vector,Matrix,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>& expression ) //Expression_SVM_SV
                : Vector( expression.getArg1().getArg2().getArg2().getColDistributionPtr(),
                          expression.getArg1().getArg2().getArg2().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A + b * y )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*A*x, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& expression )

                : Vector( expression.getArg2().getArg1().getDistributionPtr(),
                          expression.getArg2().getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * A * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: a*x*A, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Scalar,Expression<Vector,Matrix,Times>,Times>& expression )
                : Vector( expression.getArg2().getArg2().getColDistributionPtr(),
                          expression.getArg2().getArg2().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( alpha * x * A )" )
    Vector::operator=( expression );
}

// linear algebra expression: A*x, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Matrix,Vector,Times>& expression )
                : Vector( expression.getArg1().getDistributionPtr(), expression.getArg1().getContextPtr() )
{
    allocate( getDistributionPtr() );
    SCAI_LOG_INFO( logger, "Constructor( A * x )" )
    Vector::operator=( expression );
}

// linear algebra expression: x*A, inherit distribution/context from matrix A

template<typename ValueType>
DenseVector<ValueType>::DenseVector( const Expression<Vector,Matrix,Times>& expression )
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
DenseVector<ValueType>* DenseVector<ValueType>::clone() const
{
    SCAI_LOG_INFO( logger, "DenseVector<ValueType>::clone" )

    DenseVector<ValueType>* newDenseVector = new DenseVector<ValueType>();

    newDenseVector->setContextPtr( mContext );

    return newDenseVector;
}

template<typename ValueType>
DenseVector<ValueType>* DenseVector<ValueType>::clone( DistributionPtr distribution ) const
{
    SCAI_LOG_INFO( logger, "DenseVector<ValueType>::create" )

    DenseVector<ValueType>* newDenseVector = new DenseVector<ValueType>( distribution );

    newDenseVector->setContextPtr( mContext );

    // give back the new vector and its ownership

    return newDenseVector;
}

template<typename ValueType>
DenseVector<ValueType>* DenseVector<ValueType>::copy() const
{
    // create a new dense vector with the copy constructor

    return new DenseVector<ValueType>( *this );
}

template<typename ValueType>
void DenseVector<ValueType>::updateHalo( const Halo& halo ) const
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
tasking::SyncToken* DenseVector<ValueType>::updateHaloAsync( const Halo& halo ) const
{
    const IndexType haloSize = halo.getHaloSize();

    // create correct size of Halo

    {
        WriteOnlyAccess<ValueType> haloAccess( mHaloValues, mContext, haloSize );
    }

    return getDistribution().getCommunicator().updateHaloAsync( mHaloValues, mLocalValues, halo );
}

template<typename ValueType>
Scalar DenseVector<ValueType>::getValue( IndexType globalIndex ) const
{
    SCAI_LOG_TRACE( logger, *this << ": getValue( globalIndex = " << globalIndex << " )" )
    ValueType myValue = static_cast<ValueType>(0.0);
    const IndexType localIndex = getDistribution().global2local( globalIndex );

    if( localIndex != nIndex )
    {
        ContextPtr contextPtr = Context::getContextPtr( context::Host );

        ReadAccess<ValueType> localAccess( mLocalValues, contextPtr );

        SCAI_LOG_TRACE( logger, "index "<< globalIndex << " is local " << localIndex )
        myValue = localAccess[localIndex];
    }

    ValueType allValue = getDistribution().getCommunicator().sum( myValue );
    SCAI_LOG_TRACE( logger, "myValue = " << myValue << ", allValue = " << allValue )
    return Scalar( allValue );
}

template<typename ValueType>
Scalar DenseVector<ValueType>::min() const
{
    //TODO: need a interface function for this

    ContextPtr contextPtr = Context::getContextPtr( context::Host );

    ReadAccess<ValueType> readLocalValues( mLocalValues, contextPtr );
 
    IndexType n = mLocalValues.size();

    const ValueType* localValues = readLocalValues.get();

    ValueType localMin = localValues[0];
#pragma omp parallel
    {
        ValueType myLocalMin = localMin;
#pragma omp for

        for( IndexType i = 0; i < n; ++i )
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

template<typename ValueType>
Scalar DenseVector<ValueType>::max() const
{
    IndexType nnu = mLocalValues.size();

    SCAI_ASSERT_GT( nnu, 0, "no local values for max" )

    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;

    ContextPtr loc = reduce.getValidContext( mLocalValues.getValidContext() );

    ReadAccess<ValueType> localValues( mLocalValues, loc );

    ValueType localMax = reduce[loc]( localValues.get(), localValues.size(), common::reduction::MAX );

    return getDistribution().getCommunicator().max( localMax );
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
    IndexType nnu = mLocalValues.size();

    ValueType localL1Norm = 0;

    if ( nnu > 0 )
    {
        // get available kernel routines for "BLAS1.asum" 

        static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;

        // find valid context, preferred is mContext

        ContextPtr loc = asum.getValidContext( mContext );

        ReadAccess<ValueType> read( mLocalValues, loc );

        SCAI_CONTEXT_ACCESS( loc )

        localL1Norm = asum[loc]( nnu, read.get(), 1 );
    }

    return getDistribution().getCommunicator().sum( localL1Norm );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::l2Norm() const
{
    IndexType nnu = mLocalValues.size();

    ValueType localDotProduct = 0;

    if( nnu > 0 )
    {
        // get available kernel routines for "BLAS1.dot" 

        static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;

        // find valid context, preferred is mContext

        ContextPtr loc = dot.getValidContext( mContext );

        ReadAccess<ValueType> read( mLocalValues, loc );

        SCAI_CONTEXT_ACCESS( loc )

        localDotProduct = dot[loc]( nnu, read.get(), 1, read.get(), 1 );
    }

    ValueType globalDotProduct = getDistribution().getCommunicator().sum( localDotProduct );

    return sqrt( globalDotProduct );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseVector<ValueType>::maxNorm() const
{
    IndexType nnu = mLocalValues.size(); // number of local rows

    ValueType localMaxNorm = static_cast<ValueType>(0.0);

    if ( nnu > 0 )
    {
        static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;

        ContextPtr loc = reduce.getValidContext( mContext );  

        ReadAccess<ValueType> read( mLocalValues, loc );

        SCAI_CONTEXT_ACCESS( loc )

        localMaxNorm = reduce[loc]( read.get(), nnu, common::reduction::ABS_MAX );
    }

    const Communicator& comm = getDistribution().getCommunicator();

    ValueType globalMaxNorm = comm.max( localMaxNorm );

    SCAI_LOG_INFO( logger,
                   comm << ": max norm " << *this << ", local max norm of " << nnu << " elements: " << localMaxNorm 
                   << ", max norm global = " << globalMaxNorm )

    return globalMaxNorm;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::swap( Vector& other )
{
    SCAI_LOG_DEBUG( logger, "swap:" << *this << " with " << other )

    DenseVector* otherPtr = dynamic_cast<DenseVector*>( &other );

    if( !otherPtr )
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
void DenseVector<ValueType>::vectorPlusVector(
    ContextPtr prefContext,
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y )
{
    SCAI_LOG_DEBUG( logger,
                    "vectorPlusVector: result:" << result << " = " << alpha << " * x:" << x << " + " << beta << " * y:" << y )

    // get function pointers, do not use fallbacks here

    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
    static LAMAKernel<blaskernel::BLASKernelTrait::axpy<ValueType> > axpy;
    static LAMAKernel<blaskernel::BLASKernelTrait::sum<ValueType> > sum;

    ContextPtr context = sum.getValidContext( setVal, axpy, prefContext );

    const IndexType nnu = result.size();

    if( &result == &x && &result == &y ) //result = alpha * result + beta * result
    {
        //result = alpha * result + beta * result
        //=>
        //result = ( alpha + beta ) * result
        //=>
        //result *= ( alpha + beta )

        SCAI_LOG_DEBUG( logger,
                        "vectorPlusVector: x = y = result, result *= " << "alpha(" << alpha << ") + beta(" << beta << ")" )

        WriteAccess<ValueType> resultAccess( result, context, true );

        SCAI_CONTEXT_ACCESS( context )
        setVal[context]( resultAccess.get(), nnu, alpha + beta, common::reduction::MULT );
    }
    else if( &result == &x ) //result = alpha * result + beta * y
    {
        ReadAccess<ValueType> yAccess( y, context );
        WriteAccess<ValueType> resultAccess( result, context, true );

        if( beta == scai::common::constants::ZERO )
        {
            SCAI_LOG_DEBUG( logger, "vectorPlusVector: result *= alpha" )

            if( alpha != scai::common::constants::ONE ) // result *= alpha
            {
                SCAI_CONTEXT_ACCESS( context )
                setVal[context]( resultAccess.get(), nnu, alpha, common::reduction::MULT );
            }
            else
            {
                // do nothing: result = 1 * result
            }
        }
        else if( beta == scai::common::constants::ONE ) // result = alpha * result + y
        {
            SCAI_LOG_DEBUG( logger, "vectorPlusVector: result = alpha * result + y" )

            if( alpha != scai::common::constants::ONE ) // result = alpha * result + y
            {
                // result *= alpha
                SCAI_CONTEXT_ACCESS( context )
                setVal[context]( resultAccess.get(), nnu, alpha, common::reduction::MULT );
            }

            // result += y
            SCAI_CONTEXT_ACCESS( context )
            axpy[context]( nnu, static_cast<ValueType>(1.0)/*alpha*/, yAccess.get(), 1, resultAccess.get(), 1 );
        }
        else // beta != 1.0 && beta != 0.0 --> result = alpha * result + beta * y
        {
            SCAI_LOG_DEBUG( logger,
                            "vectorPlusVector: result = alpha(" << alpha << ")" << " * result + beta(" << beta << ") * y" )

            if( alpha != scai::common::constants::ONE )
            {
                SCAI_CONTEXT_ACCESS( context )
                setVal[context]( resultAccess.get(), nnu, alpha, common::reduction::MULT );
            }

            SCAI_CONTEXT_ACCESS( context )
            axpy[context]( nnu, beta, yAccess.get(), 1, resultAccess.get(), 1 );
        }
    }
    else if( &result == &y ) // result = alpha * x + beta * result
    {
        SCAI_LOG_DEBUG( logger,
                        "vectorPlusVector: result = alpha(" << alpha << ")" << " * x + beta(" << beta << ") * result" )

        // so we do here:  result = beta * result, result += alpha * x

        ReadAccess<ValueType> xAccess( x, context );
        WriteAccess<ValueType> resultAccess( result, context, true );

        if( beta != scai::common::constants::ONE ) // result = [alpha * x + ] beta * result
        {
            // result *= beta
            SCAI_CONTEXT_ACCESS( context )
            setVal[context]( resultAccess.get(), nnu, beta, common::reduction::MULT );
        }

        if( alpha != scai::common::constants::ZERO )
        {
            // result = alpha * x + result
            SCAI_CONTEXT_ACCESS( context )
            axpy[context]( nnu, alpha, xAccess.get(), 1, resultAccess.get(), 1 );
        }
    }
    else // result = alpha * x + beta * y
    {
        SCAI_LOG_DEBUG( logger,
                        "vectorPlusVector: result = alpha(" << alpha << ")" << " * x + beta(" << beta << ") * y" )

        ReadAccess<ValueType> xAccess( x, context );
        ReadAccess<ValueType> yAccess( y, context );
        // no need to keep old values of result
        WriteAccess<ValueType> resultAccess( result, context, false );

        SCAI_CONTEXT_ACCESS( context )
        sum[context]( nnu, alpha, xAccess.get(), beta, yAccess.get(), resultAccess.get() );
    }

    SCAI_LOG_INFO( logger, "vectorPlusVector done" )
}

template<typename ValueType>
tasking::SyncToken* DenseVector<ValueType>::vectorPlusVectorAsync(
    ContextPtr /*context*/,
    HArray<ValueType>& /*result*/,
    const ValueType /*alpha*/,
    const HArray<ValueType>& /*x*/,
    const ValueType /*beta*/,
    const HArray<ValueType>& /*y*/)
{
    COMMON_THROWEXCEPTION( "vectorPlusVectorAsync not implemented yet" )
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

    if( x.getDistribution() != y.getDistribution() )
    {
        COMMON_THROWEXCEPTION(
                        "distribution do not match for z = alpha * x + beta * y, z = "<< *this <<" , x = "<< x <<" , y = "<< y )
    }

    if( x.getDistribution() != getDistribution() || x.size() != size() )
    {
        resize( x.getDistributionPtr() );
    }

    if( typeid( *this ) == typeid( x ) && typeid( *this ) == typeid( y ) )
    {
        const DenseVector<ValueType>& denseX = dynamic_cast<const DenseVector<ValueType>&>( x );
        const DenseVector<ValueType>& denseY = dynamic_cast<const DenseVector<ValueType>&>( y );

        if( mLocalValues.size() != denseX.mLocalValues.size() )
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

        SCAI_LOG_DEBUG( logger, "call vectorPlusVector" )

        vectorPlusVector( mContext, mLocalValues, alpha, denseX.mLocalValues, beta, denseY.mLocalValues );
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

    if( this->getValueType() == other.getValueType() )
    {
        if( getDistribution() != other.getDistribution() )
        {
            COMMON_THROWEXCEPTION( "distribution do not match for this * other, this = "<< *this <<" , other = "<< other )
        }

        const DenseVector<ValueType>* denseOther = dynamic_cast<const DenseVector<ValueType>*>( &other );

        SCAI_ASSERT_DEBUG( denseOther, "dynamic_cast failed for other = " << other )

        SCAI_LOG_DEBUG( logger, "Calculating local dot product at " << *mContext )

        // get available kernel routines for "BLAS1.dot" 

        static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;

        // find valid context, preferred is mContext

        ContextPtr loc = dot.getValidContext( mContext );

        // Now do the dot production at location loc ( might have been changed to other location  )

        ReadAccess<ValueType> localRead( mLocalValues, loc );
        ReadAccess<ValueType> otherRead( denseOther->mLocalValues, loc );

        SCAI_CONTEXT_ACCESS( loc )

        const IndexType localSize = mLocalValues.size();

        SCAI_ASSERT_EQUAL_DEBUG( localSize, getDistribution().getLocalSize() )

        const ValueType localDotProduct = dot[loc]( localSize, localRead.get(), 1, otherRead.get(), 1 );

        SCAI_LOG_DEBUG( logger, "Calculating global dot product form local dot product = " << localDotProduct )

        ValueType dotProduct = getDistribution().getCommunicator().sum( localDotProduct );

        SCAI_LOG_DEBUG( logger, "Global dot product = " << dotProduct )

        return dotProduct;
    }

    COMMON_THROWEXCEPTION(
                    "Can not calculate a dot product of "<< typeid( *this ).name() <<" and "<< typeid( other ).name() )
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

    HArrayUtils::setScalar( mLocalValues, value, common::reduction::COPY, mContext );
}

template<typename ValueType>
void DenseVector<ValueType>::assign( const _HArray& localValues, DistributionPtr dist )
{
    SCAI_LOG_INFO( logger, "assign vector with localValues = " << localValues << ", dist = " << *dist )

    SCAI_ASSERT_EQUAL_ERROR( localValues.size(), dist->getLocalSize() )

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
    const IndexType size = mLocalValues.size();

    static LAMAKernel<UtilKernelTrait::invert<ValueType> > invert;

    const ContextPtr loc = invert.getValidContext( this->getContextPtr() );

    SCAI_CONTEXT_ACCESS( loc );

    WriteAccess<ValueType> wValues( mLocalValues, loc );

    invert[loc]( wValues.get(), size );
}

template<typename ValueType>
size_t DenseVector<ValueType>::getMemoryUsage() const
{
    // Note: memory of mHaloValues is not counted, is just a temporary

    size_t memoryUsage = sizeof(ValueType) * mLocalValues.size();
    return getDistribution().getCommunicator().sum( memoryUsage );
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseVector<ValueType>::redistribute( DistributionPtr distribution )
{
    SCAI_ASSERT_EQUAL_ERROR( size(), distribution->getGlobalSize() )

    if( getDistribution() == *distribution )
    {
        SCAI_LOG_INFO( logger, *this << " redistribute to same distribution " << *distribution )

        // we can keep local/global values, but just set dist pointer

        setDistributionPtr( distribution );
    }

    else if( getDistribution().isReplicated() )

    {
        SCAI_LOG_INFO( logger, *this << ": replicated vector" << " will be localized to " << *distribution )

        HArray<ValueType> newLocalValues;

        ContextPtr hostContext = Context::getContextPtr( context::Host );

        {
            const IndexType newSize = distribution->getLocalSize();

            ReadAccess<ValueType> rLocalValues( mLocalValues, hostContext );
            WriteOnlyAccess<ValueType> wNewLocalValues( newLocalValues, hostContext, newSize );

#pragma omp parallel for

            for( IndexType i = 0; i < size(); ++i )
            {
                if( distribution->isLocal( i ) )
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

    else if( distribution->isReplicated() )
    {
        SCAI_LOG_INFO( logger, *this << " will be replicated" )

        // replicate a distributed vector

        HArray<ValueType> globalValues;

        ContextPtr hostContext = Context::getContextPtr( context::Host );

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

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::resizeImpl()
{
    // resize array with local values

    mLocalValues.resize( getDistribution().getLocalSize() );
}

/* -- IO ------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::readVectorHeader(
    const std::string& filename,
    File::FileType& fileType,
    long& dataTypeSize )
{
    SCAI_LOG_INFO( logger, "Read Vector Header" )

    char charFileType;
    std::ifstream inFile( filename.c_str(), std::ios::in );

    if( !inFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Unable to open vector header file " + filename + "." )
    }

    IndexType numrows;

    inFile >> charFileType;
    inFile >> numrows;
    inFile >> dataTypeSize;
    inFile.close();

    // Cyclic distribution where first processor gets all

    DistributionPtr distribution( new NoDistribution( numrows ) );
    setDistributionPtr( distribution );

    switch( charFileType )
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
            COMMON_THROWEXCEPTION( "Invalid header file." )
    }

    SCAI_LOG_TRACE( logger, "Read Vector Header, size = " << size() )
}

template<typename ValueType>
void DenseVector<ValueType>::writeToFile(
    const std::string& fileBaseName,
    const File::FileType fileType/*=XDR*/,
    const common::scalar::ScalarType dataType/*=DOUBLE*/) const
{
    std::string file = fileBaseName.c_str();
    file += ".vec";

    switch( fileType )
    {
        case File::FORMATTED:
        {
            writeVectorToFormattedFile( file );
            break;
        }

        case File::BINARY:
        {
            writeVectorToBinaryFile( file, dataType );
            break;
        }

        case File::XDR:
        {
            writeVectorToXDRFile( file, dataType );
            break;
        }

        case File::MATRIX_MARKET:
        {
            writeVectorToMMFile( fileBaseName, dataType );
            break;
        }

        default:
            COMMON_THROWEXCEPTION( "Unknown file type definition." );
    } //switch(fileType)

    long dataTypeSize = getDataTypeSize<ValueType>( dataType );

    if( fileType != File::MATRIX_MARKET )
    {
		writeVectorHeader( fileBaseName + ".frv", fileType, dataTypeSize );
    }
}

template<typename ValueType>
void DenseVector<ValueType>::writeVectorHeader(
    const std::string& fileName,
    const File::FileType& fileType,
    const long dataTypeSize ) const
{
    char charFileType;

    switch( fileType )
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
            COMMON_THROWEXCEPTION( "Invalid header file." )
    }

    std::ofstream outFile( fileName.c_str(), std::ios::out );

    if( !outFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Unable to open vector header file " + fileName + "." )
    }

    outFile << charFileType << std::endl;
    outFile << size() << std::endl;
    outFile << dataTypeSize;
    outFile.close();
}

template<typename ValueType>
void DenseVector<ValueType>::writeVectorToMMFile( const std::string& filename, const common::scalar::ScalarType& dataType ) const
{
	IndexType numRows = size();

	std::string fullFilename = filename;
 
    if ( !_StorageIO::hasSuffix( filename, ".mtx" ) )
    {
        fullFilename += ".mtx";
    }

	_StorageIO::writeMMHeader( true, numRows, 1, numRows, fullFilename, dataType );

    std::ofstream ofile;
    ofile.open( fullFilename.c_str(), std::ios::out | std::ios::app );

    if( ofile.fail() )
    {
        COMMON_THROWEXCEPTION( "DenseVector<ValueType>::writeVectorToMMFile: '" + fullFilename + "' could not be reopened." )
    }

    ContextPtr hostContext = Context::getContextPtr( context::Host );

    ReadAccess<ValueType> dataRead( mLocalValues, hostContext );

    for( IndexType ii = 0; ii < numRows; ++ii )
    {

        if( dataType == common::scalar::PATTERN )
        {
            ofile << ii + 1;
        }
        else
        {
            ofile << " " << dataRead[ii];
        }

        ofile << "\n"; //std::endl;
    }

    ofile.close();
}

template<typename ValueType>
void DenseVector<ValueType>::writeVectorToBinaryFile( const std::string& file, const common::scalar::ScalarType type ) const
{
    std::fstream outFile( file.c_str(), std::ios::out | std::ios::binary );

    writeVectorDataToBinaryFile( outFile, type );

    outFile.close();
}

template<typename FileType,typename DataType>
static void writeDataToXDRFile( XDRFileStream& outFile, const DataType* data, const IndexType n )
{
    if( typeid(FileType) == typeid(DataType) )
    {
        outFile.write( data, n ); // no conversion needed
        return;
    }

    // so user data has to be converted in file type data

    scoped_array<FileType> buffer( new FileType[n] );

    for( IndexType i = 0; i < n; i++ )
    {
        buffer[i] = static_cast<FileType>( data[i] );
    }

    outFile.write( buffer.get(), n );
}

template<typename ValueType>
void DenseVector<ValueType>::writeVectorToXDRFile( const std::string& file, const common::scalar::ScalarType dataType ) const
{
    XDRFileStream outFile( file.c_str(), std::ios::out );

    IndexType numRows = size();

    long nnu = static_cast<long>( numRows );

    long dataTypeSize = getDataTypeSize<ValueType>( dataType );

    outFile.write( &nnu );
    outFile.write( &dataTypeSize );

    ContextPtr hostContext = Context::getContextPtr( context::Host );

    ReadAccess<ValueType> dataRead( mLocalValues, hostContext );

    switch ( dataType )
    {

#define LAMA_WRITE_XDR_CASE( z, I, _ )                                                                     \
        case common::TypeTraits<ARITHMETIC_HOST_TYPE_##I>::stype :                                         \
            writeDataToXDRFile<ARITHMETIC_HOST_TYPE_##I, ValueType>( outFile, dataRead.get(), numRows );   \
            break;                                                                                         \

        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_WRITE_XDR_CASE, _ )

#undef LAMA_WRITE_XDR_CASE

        case common::scalar::INTERNAL:
            writeDataToXDRFile<ValueType, ValueType>( outFile, dataRead.get(), numRows );
            break;

        case common::scalar::PATTERN:
            break;

        default:
            COMMON_THROWEXCEPTION( "unsupported file data type = " << dataType )
    }

    outFile.write( &nnu );
    outFile.write( &dataTypeSize );
    outFile.close();
}

template<typename FileType,typename DataType>
static void writeBinaryData( std::fstream& outFile, const DataType data[], const IndexType n )
{
    if( typeid(FileType) == typeid(DataType) )
    {
        // no type conversion needed

        outFile.write( reinterpret_cast<const char*>( data ), sizeof(DataType) * n );
        outFile.flush();
        return;
    }

    // allocate buffer for type conversion

    scoped_array<FileType> buffer( new FileType[n] );

    for( IndexType i = 0; i < n; i++ )
    {
        buffer[i] = static_cast<FileType>( data[i] );
    }

    outFile.write( reinterpret_cast<const char*>( buffer.get() ), sizeof(FileType) * n );
    outFile.flush();
    return;
}

template<typename ValueType>
void DenseVector<ValueType>::writeVectorDataToBinaryFile( std::fstream& outFile, const common::scalar::ScalarType type ) const
{
    IndexType numRows = size();

    ContextPtr contextPtr = Context::getContextPtr( context::Host );

    ReadAccess<ValueType> dataRead( mLocalValues, contextPtr );

    switch( type )
    {

#define LAMA_WRITE_CASE( z, I, _ )                                                                      \
        case common::TypeTraits<ARITHMETIC_HOST_TYPE_##I>::stype :                                      \
            writeBinaryData<ARITHMETIC_HOST_TYPE_##I, ValueType>( outFile, dataRead.get(), numRows );   \
            break;                                                                                      \

        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_WRITE_CASE, _ )

#undef LAMA_WRITE_CASE

        case common::scalar::INTERNAL:
            writeBinaryData<ValueType,ValueType>( outFile, dataRead.get(), numRows );
            break;

        default:
            COMMON_THROWEXCEPTION( "unsupported file data type = " << type )
    }
}

template<typename ValueType>
void DenseVector<ValueType>::writeVectorToFormattedFile( const std::string& file ) const
{
    ContextPtr hostContext = Context::getContextPtr( context::Host );

    std::fstream outFile( file.c_str(), std::ios::out );

    ReadAccess<ValueType> dataRead( mLocalValues, hostContext );

    for( IndexType i = 0; i < size(); ++i )
    {
        outFile << dataRead[i] << std::endl;
    }

    outFile.close();
}

template<typename ValueType>
void DenseVector<ValueType>::readVectorFromFormattedFile( const std::string& fileName )
{
    ContextPtr hostContext = Context::getContextPtr( context::Host );

    std::ifstream inFile( fileName.c_str(), std::ios::in );

    if( !inFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Could not open formatted ascii vector file." )
    }

    const IndexType n = size();

    WriteOnlyAccess<ValueType> dataWrite( mLocalValues, hostContext, n );

    for( IndexType i = 0; i < n; ++i )
    {
        inFile >> dataWrite[i];
    }

    inFile.close();
    SCAI_LOG_TRACE( logger, "read Vector From Formatted File, " << size() << " values" )
}

template<typename ValueType>
void DenseVector<ValueType>::readVectorFromBinaryFile( const std::string& fileName, const common::scalar::ScalarType type )
{
    std::fstream inFile( fileName.c_str(), std::ios::in | std::ios::binary );

    if( !inFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Could not open binary vector file " << fileName << "." )
    }

    readVectorDataFromBinaryFile( inFile, type );

    inFile.close();
}

template<typename FileDataType, typename UserDataType>
static void readXDRData( XDRFileStream& inFile, UserDataType data[], const IndexType n )
{
    if( typeid(FileDataType) == typeid(UserDataType) )
    {
        // no type conversion needed

        inFile.read( data, n );
        return;
    }

    // allocate a temporary buffer for n values of FileDataType to read the data

    scoped_array<FileDataType> buffer( new FileDataType[n] );

    inFile.read( buffer.get(), n );

    for( IndexType i = 0; i < n; i++ )
    {
        data[i] = static_cast<UserDataType>( buffer[i] );
    }
}

template<typename ValueType>
void DenseVector<ValueType>::readVectorFromMMFile( const std::string& fileName )
{
    bool isSymmetric, isPattern;
    IndexType numRows, numColumns, numValues;

    _StorageIO::readMMHeader( numRows, numColumns, numValues, isPattern, isSymmetric, fileName );

    SCAI_ASSERT_EQUAL_ERROR( numColumns, 1 )

    std::ifstream ifile;
    ifile.open( fileName.c_str(), std::ios::in );

    if( ifile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not reopen file '" << fileName << "'." )
    }

    CommunicatorPtr comm = Communicator::get();
    DistributionPtr dist( new CyclicDistribution( numRows, numRows, comm ) );

    allocate( dist );

    // First reading in the beginning of the rows
    // then reading in the values and columns of the rows
    //Jump to the beginning of the Values
    char c = '%';

    while( c == '%' )
    {
        ifile >> c;
        ifile.ignore( 1024, '\n' );
    }

    IndexType lines = numValues;

    IndexType i;
    ValueType val;

    std::string line;

    WriteOnlyAccess<ValueType> vector( mLocalValues, numValues );
    ValueType* vPtr = vector.get();

    for( int l = 0; l < lines && !ifile.eof(); ++l )
    {
        std::getline( ifile, line );
        std::istringstream reader( line );

        if( isPattern )
        {
            reader >> i;
            val = 1.0;
            i--;
        }
        else
        {
            reader >> val;
            i = l;
        }

        vPtr[i] = val;

    }

    if( ifile.eof() )
    {
    	COMMON_THROWEXCEPTION( "'" << fileName << "': reached end of file, before having read all data." )
    }

    ifile.close();
    ifile.close();
    SCAI_LOG_INFO( logger, "construct vector " << numRows )
}

template<typename ValueType>
void DenseVector<ValueType>::readVectorFromXDRFile( const std::string& fileName, const long dataTypeSizeHeader )
{
    XDRFileStream inFile( fileName.c_str(), std::ios::in | std::ios::binary );

    if( !inFile.is_open() )
    {
        COMMON_THROWEXCEPTION( "Could not open XDR vector file." )
    }

    // Number of elements
    long nnuLong = 0;
    inFile.read( &nnuLong );

    IndexType nnu = static_cast<IndexType>( nnuLong );

    if( size() != nnu )
    {
        COMMON_THROWEXCEPTION( "Header file doesn't fit to vector data file. Unequal nnu value." )
    }

    // double or flaot vector data
    long dataTypeSize = 0;
    inFile.read( &dataTypeSize );

    if( dataTypeSizeHeader != dataTypeSize )
    {
        COMMON_THROWEXCEPTION( "Header file doesn't fit to vector data file. Unequal data type size." )
    }

    // Attention: determination of file type by size is ambiguous, e.g. Complex and Double
    //            have same size. If ambiguous, own ValueType is the preferred one

    common::scalar::ScalarType fileType = getDataType<ValueType>( dataTypeSize );

    WriteOnlyAccess<ValueType> writeData( mLocalValues, nnu );

    switch( fileType )
    {
        case common::scalar::INTERNAL:
            readXDRData<ValueType,ValueType>( inFile, writeData.get(), nnu );
            break;

#define LAMA_READ_XDR_CASE( z, I, _ )                                                          \
        case common::TypeTraits<ARITHMETIC_HOST_TYPE_##I>::stype:                              \
            readXDRData<ARITHMETIC_HOST_TYPE_##I, ValueType>( inFile, writeData.get(), nnu );  \
            break;                                                                             \

        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_READ_XDR_CASE, _ )

#undef LAMA_READ_XDR_CASE

        case common::scalar::PATTERN:
            // that might be okay
            break;

        default:
            COMMON_THROWEXCEPTION( "Invalid data type size of vector data." )
    }

    // Validate Header

    nnuLong = 0;
    inFile.read( &nnuLong );

    if( size() != static_cast<IndexType>( nnuLong ) )
    {
        COMMON_THROWEXCEPTION( "Invalid header of the vector file. Unequal nnu." )
    }

    long checkDataType = 0;
    inFile.read( &checkDataType );

    if( checkDataType != dataTypeSize )
    {
        COMMON_THROWEXCEPTION( "Invalid header of the vector file. Unequal data type size." )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseVector<ValueType>::readVectorDataFromBinaryFile( std::fstream &inFile, const common::scalar::ScalarType type )
{
    IndexType n = size();

    SCAI_LOG_INFO( logger,
                   "read DenseVector<" << TypeTraits<ValueType>::id() << "> from binary file, size = " << n << ", dataType = " << ( ( common::scalar::ScalarType ) type ) )

    WriteOnlyAccess<ValueType> writeData( mLocalValues, n );

    switch( type )
    {

#define LAMA_READ_BIN_CASE( z, I, _ )                                                                   \
        case common::TypeTraits<ARITHMETIC_HOST_TYPE_##I>::stype :                                      \
            FileIO::readBinaryData<ARITHMETIC_HOST_TYPE_##I, ValueType>( inFile, writeData.get(), n );  \
            break;                                                                                      \

        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_READ_BIN_CASE, _ )

#undef LAMA_READ_BIN_CASE

        case common::scalar::INTERNAL:
            FileIO::readBinaryData<ValueType,ValueType>( inFile, writeData.get(), n );
            break;

        default:
            COMMON_THROWEXCEPTION( "unsupported file data type = " << type )
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Vector* DenseVector<ValueType>::create()
{
    return new DenseVector<ValueType>();
}

template<typename ValueType>
std::pair<VectorKind, common::scalar::ScalarType> DenseVector<ValueType>::createValue()
{
    common::scalar::ScalarType skind = common::getScalarType<ValueType>();
    return std::pair<VectorKind, common::scalar::ScalarType> ( DENSE, skind );
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

#define LAMA_DENSE_VECTOR_INSTANTIATE(z, I, _)                                     \
    template class COMMON_DLL_IMPORTEXPORT DenseVector<ARITHMETIC_HOST_TYPE_##I> ;

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_DENSE_VECTOR_INSTANTIATE, _ )

#undef LAMA_DENSE_VECTOR_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
