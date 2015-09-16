/**
 * @file CRTPMatrix.hpp
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
 * @brief Template class for common methods of SparseMatrix<ValueType> and DenseMatrix<ValueType>
 *        to deal with polymorphism
 * @author Thomas Brandes
 * @date 09.08.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/Matrix.hpp>

// local library
#include <scai/lama/expression/all.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/Assert.hpp>

namespace scai
{

namespace lama
{

/** This template class supports static polymorphism to define
 *  common routines for base classes SparseMatrix<ValueType> and DenseMatrix<ValueType>.
 *
 *  Therefore it uses the Curiously Recurring Template Pattern (CRTP)
 *  as a C++ idiom where a class X derived from this template is
 *  itself a template argument.
 *
 *  @tparam Derived is the derived class of Matrix
 *  @tparam ValueType specifies the type of the matrix values
 *
 *  @todo: create, copy should also be defined here
 */

template<class Derived,typename ValueType>
class CRTPMatrix: public Matrix
{
public:

    /** Default constructor. */

    CRTPMatrix()
        : Matrix()
    {
    }

    /** Constructor for a replicated zero-matrix of a given shape. */

    CRTPMatrix( const IndexType numRows, const IndexType numColumns )
        : Matrix( numRows, numColumns )
    {
    }

    /** Constructor for a distributed matrix. */

    CRTPMatrix( DistributionPtr rowDistribution, DistributionPtr colDistribution )
        : Matrix( rowDistribution, colDistribution )
    {
    }

    void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        SCAI_REGION( "Mat.timesVector" )

        SCAI_LOG_INFO( logger, result << " = " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y )

        if( &result == &y )
        {
            SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
        }
        else if( &result == &x )
        {
            COMMON_THROWEXCEPTION( "alias: result = x is not handled, use temporary" )
        }
        else
        {
            // we inherit the row distribution of this matrix to result

            result.resize( getDistributionPtr() );

            // no more to check: result.size() == mNumRows, getDistirubtion() == result.getDistribution()
        }

        SCAI_ASSERT_EQUAL_ERROR( x.getDistribution(), getColDistribution() )
        SCAI_ASSERT_EQUAL_ERROR( y.getDistribution(), getDistribution() )

        const DenseVector<ValueType>* denseX = dynamic_cast<const DenseVector<ValueType>*>( &x );
        const DenseVector<ValueType>* denseY = dynamic_cast<const DenseVector<ValueType>*>( &y );
        DenseVector<ValueType>* denseResult = dynamic_cast<DenseVector<ValueType>*>( &result );

        SCAI_ASSERT( denseX, x << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

        // Note: in case of beta == 0, we might skip this test

        SCAI_ASSERT( denseY, y << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

        SCAI_ASSERT( denseResult, result << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

        static_cast<const Derived*>( this )->matrixTimesVectorImpl( *denseResult, alpha.getValue<ValueType>(), *denseX,
                beta.getValue<ValueType>(), *denseY );
    }

    void vectorTimesMatrix(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        SCAI_REGION( "Mat.vectorTimes" )

        SCAI_LOG_INFO( logger, result << " = " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y )

        if( &result == &y )
        {
            SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
        }
        else if( &result == &x )
        {
            COMMON_THROWEXCEPTION( "alias: result = x is not handled, use temporary" )
        }
        else
        {
            // we inherit the row distribution of this matrix to result

            result.resize( getDistributionPtr() );

            // no more to check: result.size() == mNumRows, getDistirubtion() == result.getDistribution()
        }

        SCAI_ASSERT_EQUAL_ERROR( x.getDistribution(), getDistribution() )
        SCAI_ASSERT_EQUAL_ERROR( y.getDistribution(), getColDistribution() )

        const DenseVector<ValueType>* denseX = dynamic_cast<const DenseVector<ValueType>*>( &x );
        const DenseVector<ValueType>* denseY = dynamic_cast<const DenseVector<ValueType>*>( &y );
        DenseVector<ValueType>* denseResult = dynamic_cast<DenseVector<ValueType>*>( &result );

        SCAI_ASSERT( denseX, x << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

        // Note: in case of beta == 0, we might skip this test

        SCAI_ASSERT( denseY, y << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

        SCAI_ASSERT( denseResult, result << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

        static_cast<const Derived*>( this )->vectorTimesMatrixImpl( *denseResult, alpha.getValue<ValueType>(), *denseX,
                beta.getValue<ValueType>(), *denseY );
    }

    using Matrix::setIdentity;
    using Matrix::operator=;

protected:
#ifndef SCAI_LOG_LEVEL_OFF
    using Matrix::logger;
#endif
};

} /* end namespace lama */

} /* end namespace scai */
