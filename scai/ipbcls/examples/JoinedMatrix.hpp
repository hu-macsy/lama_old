/**
 * @file JoinedMatrix.hpp
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
 * @brief Abstract Matrix class that encapsulates A and an explicit transposed A
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#pragma once

#include <scai/lama/matrix/OperatorMatrix.hpp>
#include "JoinedVector.hpp"

namespace scai

{

namespace lama

{

/** An object of this class encapsulates two matrices [ A1; A2 }
 */

template<typename ValueType>
class JoinedMatrix : public OperatorMatrix<ValueType>
{

public:

    JoinedMatrix ( const Matrix<ValueType>& A1, const Matrix<ValueType>& A2 ) :

        OperatorMatrix<ValueType>( dmemo::DistributionPtr( new dmemo::JoinedDistribution( A1.getRowDistributionPtr(),
                                   A2.getRowDistributionPtr() ) ),
                                   A1.getColDistributionPtr() ),
        mA1( A1 ),
        mA2( A2 ),
        joinedDim( 0 )
    {
        SCAI_ASSERT_EQ_ERROR( A1.getNumColumns(), A2.getNumColumns(), "joined matrices must have same number of columns" )
        SCAI_ASSERT_EQ_ERROR( A1.getValueType(), A2.getValueType(), "joined matrices must have same value type" )
    }

    ~JoinedMatrix()
    {
    }

    virtual void matrixTimesVector(
        Vector<ValueType>& result,
        const ValueType alpha,
        const Vector<ValueType>& x,
        const ValueType beta,
        const Vector<ValueType>* y,
        const common::MatrixOp op ) const
    {
        // SCAI_ASSERT result is JoinedVector
        // SCAI_ASSERT if beta != 0 y is JoinedVector

        SCAI_LOG_INFO( logger, "joined matrix, result = " << alpha << " * A * x + " << beta << " * y " )

        if ( !common::isTranspose( op ) )
        {
            SCAI_ASSERT_EQ_ERROR( result.getVectorKind(), VectorKind::JOINED, "joined vector expected" )
            JoinedVector<ValueType>& jResult = reinterpret_cast<JoinedVector<ValueType>&>( result );

            if ( y == NULL )
            {
                mA1.matrixTimesVector( jResult.first(), alpha, x, beta, NULL, op );
                mA2.matrixTimesVector( jResult.second(), alpha, x, beta, NULL, op );
            }
            else
            {
                SCAI_ASSERT_EQ_ERROR( y->getVectorKind(), VectorKind::JOINED, "joined vector expected" )
                const JoinedVector<ValueType>* jY = reinterpret_cast<const JoinedVector<ValueType>*>( y );

                mA1.matrixTimesVector( jResult.first(), alpha, x, beta, &jY->first(), op );
                mA2.matrixTimesVector( jResult.second(), alpha, x, beta, &jY->second(), op );
            }
        }
        else
        {
            SCAI_ASSERT_EQ_ERROR( x.getVectorKind(), VectorKind::JOINED, "joined vector expected" )
            const JoinedVector<ValueType>& jX = reinterpret_cast<const JoinedVector<ValueType>&>( x );

            if ( y == NULL )
            {
                mA1.matrixTimesVector( result, alpha, jX.first(), beta, NULL, op );
            }
            else
            {
                mA1.matrixTimesVector( result, alpha, jX.first(), beta, y, op );
            }

            mA2.matrixTimesVector( result, alpha, jX.second(), 1, &result, op );
        }
    }

    /** Implementation of pure method Matrix<ValueType>::matrixTimesVectorDense
     *
     *  This method will never be called as we have already overridden matrixTimesVectorDense
     */
    virtual void matrixTimesVectorDense(
        DenseVector<ValueType>&,
        const ValueType,
        const DenseVector<ValueType>&,
        const ValueType,
        const DenseVector<ValueType>*,
        const common::MatrixOp ) const
    {
        COMMON_THROWEXCEPTION( "Joined matrix does not support matrixTimesVectorDense" )
    }

    virtual void reduce(
        Vector<ValueType>& v,
        const IndexType dim,
        const common::BinaryOp reduceOp,
        const common::UnaryOp elemOp ) const
    {
        if ( dim == 0 )
        {
            COMMON_THROWEXCEPTION( "unsupported" )
        }
        else
        {
            mA1.reduce( v, dim, reduceOp, elemOp );
            DenseVector<ValueType> tmpV;
            mA2.reduce( tmpV, dim, reduceOp, elemOp );
            v += tmpV;
        }
    }

    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA1.getContextPtr();
    }

    /** This method must be provided so that solvers can decide about the type of additional runtime vectors. */

    virtual common::ScalarType getValueType() const
    {
        return mA1.getValueType();
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "JoinedMatrix( " << mA1 << ", " << mA2 << " )";
    }

    /** Joined matrix has joined row distribution, so we need a new method for temporary row vector. */

    virtual Vector<ValueType>* newTargetVector() const
    {
        return new JoinedVector<ValueType>( mA1.newTargetVector(), mA2.newTargetVector() );
    }

private:

    using Matrix<ValueType>::logger;

    const Matrix<ValueType>& mA1;
    const Matrix<ValueType>& mA2;

    IndexType joinedDim;
};

}

}
