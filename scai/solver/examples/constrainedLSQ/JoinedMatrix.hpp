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

#include <scai/lama/matrix/AbstractMatrix.hpp>
#include "JoinedVector.hpp"

namespace scai

{

namespace lama

{

/** An object of this class encapsulates two matrices [ A1; A2 }
 */

class JoinedMatrix : public AbstractMatrix 
{

public:

    JoinedMatrix ( const Matrix& A1, const Matrix& A2 ) :

        AbstractMatrix( dmemo::DistributionPtr( new dmemo::JoinedDistribution( A1.getRowDistributionPtr(),
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

    /** Override default implementation of Matrix as we need here a joined vector */

    virtual Vector* newVector( const IndexType dim ) const
    {
        if ( dim == joinedDim )
        {
            return new JoinedVector( mA1.newVector( joinedDim ), mA2.newVector( joinedDim ) );
        }
        else
        {
            return mA1.newVector( dim );
        }
    }

    virtual void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        // SCAI_ASSERT result is JoinedVector
        // SCAI_ASSERT if beta != 0 y is JoinedVector

        SCAI_LOG_INFO( logger, "result = " << alpha << " * A * x + " << beta << " * y " )
    
        SCAI_ASSERT_EQ_ERROR( result.getVectorKind(), Vector::JOINED, "joined vector expected" )
        JoinedVector& jResult = reinterpret_cast<JoinedVector&>( result );

        SCAI_ASSERT_EQ_ERROR( y.getVectorKind(), Vector::JOINED, "joined vector expected" )
        const JoinedVector& jY = reinterpret_cast<const JoinedVector&>( y );

        mA1.matrixTimesVector( jResult.first(), alpha, x, beta, jY.first() );
        mA2.matrixTimesVector( jResult.second(), alpha, x, beta, jY.second() );
    }
    
    virtual void vectorTimesMatrix(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        SCAI_ASSERT_EQ_ERROR( x.getDistribution(), getRowDistribution(), "distribution mismatch" )

        SCAI_ASSERT_EQ_ERROR( x.getVectorKind(), Vector::JOINED, "joined vector expected" )
        const JoinedVector& jX = reinterpret_cast<const JoinedVector&>( x );

        mA1.vectorTimesMatrix( result, alpha, jX.first(), beta, y );
        mA2.vectorTimesMatrix( result, alpha, jX.second(), 1, result );
    }
    
    virtual void reduce(
        Vector& v,
        const IndexType dim,
        const common::binary::BinaryOp reduceOp,
        const common::unary::UnaryOp elemOp ) const
    {
        if ( dim == 0 )
        {
            COMMON_THROWEXCEPTION( "unsupported" )
        }
        else
        {
            mA1.reduce( v, dim, reduceOp, elemOp );
            VectorPtr tmpV( v.copy() );
            mA2.reduce( *tmpV, dim, reduceOp, elemOp );
            v += *tmpV;
        }
    }

    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA1.getContextPtr();
    }   
    
    /** This method must be provided so that solvers can decide about the type of additional runtime vectors. */

    virtual common::scalar::ScalarType getValueType() const
    {
        return mA1.getValueType();
    }   

private:

    const Matrix& mA1;
    const Matrix& mA2;

    IndexType joinedDim;
};

}

}
