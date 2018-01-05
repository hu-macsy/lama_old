/**
 * @file MatrixWithT.hpp
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
 * @brief Operator matrix class that encapsulates A and an explicit (conj) transposed A
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#pragma once

#include <scai/lama/matrix/OperatorMatrix.hpp>
#include <scai/common/ScalarType.hpp>

namespace scai

{

namespace lama

{

/** An object of this class encapsulates a matrix A where the (conj) transposed matrix
 *  is constructed explicitly.
 *
 *  The disadvantage of additional memory allocation and the construction of the
 *  transposed matrix is mostly compensated by the benefits of much faster 
 *  multiplication A' * x.
 */

template<typename ValueType>
class MatrixWithT : public OperatorMatrix<ValueType>
{

public:

    /**
     *  Encapsulate matrices explicilty.
     *
     *  @param A  is matrix that is used in all operations alpha * M * x
     *  @param AT is matrix that is used in all operations alpha * x * M 
     */
    MatrixWithT ( const Matrix<ValueType>& A, const Matrix<ValueType>& AT ) :

        OperatorMatrix<ValueType>( A.getRowDistributionPtr(), A.getColDistributionPtr() ),
        mA( A ), 
        mAT( AT )

    {
    }

    /**
     *  Encapsulate matrices explicilty.
     *
     *  @param A  is matrix that is used in all operations alpha * A * x
     *
     *  Builds explicitly the matrix AT as AT = transpose(A) or AT = conjTranspose(A)
     */
    MatrixWithT ( const Matrix<ValueType>& A, bool conjFlag = true ) :

        OperatorMatrix<ValueType>( A.getRowDistributionPtr(), A.getColDistributionPtr() ),
        mATPtr( A.newMatrix() ),
        mA( A ),
        mAT( *mATPtr )

    {
        SCAI_ASSERT_EQ_ERROR( A.getCommunicationKind(), mAT.getCommunicationKind(), "bad new matrix" )
        SCAI_ASSERT_EQ_ERROR( A.getContextPtr(), mAT.getContextPtr(), "bad new matrix" )
        mATPtr->assignTranspose( A );
        if ( conjFlag && common::isComplex( A.getValueType() ) )
        {
            mATPtr->conj();
        }
    }

    ~MatrixWithT()
    {
    }

    virtual void matrixTimesVectorImpl(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>& y ) const
    {
        mA.matrixTimesVector( result, alpha, x, beta, y );
    }
    
    virtual void vectorTimesMatrixImpl(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>& y ) const
    {
        mAT.matrixTimesVector( result, alpha, x, beta, y );
    }
    
    virtual void reduce(
        _Vector& v,
        const IndexType dim,
        const common::BinaryOp reduceOp,
        const common::UnaryOp elemOp ) const
    {
        if ( dim == 1 )
        {
            // reduce is much faster on the transposed matrix

            mAT.reduce( v, 0, reduceOp, elemOp );
        }
        else
        {
            mA.reduce( v, dim, reduceOp, elemOp );
        }
    }

    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA.getContextPtr();
    }   
    
private:

    MatrixPtr<ValueType> mATPtr;   // transposed matrix might also be allocated here

    const Matrix<ValueType>& mA;
    const Matrix<ValueType>& mAT;
};

}

}
