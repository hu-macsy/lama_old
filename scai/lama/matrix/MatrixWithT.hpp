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
 * @brief Operator matrix class that encapsulates A and an explicit transposed A
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#pragma once

#include <scai/lama/matrix/OperatorMatrix.hpp>

namespace scai

{

namespace lama

{

/** An object of this class encapsulates a matrix A where the transposed matrix
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
     *  @brief Construct the operator matrix by providing A and At expliclitly.
     *
     *  @param[in]  A is the original matrix
     *  @param[in]  AT is transpose(A)
     *  
     */
    MatrixWithT ( const Matrix<ValueType>& A, const Matrix<ValueType>& AT ) :

        OperatorMatrix<ValueType>( A.getRowDistributionPtr(), A.getColDistributionPtr() ),
        mA( A ), 
        mAT( AT )
    {
        // we cannot verify that AT is really the transpose of A, but we check at least distributions

        SCAI_ASSERT_EQ_ERROR(  A.getRowDistribution(), AT.getColDistribution(), "dist mismatch" )
        SCAI_ASSERT_EQ_ERROR(  A.getColDistribution(), AT.getRowDistribution(), "dist mismatch" )
    }

    /** 
     *  @brief Construct the operator matrix by providing A only
     *
     *  @param[in]  A is the original matrix
     *
     *  This constructor builds the tranposed matrix of A explicitly.
     */
    MatrixWithT ( const Matrix<ValueType>& A ) :

        OperatorMatrix<ValueType>( A.getRowDistributionPtr(), A.getColDistributionPtr() ),
        mATPtr( A.newMatrix() ),
        mA( A ),
        mAT( *mATPtr )

    {
        SCAI_ASSERT_EQ_ERROR( A.getCommunicationKind(), mAT.getCommunicationKind(), "bad new matrix" )
        SCAI_ASSERT_EQ_ERROR( A.getContextPtr(), mAT.getContextPtr(), "bad new matrix" )
        mATPtr->assignTranspose( A );
    }

    ~MatrixWithT()
    {
    }

    /**
     *  Implementation of pure method Matrix<ValueType>::matrixTimesVectorDense
     */
    virtual void matrixTimesVectorDense(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>* y,
        const common::MatrixOp op ) const
    {
        if ( common::MatrixOp::NORMAL == op )
        {
            // no tranpose, use orginal matrix
            mA.matrixTimesVector( result, alpha, x, beta, y, op );
        }
        else if ( common::MatrixOp::TRANSPOSE == op )
        {
            // use transposed matrix and unset transposeFlag
            mAT.matrixTimesVector( result, alpha, x, beta, y, common::MatrixOp::NORMAL );
        }
        else
        {
            COMMON_THROWEXCEPTION( "matrix op " << op << " not supported yet" )
        }
    }
    
    /**
     *  Implement Matrix<ValueType>::reduce with more efficient solution.
     */
    virtual void reduce(
        Vector<ValueType>& v,
        const IndexType dim,
        const common::BinaryOp reduceOp,
        const common::UnaryOp elemOp ) const
    {
        if ( dim == 1 )
        {
            // reduce along columns is much faster on the transposed matrix

            mAT.reduce( v, 0, reduceOp, elemOp );
        }
        else
        {
            // reduce along rows efficient as entries are grouped by rows

            mA.reduce( v, dim, reduceOp, elemOp );
        }
    }

    virtual const MatrixStorage<ValueType>& getLocalStorage() const
    {
        return mA.getLocalStorage();
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
