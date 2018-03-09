/**
 * @file CentralPathHessian.hpp
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
 * @brief CentralPathHessian matrix as used in LSQ with box boundary conditions
 * @author Thomas Brandes, Andreas Borgen Langva
 * @date 21.07.2017
 */

#pragma once

#include <scai/lama/matrix/OperatorMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>

namespace scai

{

namespace lama

{

/** The CentralPathHessian matrix stands for H = 2 * tau transpose( A ) * A + D 
 *
 *  This matrix is used for a solver that computes the search direction in each 
 *  iteration step of an iterative method in least square with box constraints.
 *
 *   - While the matrix A remains unchanged, the values of tau and D
 *     will be updated in each iteration step
 *
 */
template<typename ValueType>
class CentralPathHessian : public OperatorMatrix<ValueType>
{

public:

    CentralPathHessian ( const Matrix<ValueType>& A ) :

        OperatorMatrix<ValueType>( A.getColDistributionPtr(), A.getColDistributionPtr() ),

        mA( A ),
        mD( NULL ),
        mTau( 0 )
    {
        // create a tmp vector that contains result A * x

        tmpTargetPtr.reset( A.newTargetVector() );
        tmpSourcePtr.reset( A.newSourceVector() );
    }

    ~CentralPathHessian()
    {
    }

    void update( const Vector<ValueType>& D, const ValueType tau )
    {
        SCAI_ASSERT_EQ_ERROR( D.size(), this->getNumColumns(), "illegal size for vector D" );

        mD = &D;
        mTau = tau;
    }

    virtual void matrixTimesVector(
        Vector<ValueType>& result,
        const ValueType alpha,
        const Vector<ValueType>& x,
        const ValueType beta,
        const Vector<ValueType>* y,
        const common::MatrixOp op ) const
    {
        SCAI_ASSERT_EQ_ERROR( op, common::MatrixOp::NORMAL, "no conj/transpose supported" )

        Vector<ValueType>& tmpSource = *tmpSourcePtr;
        Vector<ValueType>& tmpTarget = *tmpTargetPtr;

        tmpSource = alpha * *mD * x; // elmentwise multiplication

        if ( y != NULL )
        {
            result = beta * *y;
            result += tmpSource;
        }
        else 
        {
            result = tmpSource;
        }

        tmpTarget = mA * x;
        result += ( 2 * mTau * alpha ) * transpose( mA ) * tmpTarget;
    }
    
    virtual void matrixTimesVectorDense(
        DenseVector<ValueType>&,
        const ValueType,
        const DenseVector<ValueType>&,
        const ValueType,
        const DenseVector<ValueType>*,
        const common::MatrixOp ) const
    {
        COMMON_THROWEXCEPTION( "should not be called" )
    }

    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA.getContextPtr();
    }   
    
    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "CentralPathHessian( 2 * " << mTau << " * A' * A + D" << " )";
    }

private:

    using _Matrix::logger;

    const Matrix<ValueType>& mA;   // reference kept the whole lifetime of this object
    const Vector<ValueType>* mD;   // changing reference to the actual vector D

    ValueType mTau;        // value tau

    // avoid reallocation of temporary vector required in matrixTimesVector

    std::unique_ptr<Vector<ValueType>> tmpTargetPtr;
    std::unique_ptr<Vector<ValueType>> tmpSourcePtr;
};

}

}
