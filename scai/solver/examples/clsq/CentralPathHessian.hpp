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

#include <scai/lama/matrix/AbstractMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

namespace scai

{

namespace lama

{

/** The CentralPathHessian matrix stands for H = 2 * tau A' * A + D 
 *
 *  This matrix is used for a solver that computes the search direction in each 
 *  iteration step of an iterative method in least square with box constraints.
 *
 *   - While the matrix A remains unchanged, the values of tau and D
 *     will be updated in each iteration step
 *
 */
class CentralPathHessian : public AbstractMatrix 
{

public:

    CentralPathHessian ( const CSRSparseMatrix<double>& A ) :

        AbstractMatrix( A.getColDistributionPtr(), A.getColDistributionPtr() ),

        mA( A ),
        mAT( NULL ),
        mD( NULL ),
        mTau( 0 )
    {
    }

    void setTransposed( const CSRSparseMatrix<double>& AT )
    {
        mAT = &AT;
    }

    ~CentralPathHessian()
    {
    }

    void update( const Vector& D, const Scalar tau )
    {
        SCAI_ASSERT_EQ_ERROR( D.size(), getNumColumns(), "illegal size for vector D" );

        mD = &D;
        mTau = tau;
    }

    virtual void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        SCAI_LOG_INFO( logger, "matrixTimesVector, mA = " << mA )

        result = *mD * x;

        result += beta * y;

        DenseVector<double> tmp( mA * x );

        if ( mAT )
        {
            result += ( 2 * mTau * alpha ) * (*mAT) * tmp;
        }
        else
        {
            result += ( 2 * mTau * alpha ) * tmp * mA;
        }
    }
    
    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA.getContextPtr();
    }   
    
    /** This method must be provided so that solvers can decide about the type of additional runtime vectors. */

    virtual common::ScalarType getValueType() const
    {
        return mA.getValueType();
    }   

private:

    const CSRSparseMatrix<double>& mA;
    const CSRSparseMatrix<double>* mAT;
    const Vector* mD;
    Scalar mTau;
};

}

}
