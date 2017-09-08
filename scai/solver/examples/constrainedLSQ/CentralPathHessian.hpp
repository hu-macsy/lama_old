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

    CentralPathHessian ( const Matrix& A ) :

        AbstractMatrix( A.getColDistributionPtr(), A.getColDistributionPtr() ),

        mA( A ),
        mD( NULL ),
        mTau( 0 )
    {
        // create a tmp vector that contains result A * x

        tmpVectorPtr.reset( A.newVector( 0 ) );
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
        Vector& tmp = *tmpVectorPtr;

        result = *mD * x;          // elmentwise multiplication
        result += beta * y;
        tmp = mA * x + 0 * tmp;
        result += ( 2 * mTau * alpha ) * tmp * mA;
    }
    
    /** This method must be provided so that solvers can decide about context of operations. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mA.getContextPtr();
    }   
    
    /** This method must be provided so that solvers can decide about the type of additional runtime vectors. */

    virtual common::scalar::ScalarType getValueType() const
    {
        return mA.getValueType();
    }   

private:

    const Matrix& mA;   // reference kept the whole lifetime of this object
    const Vector* mD;   // reference to the actual vector D
    Scalar mTau;        // value tau

    VectorPtr tmpVectorPtr;   // avoid reallocation of temporary vector required in matrixTimesVector
};

}

}
