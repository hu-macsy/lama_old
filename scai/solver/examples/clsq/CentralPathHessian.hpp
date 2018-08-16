/**
 * @file CentralPathHessian.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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
