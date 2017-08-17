/**
 * @file ConstrainedLeastSquares.hpp
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
 * @brief Solver for Least Square with Boundary Conditions
 * @author Thomas Brandes, Andreas Borgen Langva
 * @date 21.07.2017
 */

#pragma once

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/common/unique_ptr.hpp>

#include <scai/solver/logger/CommonLogger.hpp>

#include "CentralPathHessian.hpp" 

namespace scai

{

namespace solver

{

class ConstrainedLeastSquares

{

public:

    /** Constructor of a solver
     *
     *  - allocates all runtime data
     *  
     *  @param[in] A is the matrix for which a least square problem is to be solved.
     *
     *  ATTENTION: The solver object only keeps a reference to the matrix A, so it
     *             must not be destroyed during the lifetime of this solver object
     */
    ConstrainedLeastSquares( const lama::Matrix& A );

    void useTranspose();

    void setTolerance( lama::Scalar eps );

    void setMaxIter( IndexType maxIter );

    void solve ( lama::Vector& x, 
                 const lama::Vector& b,
                 const lama::Vector& lb,
                 const lama::Vector& ub );

private:

    void dualityGap(
        lama::Scalar& gap,
        lama::Scalar& dualObj,
        const lama::Vector& b,
        const lama::Vector& x,
        const lama::Vector& lb,
        const lama::Vector& ub );

    void computeSearchDirection(
        lama::Vector& dx,
        const lama::Vector& b,
        const lama::Vector& x,
        const lama::Scalar tau,
        const lama::Vector& lb,
        const lama::Vector& ub,
        const lama::Scalar gap );

    lama::Scalar stepSize(
        const lama::Vector& b,
        const lama::Vector& x,
        const lama::Scalar tau,
        const lama::Vector& lb,
        const lama::Vector& ub,
        const lama::Vector& dx,
        const lama::Scalar alpha,
        const lama::Scalar beta );

    lama::Scalar centralPathObjective (
        const lama::Vector& b,
        const lama::Vector& x,
        const lama::Scalar tau,
        const lama::Vector& lb,
        const lama::Vector& ub );
    
    // member variables of this class

    const lama::Matrix& mA;   

    lama::CentralPathHessian mCentralPathHessian;

    lama::Scalar mTolerance;

    IndexType mMaxIter;

    LoggerPtr mSolverLogger;   // used as logger during the solver phase

    // Runtime data, used during solve step

    lama::DenseVector<double> mDiagATA;

    lama::CSRSparseMatrix<double> mDiagonalMatrix;

    // static logger for this class

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    lama::VectorPtr resPtr;   // keep temporary for residual, A * x - b
    lama::VectorPtr kappaPtr;  // 

    lama::VectorPtr d1Ptr;   // keep temporary for x - lb
    lama::VectorPtr d2Ptr;   // keep temporary for ub - x
    lama::VectorPtr dxPtr;   // keep temporary for ub - x

    lama::VectorPtr dPtr;   
    lama::VectorPtr DPtr;   
    lama::VectorPtr gPtr;   
    lama::VectorPtr mPtr;   
    lama::VectorPtr xNewPtr;   

};

}

}
