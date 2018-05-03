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

#include <scai/solver/logger/CommonLogger.hpp>

#include "CentralPathHessian.hpp"

#include <memory>

namespace scai

{

namespace solver

{

class ConstrainedLeastSquares

{

public:

    ConstrainedLeastSquares( const lama::CSRSparseMatrix<double>& A );

    void useTranspose();

    void setTolerance( double eps );

    void setMaxIter( IndexType maxIter );

    void solve ( lama::DenseVector<double>& x, 
                 const lama::DenseVector<double>& b,
                 const lama::DenseVector<double>& lb,
                 const lama::DenseVector<double>& ub );

private:

    void dualityGap(
        double& gap,
        double& dualObj,
        const lama::DenseVector<double>& b,
        const lama::DenseVector<double>& x,
        const lama::DenseVector<double>& lb,
        const lama::DenseVector<double>& ub );

    void computeSearchDirection(
        lama::DenseVector<double>& dx,
        const lama::DenseVector<double>& b,
        const lama::DenseVector<double>& x,
        const double tau,
        const lama::DenseVector<double>& lb,
        const lama::DenseVector<double>& ub,
        const double gap );

    double stepSize(
        const lama::DenseVector<double>& b,
        const lama::DenseVector<double>& x,
        const double tau,
        const lama::DenseVector<double>& lb,
        const lama::DenseVector<double>& ub,
        const lama::DenseVector<double>& dx,
        const double alpha,
        const double beta );

    double centralPathObjective (
        const lama::DenseVector<double>& b,
        const lama::DenseVector<double>& x,
        const double tau,
        const lama::DenseVector<double>& lb,
        const lama::DenseVector<double>& ub );
    
    // member variables of this class

    const lama::CSRSparseMatrix<double>& mA;

    lama::CentralPathHessian mCentralPathHessian;

    std::unique_ptr<lama::CSRSparseMatrix<double> > mAT;

    double mTolerance;

    double mMaxIter;

    LoggerPtr mSolverLogger;   // used as logger during the solver phase

    // Runtime data, used during solve step

    lama::DenseVector<double> mDiagATA;

    lama::CSRSparseMatrix<double> mDiagonalMatrix;

    // static logger for this class

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

}

}
