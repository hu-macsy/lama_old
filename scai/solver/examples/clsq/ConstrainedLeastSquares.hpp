/**
 * @file ConstrainedLeastSquares.hpp
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
