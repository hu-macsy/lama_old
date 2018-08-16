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
#include <scai/lama/Vector.hpp>

#include <scai/common/NonCopyable.hpp>

#include <scai/solver/logger/CommonLogger.hpp>

#include <memory>

enum class InnerSolverType
{
    StandardCG,
    NewtonStepCG
};

namespace scai
{

namespace ipbcls
{

template<typename ValueType>
class ConstrainedLeastSquares : private common::NonCopyable
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
    ConstrainedLeastSquares( const lama::Matrix<ValueType>& A );

    /**
     * Set the tolerance for the relative objective value of the optimization problem.
     *
     * Specifically, this means that the algorithm terminates if it finds an approximation to the
     * solution of the bound-constrained least squares problem that satisfies
     *
     *     (|| r || - || r* ||) / || r* || <= tol,
     *
     * where || r* || denotes the exact minimum (Euclidean) norm for the least squares problem and r denotes the residual of the
     * current approximation.
     *
     * Note that the algorithm terminates if *either* the objective tolerance or the residual tolerance is satisfied.
     */
    void setObjectiveTolerance( ValueType tol );

    /**
     * Set the tolerance for the relative residual norm.
     *
     * Specifically, this means that the algorithm terminates if it finds an approximation to the
     * solution of the bound-constrained least squares problem that satisfies
     *
     *     ||r|| <= tol * ( ||A|| * ||x|| + ||b|| )
     *
     * where r = Ax - b and ||A|| is a lower bound for the matrix 2-norm of A, which corresponds to the
     * largest singular value of A. A practical estimate that is commonly used is the
     * largest absolute element of A.
     *
     * By default, the estimate for ||A|| is set to zero, which normally works well for well-behaved problems,
     * but for very ill-conditioned A it may be difficult to satisfy.
     *
     * For most problems, the objective tolerance is far more practical. However, for problems where
     * the exact minimum ||r*|| is zero or very close to zero, the objective tolerance criterion breaks down.
     * In these cases, the residual tolerance is preferred.
     *
     * Note that the algorithm terminates if *either* the objective tolerance or the residual tolerance is satisfied.
     */
    void setResidualTolerance( ValueType tol, ValueType matrixNormLowerBound = ValueType( 0 ) );

    /**
     * Returns the objective tolerance.
     *
     * See setObjectiveTolerance( ValueType ) for details.
     *
     */
    ValueType getObjectiveTolerance() const;

    /**
     * Returns the residual tolerance.
     *
     * See setResidualTolerance( ValueType, ValueType ) for details.
     *
     */

    ValueType getResidualTolerance() const;

    void setInnerSolverType ( InnerSolverType type );

    void setMaxIter( IndexType maxIter );

    void solve ( lama::Vector<ValueType>& x,
                 const lama::Vector<ValueType>& b,
                 const lama::Vector<ValueType>& lb,
                 const lama::Vector<ValueType>& ub ) const;

private:

    const lama::Matrix<ValueType>& mA;

    ValueType mObjTolerance;
    ValueType mResTolerance;
    ValueType mMatrixNormLowerBound;
    IndexType mMaxIter;
    InnerSolverType mInnerSolverType;

    // static logger for this class

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

}

}
