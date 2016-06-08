/**
 * @file CGNR.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief CGNR.cpp
 * @author David schissler
 * @date 27.05.2015
 */

// hpp
#include <scai/solver/CGNR.hpp>

// local library
#include <scai/solver/mepr/SolverEps.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/DenseVector.hpp>

// std
#include <limits>

namespace scai
{

namespace solver
{

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

SCAI_LOG_DEF_LOGGER( CGNR::logger, "Solver.CGNR" )

CGNR::CGNR( const std::string& id )
    : IterativeSolver( id ) {}


CGNR::CGNR( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id , logger ) {}

CGNR::CGNR( const CGNR& other )
    : IterativeSolver( other ) {}



CGNR::CGNRRuntime::CGNRRuntime()
    : IterativeSolverRuntime() {}

CGNR::~CGNR() {}

CGNR::CGNRRuntime::~CGNRRuntime() {}


void CGNR::initialize( const Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )

    IterativeSolver::initialize( coefficients );
    CGNRRuntime& runtime = getRuntime();

    runtime.mEps = mepr::SolverEps<SCAI_ARITHMETIC_HOST_LIST>::get( coefficients.getValueType() ) * 3.0;

    runtime.mTransposedMat.reset( coefficients.newMatrix() );

    runtime.mVecD.reset( coefficients.newDenseVector() );
    runtime.mVecW.reset( coefficients.newDenseVector() );
    runtime.mVecZ.reset( coefficients.newDenseVector() );
    runtime.mResidual2.reset( coefficients.newDenseVector() );

    runtime.mTransposedMat->assignTranspose( coefficients );
    runtime.mTransposedMat->conj();
}

void CGNR::solveInit( Vector& solution, const Vector& rhs )
{
    CGNRRuntime& runtime = getRuntime();

    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;

    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), rhs.size(), "mismatch: #rows of matrix, rhs" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "mismatch: #cols of matrix, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), rhs.getDistribution(), "mismatch: matrix row dist, rhs dist" )

    // Initialize
    this->getResidual();
    *runtime.mResidual2 = (*runtime.mTransposedMat) * (*runtime.mResidual);

    if(mPreconditioner != NULL)
    {
        *runtime.mVecZ = Scalar(0.0);
        mPreconditioner->solve(*runtime.mVecZ,*runtime.mResidual2);
    }
    else *runtime.mVecZ = *runtime.mResidual2;

    *runtime.mVecD = *runtime.mVecZ;
    runtime.mSolveInit = true;
}

void CGNR::iterate()
{
    CGNRRuntime& runtime = getRuntime();

    const Matrix& A = *runtime.mCoefficients;
    const Matrix& transposedA = *runtime.mTransposedMat;
    Vector& vecW = *runtime.mVecW;
    Vector& vecD = *runtime.mVecD;
    Vector& vecZ = *runtime.mVecZ;
    Vector& residual = *runtime.mResidual;
    Vector& residual2 = *runtime.mResidual2;
    Vector& solution = *runtime.mSolution;
    Scalar alpha;
    Scalar beta;

    lama::L2Norm norm;
    const Scalar& eps = runtime.mEps;

    vecW = A * vecD;
    Scalar normVecW = norm.apply( vecW );
    Scalar scalarProduct = vecZ.dotProduct(residual2);

    if ( normVecW < eps || 1.0/normVecW < eps )           //norm is small
    {
        alpha = 0.0;
    }
    else
    {
        alpha = scalarProduct / ( normVecW * normVecW );
    }

    solution = solution + alpha * vecD;
    residual = residual - alpha * vecW;
    residual2 = transposedA * residual;

    // PRECONDITIONING
    if(mPreconditioner != NULL) mPreconditioner->solve(vecZ,residual2);
    else vecZ = residual2;

    if ( abs(scalarProduct) < eps || 1.0/abs(scalarProduct) < eps )        //norm is small
    {
        beta = 0.0;
    }
    else
    {
        beta = vecZ.dotProduct(residual2) / scalarProduct;
    }

    vecD = vecZ + beta * vecD;
    //CGNR Implementation End
    mCGNRRuntime.mSolution.setDirty( false );
}

SolverPtr CGNR::copy()
{
    return SolverPtr( new CGNR( *this ) );
}

CGNR::CGNRRuntime& CGNR::getRuntime()
{
    return mCGNRRuntime;
}

const CGNR::CGNRRuntime& CGNR::getConstRuntime() const
{
    return mCGNRRuntime;
}

std::string CGNR::createValue()
{
    return "CGNR";
}

Solver* CGNR::create( const std::string name )
{
    return new CGNR( name );
}

void CGNR::writeAt( std::ostream& stream ) const
{
    stream << "CGNR ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

} /* end namespace solver */

} /* end namespace scai */
