/**
 * @file solver/examples/myJacobi/myJacobi.cpp
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
 * @brief Contains the implementation of the class MyJacobi.cpp
 * @author Alexander Büchel, Matthias Makulla
 * @date 22.02.2012
 */

#include "MyJacobiModule.hpp"

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/norm/L2Norm.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>

#include <cmath>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

int main( int , char** )
{
    typedef RealType ValueType;
    int dim = 3;
    int stencil = 7;
    int size = 100;
    CSRSparseMatrix<ValueType> matrix;
    MatrixCreator::buildPoisson( matrix, dim, stencil, size, size, size );
    int vectorSize = static_cast<int>( std::pow( size, dim ) );
    DenseVector<ValueType> exactSolution( vectorSize, 1.0 );
    DenseVector<ValueType> rhs = matrix * exactSolution;
    DenseVector<ValueType> solution( vectorSize, 0.0 );
    LoggerPtr slogger( new CommonLogger( "MyJacobiLogger:", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly ) );
    MyJacobi jacobiSolver( "MyJacobi", slogger );
    jacobiSolver.initialize( matrix );
    CriterionPtr criterion( new IterationCount( 100 ) );
    jacobiSolver.setStoppingCriterion( criterion );
    jacobiSolver.solve( solution, rhs );
    DenseVector<ValueType> diff( solution - exactSolution );
    L2Norm l2Norm;
    Scalar norm = l2Norm( diff );
}
