/**
 * @file solver/examples/myJacobi/myJacobi.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Contains the implementation of the class MyJacobi.
 * @author Matthias Makulla
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
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <cmath>

using namespace scai;
using namespace solver;
using namespace lama;
using namespace hmemo;

int main( int , char** )
{
    typedef DefaultReal ValueType;
    int dim = 3;
    int stencil = 7;
    int size = 100;
    CSRSparseMatrix<ValueType> matrix;
    MatrixCreator::buildPoisson( matrix, dim, stencil, size, size, size );
    int vectorSize = static_cast<int>( std::pow( size, dim ) );
    const auto exactSolution = fill<DenseVector<ValueType>>( vectorSize, 1 );
    const auto rhs = eval<DenseVector<ValueType>>( matrix * exactSolution );
    auto solution = fill<DenseVector<ValueType>>( vectorSize, 0 );
    LoggerPtr slogger( new CommonLogger( "MyJacobiLogger:", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly ) );
    MyJacobi<ValueType> jacobiSolver( "MyJacobi", slogger );
    jacobiSolver.initialize( matrix );
    CriterionPtr<ValueType> criterion( new IterationCount<ValueType>( 100 ) );
    jacobiSolver.setStoppingCriterion( criterion );
    jacobiSolver.solve( solution, rhs );
    const auto diff = eval<DenseVector<ValueType>>( solution - exactSolution );
    L2Norm<ValueType> l2Norm;
    ValueType norm = l2Norm( diff );
    std::cout << "Final norm = " << norm << std::endl;
}
