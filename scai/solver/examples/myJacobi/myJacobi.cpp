/**
 * @file MyJacobiTest.cpp
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
 * @brief Contains the implementation of the class MyJacobi.cpp
 * @author: Alexander Büchel, Matthias Makulla
 * @date 22.02.2012
 * @since 1.0.0
 **/

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
   typedef double ValueType;

    int dim = 3;
    int stencil = 7;
    int size = 100;

    CSRSparseMatrix<ValueType> matrix;
    MatrixCreator<double>::buildPoisson( matrix, dim, stencil, size, size, size );

    int vectorSize = static_cast<int>( std::pow(size,dim) );

    DenseVector<ValueType> exactSolution( vectorSize, 1.0 );
    DenseVector<ValueType> rhs = matrix * exactSolution;

    DenseVector<ValueType> solution( vectorSize, 0.0 );
    
    LoggerPtr slogger( new CommonLogger( "MyJacobiLogger:", LogLevel::convergenceHistory, LoggerWriteBehaviour::toConsoleOnly) );
    MyJacobi jacobiSolver( "MyJacobi", slogger );
    jacobiSolver.initialize( matrix );

    CriterionPtr criterion( new IterationCount( 100 ) );
    jacobiSolver.setStoppingCriterion( criterion );
    jacobiSolver.solve( solution, rhs );

    DenseVector<ValueType> diff( solution - exactSolution );

    L2Norm l2Norm;
    Scalar norm = l2Norm( diff );
}
