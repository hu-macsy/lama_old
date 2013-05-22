/**
 * @file solver.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief solver.cpp
 * @author kbuschulte
 * @date 12.09.2012
 * @since 1.0.0
 */

#include <lama.hpp>

// Matrix & vector related includes
#include <lama/DenseVector.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matutils/MatrixCreator.hpp>

// Solver related includes
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/solver/logger/Timer.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/CG.hpp>

using namespace lama;

int main()
{
   typedef double ValueType;

   DenseVector<ValueType> solution( 10000, 2.0 );
   const DenseVector<ValueType> rhs( solution.getDistributionPtr(), 1.0 );

   CSRSparseMatrix<ValueType> matrix;
   MatrixCreator<ValueType>::buildPoisson2D( matrix, 9, 100, 100 );

   LoggerPtr logger(
      new CommonLogger(
         "<CG>: ",
         LogLevel::solverInformation,
         LoggerWriteBehaviour::toConsoleOnly,
         std::auto_ptr<Timer>( new Timer() ) ) );

   CriterionPtr criterion( new IterationCount( 10 ) );

   CG solver( "mySolverName", logger );

   solver.setStoppingCriterion( criterion );
   solver.initialize( matrix );
   solver.solve( solution, rhs );
}
