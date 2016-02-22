/**
 * @file solver.cpp
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
 * @brief solver.cpp
 * @author kbuschulte
 * @date 12.09.2012
 * @since 1.0.0
 */

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

// Solver related includes
#include <scai/lama/solver/logger/CommonLogger.hpp>
#include <scai/lama/solver/logger/Timer.hpp>
#include <scai/lama/solver/criteria/IterationCount.hpp>
#include <scai/lama/solver/CG.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai::lama;

int main()
{
   //
   // Define the ValueType used for the vector
   //
   typedef double ValueType;

   //
   // Create DenseVectors for solution and right hand side
   //
   DenseVector<ValueType> solution( 10000, 2.0 );
   const DenseVector<ValueType> rhs( solution.getDistributionPtr(), 1.0 );

   //
   // Create a 2D Poisson Matrix with 9Points and dimension 100 in every direction
   //
   CSRSparseMatrix<ValueType> matrix;
   MatrixCreator<ValueType>::buildPoisson2D( matrix, 9, 100, 100 );

   //
   // Set up CG solver
   //

   // Logger (optional) prints out some solver information to console
   LoggerPtr logger(
      new CommonLogger(
         "<CG>: ",
         LogLevel::solverInformation,
         LoggerWriteBehaviour::toConsoleOnly,
         scai::common::shared_ptr<Timer>( new Timer() ) ) );

   // stopping criterion for solver: stopps after 10 iterations
   CriterionPtr criterion( new IterationCount( 10 ) );

   // solver itself
   CG solver( "mySolverName", logger );

   // initialization
   solver.setStoppingCriterion( criterion );
   solver.initialize( matrix );

   // solution phase
   solver.solve( solution, rhs );

   solution.writeToFile( "solution", File::FORMATTED );

   std::cout << "Solution vector is written to 'solution.frm/.vec'" << std::endl;

   //
   //  That's it.
   //
   std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

   return EXIT_SUCCESS;
}
