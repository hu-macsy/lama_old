/**
 * @file solver.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief solver.cpp
 * @author kbuschulte
 * @date 12.09.2012
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
