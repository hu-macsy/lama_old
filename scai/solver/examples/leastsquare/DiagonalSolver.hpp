/**
 * @file DiagonalSolver.hpp
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
 * @brief DiagonalSolver is a solver where the matrix has only a diagonal.
 * @author Thomas Brandes
 * @date 21.07.2017
 */
#pragma once

#include <scai/common/config.hpp>
#include <scai/solver/Solver.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{

class DiagonalSolver: public Solver

{
public:
    /**
     * @brief Creates an DiagonalSolver with the specified id.
     *
     * @param[in] id The id of this solver
     */
    DiagonalSolver( const std::string& id ) : Solver( id )
    {
        SCAI_LOG_INFO( Solver::logger, "DiagonalSolver, id = " << id )
    }

    /**
     * @brief Creates an DiagonalSolver with the specified id and logger.
     *
     * @param[in] id     The id of this solver
     * @param[in] logger The logger used by this solver.
     */
    DiagonalSolver( const std::string& id, LoggerPtr logger ) : Solver( id, logger )
    {
        SCAI_LOG_INFO( Solver::logger, "DiagonalSolver, id = " << id )
    }

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    DiagonalSolver( const DiagonalSolver& other ) : Solver( other )
    {
        COMMON_THROWEXCEPTION( "copy not supported for diagonal solver" )
    }

    virtual ~DiagonalSolver()
    {
        SCAI_LOG_INFO( Solver::logger, "~DiagonalSolver" )
    }

    void initialize( const lama::Vector& diagonal )
    {
        getRuntime().mInvDiagonal.reset( diagonal.copy() );
        getRuntime().mInvDiagonal->invert();
        getRuntime().mInitialized = true;
    }

    /**
     * @brief Initializes the solver by inverting and storing the given matrix.
     *
     * @param[in] coefficients  The matrix A from A * u = f.
     */
    virtual void initialize( const lama::Matrix& )
    {
        COMMON_THROWEXCEPTION( "DiagonalSolver cannot be initialized by a matrix." )
    }

    /**
     * @brief Solves the equation system with the given rhs and stores the
     *        result in the given vector.
     *
     * Solves the equation system with the given rhs. Must be initialized first.
     */
    virtual void solveImpl()
    {
        DiagonalSolverRuntime& runtime = getRuntime();
        *runtime.mSolution =  *runtime.mInvDiagonal * *runtime.mRhs;
    }

    struct DiagonalSolverRuntime: SolverRuntime
    {
        DiagonalSolverRuntime()
        {
        }

        virtual ~DiagonalSolverRuntime()
        {
        }
    
        std::shared_ptr<lama::Vector> mInvDiagonal;
    };
    
    virtual SolverPtr copy()
    {
        COMMON_THROWEXCEPTION( "copy unsupported here" )
        return SolverPtr();
    }

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual DiagonalSolverRuntime& getRuntime()
    {
        return mDiagonalSolverRuntime;
    }

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const DiagonalSolverRuntime& getConstRuntime() const
    {
        return mDiagonalSolverRuntime;
    }

protected:

    DiagonalSolverRuntime mDiagonalSolverRuntime;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "DiagonalSolver( id = " << mId << " )";
    }

};

} /* end namespace solver */

} /* end namespace scai */
