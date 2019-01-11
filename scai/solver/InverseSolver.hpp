/**
 * @file InverseSolver.hpp
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
 * @brief Contains the class InversteSolver.
 * @author Thomas Brandes
 * @date 10.04.2013
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>

// internal scai libraries
#include <scai/lama/matrix/DenseMatrix.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief Solver class that uses matrix inverse to solve an equation system.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT InverseSolver:
    public Solver<ValueType>,
    public _Solver::Register<InverseSolver<ValueType> >
{
public:
    /**
     * @brief Creates an InverseSolver with the specified id.
     *
     * @param[in] id The id of this solver
     */
    InverseSolver( const std::string& id );

    /**
     * @brief Creates an InverseSolver with the specified id and logger.
     *
     * @param[in] id     The id of this solver
     * @param[in] logger The logger used by this solver.
     */
    InverseSolver( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    InverseSolver( const InverseSolver& other );

    virtual ~InverseSolver();

    /**
     * @brief Initializes the solver by inverting and storing the given matrix.
     *
     * @param[in] coefficients  The matrix A from A*u=f.
     */
    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
     * @brief Solves the equation system with the given rhs and stores the
     *        result in the given vector.
     *
     * Solves the equation system with the given rhs. Must be initialized first.
     */
    virtual void solveImpl();

    /** This method returns the inverse of the coefficient matrix.
     *
     *  This routine must not be called before having called 'initialize'.
     */
    const lama::Matrix<ValueType>& getInverse() const;

    /** Runtime data for inverse solver, will contain the inverse explicitly. */

    struct InverseSolverRuntime: Solver<ValueType>::SolverRuntime
    {
        lama::MatrixPtr<ValueType> mInverse;
    };

    virtual InverseSolver<ValueType>* copy();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual InverseSolverRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const InverseSolverRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

protected:

    InverseSolverRuntime mInverseSolverRuntime;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    using Solver<ValueType>::mLogger;

private:

    void logStartSolve();
    void logEndSolve();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace solver */

} /* end namespace scai */
