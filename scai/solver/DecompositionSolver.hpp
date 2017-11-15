/**
 * @file DecompositionSolver.hpp
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
 * @brief Contains the class DecompositionSolver.
 * @author Lauretta Schubert
 * @date 20.07.2016
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
class COMMON_DLL_IMPORTEXPORT DecompositionSolver:

    public Solver<ValueType>,
    public _Solver::Register<DecompositionSolver<ValueType> >
{
public:
    /**
     * @brief Creates an DecompositionSolver with the specified id.
     *
     * @param[in] id The id of this solver
     */
    DecompositionSolver( const std::string& id );

    /**
     * @brief Creates an DecompositionSolver with the specified id and logger.
     *
     * @param[in] id     The id of this solver
     * @param[in] logger The logger used by this solver.
     */
    DecompositionSolver( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    DecompositionSolver( const DecompositionSolver& other );

    virtual ~DecompositionSolver();

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

    void solveImplTyped( const lama::SparseMatrix<ValueType>& coefficients );

    struct DecompositionSolverRuntime: Solver<ValueType>::SolverRuntime
    {
        bool mIsSymmetric;
    };

    virtual DecompositionSolver<ValueType>* copy();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual DecompositionSolverRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const DecompositionSolverRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

protected:

    DecompositionSolverRuntime mDecompositionSolverRuntime;

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
