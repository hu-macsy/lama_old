/**
 * @file OmegaSolver.hpp
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
 * @endlicense
 *
 * @brief OmegaSolver.hpp
 * @author Kai Buschulte
 * @date 10.08.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/IterativeSolver.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace solver
{

class OmegaSolver;
typedef common::shared_ptr<OmegaSolver> OldSolutionHandlerPtr;

/**
 * @brief The OldSolutionHandler class only manages the omega parameter
 * For solvers like a Jacobi.
 */
class COMMON_DLL_IMPORTEXPORT OmegaSolver: public IterativeSolver
{
public:

    /**
     * @brief Creates a solver with the given id.
     *
     * @param[in] id The id of the solver.
     */
    OmegaSolver( const std::string& id );

    /**
     * @brief Creates a solver with the given id and omega.
     *
     * @param[in] id    The id of the solver.
     * @param[in] omega The omega parameter which is used by the jacobi solver.
     */
    OmegaSolver( const std::string& id, const lama::Scalar omega );

    /**
     * @brief Creates a solver with the given id and logger.
     *
     * @param[in] id     The id of the solver.
     * @param[in] logger The logger which is used by this solver.
     */
    OmegaSolver( const std::string& id, LoggerPtr logger );

    /**
     * @brief Creates a solver with the given id, omega and logger.
     *
     * @param[in] id     The id of the solver.
     * @param[in] omega  The omega parameter which is used by the jacobi solver.
     * @param[in] logger The logger used by the solver.
     */
    OmegaSolver( const std::string& id, const lama::Scalar omega, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    OmegaSolver( const OmegaSolver& other );

    /**
     * @brief Destructor.
     */
    virtual ~OmegaSolver();

    /**
     * @brief This abstract method is used by derived solvers to initialize a
     *        omega solver.
     */
    virtual void initialize( const lama::Matrix& coefficients );

    /**
     * @brief Sets the omega parameter of this.
     *
     * @param[in] omega The omega parameter of the omega solver.
     */
    void setOmega( const lama::Scalar omega );

    /**
     * @brief Returns omega.
     *
     * @return Omega.
     */
    lama::Scalar getOmega() const;

    struct OmegaSolverRuntime: IterativeSolverRuntime
    {
        OmegaSolverRuntime();
        virtual ~OmegaSolverRuntime();
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual OmegaSolverRuntime& getRuntime() =0;

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const OmegaSolverRuntime& getConstRuntime() const =0;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy() =0;

protected:
    /**
     * @brief This abstract method is used by derived solvers to represent a
     *        solver iteration.
     */
    virtual void iterate() = 0;

    lama::Scalar mOmega;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace solver */

} /* end namespace scai */
