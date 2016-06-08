/**
 * @file IterativeSolver.hpp
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
 * @brief IterativeSolver.h
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>

// local library
#include <scai/solver/criteria/Criterion.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{

class IterativeSolver;
typedef common::shared_ptr<IterativeSolver> IterativeSolverPtr;

/**
 * @brief Uses iterative methods to solve the equation system.
 */
class COMMON_DLL_IMPORTEXPORT IterativeSolver:
    public Solver
{
public:
    /**
     * @brief Creates a solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    IterativeSolver( const std::string& id );

    /**
     * @brief Constructs a new solver with the given id and logger.
     *
     * @param id The name or id of the solver.
     * @param logger The logger which the solver shall use
     */
    IterativeSolver( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    IterativeSolver( const IterativeSolver& other );

    /**
     * @brief Destructor for the IterativeSolver class.
     */
    virtual ~IterativeSolver();

    /**
     * @brief Used to initialize a solver with a certain matrix A
     *        from A*u=f.
     *
     * This method initializes a solver with a certain coefficient-matrix.
     * The only thing it does is storing the matrix pointer as a member
     * for derived solver classes to use it. The caller delegates the
     * property of the pointer to the Solver instance.
     *
     * This method may be overwritten by derived classes which desire
     * more complex initialization.
     *
     * @param coefficients The matrix A from A*u=f.
     */
    virtual void initialize( const lama::Matrix& coefficients );

    /**
     * @brief Solves the equation system. Rhs and starting solution have to
     * be initialized first! (call solveInit( rhs, solution ) ).
     * The solver needs to be initialized first with the matrix from
     * the equation to solve, e.g.
     *     A from A*u=f (call solver::initialize(A) for example)
     *
     * This method solves the equation system by using the given rhs and
     * solution. For most iterative solvers the solution-vector is used as
     * a starting solution for the solve process. This class does not take
     * responsibility for deleting the vectors after the solver! Make sure
     * you do not delete the vectors during the solver process.
     */
    virtual void solveImpl();

    /**
     * @brief set a new StoppingCriterion to the solver.
     *
     * @param[in] criterion the new criterion.
     */
    void setStoppingCriterion( const CriterionPtr criterion );

    /**
     * @brief Sets the preconditioner of this solver.
     *
     * Preconditioner should be set before initializing the this solver, otherwise the
     * preconditioner initialization should be executed manually
     *
     * @param conditioner The preconditioner
     */
    void setPreconditioner( SolverPtr const conditioner );

    /**
     * @brief returns the preconditioner of this solver
     *
     * @return the preconditioner
     */
    const SolverPtr getPreconditioner() const;

    /**
     * @brief returns the number of iterations, this solver has done so far.
     *
     */
    int getIterationCount() const;

    /**
     * @brief Copies the status indepedent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy() = 0;

    struct IterativeSolverRuntime: SolverRuntime
    {
        IterativeSolverRuntime();
        virtual ~IterativeSolverRuntime();

        /**
         * @brief The number of iterations the solver currently has performed.
         *
         * This is needed by the iteraton count stopping criteria and is
         * maintained in the solver() method.
         */
        IndexType mIterations;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual IterativeSolverRuntime& getRuntime() = 0;

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const IterativeSolverRuntime& getConstRuntime() const = 0;

    static IterativeSolver* create( const std::string type, const std::string name );

protected:

    /**
     * @brief Checks if all of the stopping criteria are satisfied.
     *
     * return Whether the criteria are satisfied or not
     */
    bool criteriaAreSatisfied() const;

    /**
     * @brief Represents one iteration step of the solver. Use this
     *        abstract method to implement new iterative solvers
     *        (e.g. Jacobi)
     *
     * This method represents one solver iteration. Idealically a derived
     * iterative solver only overwrites this method. It is only used by the
     * solve() method of this class. Use the member variables m_rhs, mSolution
     * and m_coefficients during the iteration process.
     */
    virtual void iterate() = 0;

    /**
     * @brief The preconditioner of this solver.
     */
    SolverPtr mPreconditioner;

    /**
     * @brief The root stopping criterion
     * evaluated every iteration in the solve method
     */
    CriterionPtr mCriterionRootComponent;

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    /**
     * Logging methods to maintain code-readability
     */
    void logStartSolve();
    void logEndSolve();
    void logIterationEndAndResidual();
    void logIterationStart();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace solver */

} /* end namespace scai */
