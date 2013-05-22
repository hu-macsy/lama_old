/**
 * @file IterativeSolver.hpp
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
 * @brief IterativeSolver.h
 * @author Matthias Makulla
 * @date 06.04.2011
 * $Id$
 */

#ifndef LAMA_ITERATIVESOLVER_HPP_
#define LAMA_ITERATIVESOLVER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/Solver.hpp>

// others
#include <lama/solver/criteria/Criterion.hpp>

// logging
#include <logging/Logger.hpp>

namespace lama
{

/**
 * @brief Uses iterative methods to solve the equation system.
 */
class LAMA_DLL_IMPORTEXPORT IterativeSolver: public Solver
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
    virtual void initialize( const Matrix& coefficients );

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
    virtual SolverPtr copy() =0;

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
    virtual IterativeSolverRuntime& getRuntime() =0;

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const IterativeSolverRuntime& getConstRuntime() const =0;

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
     * Logging methods to maintain code-readability
     */
    void logStartSolve();
    void logEndSolve();
    void logIterationEndAndResidual();
    void logIterationStart();

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace lama

#endif // LAMA_ITERATIVESOLVER_HPP_
