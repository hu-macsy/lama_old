/**
 * @file Solver.hpp
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
 * @brief Solver superclass. Direct solvers can be derived from here,
 * iterative sovlers should use the IterativeSolver class.
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */
#ifndef LAMA_SOLVER_HPP_
#define LAMA_SOLVER_HPP_

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>
#include <scai/common/Printable.hpp>

// others
#include <scai/lama/Vector.hpp>

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/solver/SolutionProxy.hpp>
#include <scai/lama/solver/logger/Logger.hpp>

// logging
#include <scai/logging.hpp>

#include <string>
#include <memory>

namespace lama
{

class Solver;
typedef common::shared_ptr<Solver> SolverPtr;

/**
 * @brief Superclass for all solvers.
 *
 * This class acts as a superclass for all solver. It offers basic
 * functionality for coefficient, rhs and solution storing, provides
 * a custom ID for a solver and a residual calculation capabilities.
 */
class COMMON_DLL_IMPORTEXPORT Solver: public Printable
{
public:
    /**
     * @brief Creates a solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    Solver( const std::string& id );

    /**
     * @brief Create a solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    Solver( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    Solver( const Solver& other );

    /**
     * @brief Solver destructor.
     */
    virtual ~Solver();

    /**
     * @brief Used to initialize a solver with a certain matrix A
     *        from A*u=f.
     *
     * This method initializes a solver with a certain coefficient-matrix.
     * The only thing it does is storing the matrix pointer as a member
     * for derived solver classes to use it. The caller delegates the
     * property of the pointer to the Solver instance.
     *
     * This method may be overwritten by base classes which desire
     * more complex initialization.
     *
     * @param coefficients The matrix A from A*u=f.
     */
    virtual void initialize( const Matrix& coefficients );

    /**
     * @brief Solves the equation system based on the given rhs.
     *
     * The solver needs to be initialized first with the matrix from
     * the equation to solve, e.g.
     *     A from A*u=f (call solver::initialize(A) for example)
     * This method is abstract. It has to be implemented by a class which inherits
     * from this class.
     *
     * @param rhs       The right hand side of A*u=f.
     * @param solution  The solution from A*u=f. Mostly used as starting
     *                  solution for an IterativeSolver.
     */
    virtual void solve( Vector& solution, const Vector& rhs );

    /**
     * @brief Initializes the solver with rhs and solution.
     *
     * @param[in]  rhs      The right hand side of the system of equations
     * @param[out] solution The allocated memory and starting solution for the system
     */
    virtual void solveInit( Vector& solution, const Vector& rhs );

    /**
     * @brief Solves the equation system. Rhs and starting solution have to
     * be initialized first! (call solveInit( rhs, solution ) ).
     * The solver needs to be initialized first with the matrix from
     * the equation to solve, e.g.
     *     A from A*u=f (call solver::initialize(A) for example)
     *
     * This method solves the equation system by using the given rhs and
     * solution. For most solvers the solution-vector is used as
     * a starting solution for the solve process. This class does not take
     * responsibility for deleting the vectors after the solver! Make sure
     * you do not delete the vectors during the solver process.
     */
    virtual void solveImpl() = 0;

    /**
     * @brief Finalizes the solving process.
     */
    virtual void solveFinalize();

    /**
     * @brief Returns the ID of this solver
     *
     * @return The ID of this solver.
     */
    const std::string& getId() const;

    /**
     * @brief Contingently calculates the current residual based on the
     *        coefficients, rhs and solution currently associated with this.
     *
     * Should be used internally only, because the three vectors
     * mentioned above have to be initialized.
     *
     * @return The current residual
     */
    const Vector& getResidual() const;

    /**
     * @brief Gets the matrix A from A*u=f.
     *
     * @return The coefficient matrix A.
     */
    const Matrix& getCoefficients() const;

    /**
     * @brief Redefines mLogger
     *
     */
    void setLogger( LoggerPtr logger );

    /**
     * @brief Switches the loglevel of mLogger
     *
     * @return
     */
    void setLogLevel( LogLevel::LogLevel level );

    /**
     * @brief Sets the context where this solver should be executed.
     *
     * Sets the context where this solver should be executed. Caution: This overrides
     * the context of the coefficients matrix A from A * u = f used to inializ this
     * solver.
     *
     * @param[in] context   the context where this solver should be executed.
     */
    virtual void setContext( ContextPtr context );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy() =0;

    /**
     * @brief Status independent solver informations
     */
    struct SolverRuntime: public common::NonCopyable
    {
        SolverRuntime();
        virtual ~SolverRuntime();

        /**
         * @brief The coefficient matrix A.
         */
        const Matrix* mCoefficients;

        /**
         * @brief The right-hand-side f.
         */
        const Vector* mRhs;

        /**
         * @brief The solution u (using the SolutionProxy).
         */
        mutable SolutionProxy mSolution;

        /**
         * @brief The residual.
         */
        mutable common::shared_ptr<Vector> mResidual;

        /**
         * @brief Flag for initialization status of solver.
         */
        bool mInitialized;

        /**
         * @brief Flag for initialization status of solve-routine.
         */
        bool mSolveInit;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual SolverRuntime& getRuntime() =0;

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const SolverRuntime& getConstRuntime() const =0;

protected:

    /**
     * @brief The ID of this solver.
     */
    std::string mId;

    /**
     * @brief The solver logger.
     *
     * May be the NullLogger if no logger has been specified.
     */
    LoggerPtr mLogger;

    /**
     * @brief For forcing the context solver dependent
     *
     * If the context for a solver is set, the context of the input matrix will be ignored
     */
    ContextPtr mContext;

    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

}
// namespace lama

#endif // LAMA_SOLVER_HPP_
