/**
 * @file Solver.hpp
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
 * @brief Solver superclass. Direct solvers can be derived from here,
 *        iterative sovlers should use the IterativeSolver class.
 * @author The LAMA Team
 * @date 06.04.2011
 */

#pragma once

// base classes

#include <scai/solver/_Solver.hpp>
#include <scai/solver/SolutionProxy.hpp>

namespace scai
{

namespace lama
{
    template<typename ValueType> class Matrix;
    template<typename ValueType> class Vector;
}

/** @brief Namespace for all solver classes, used for all stuff of project solver */

namespace solver
{

/* ========================================================================= */
/*     Solver<ValueType>                                                     */
/* ========================================================================= */

/**
 * @brief Base class for all typed solvers.
 *
 * This class acts as a superclass for all solver. It offers basic
 * functionality for coefficient, rhs and solution storing, provides
 * a custom ID for a solver and a residual calculation capabilities.
 */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Solver: public _Solver
{
public:
    /**
     * @brief Creates a solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    Solver( const std::string& id );

    Solver( const std::string& id, LoggerPtr logger );
 
    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    Solver( const Solver& other );

    /**
     * @brief Solver destructor.
     */
    virtual ~Solver();

    /** Implementation of pure method _Solver::getValueType */

    virtual common::ScalarType getValueType() const;

    /**
     * @brief Create a new solver of a certain type.
     *
     */
    static Solver* getSolver( const std::string& solverType );

    /**
     *  Get all solver types available for this value type by using _Solver::createValues 
     */
    static void getCreateValues( std::vector<std::string>& values );

    /**
     * @brief Query if a solver can be generated
     *
     */
    static bool canCreate( const std::string& solverType );
     
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
    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
     * @brief Solves the equation system based on the given rhs.
     *
     * The solver needs to be initialized first with the matrix from
     * the equation to solve, e.g.
     *     A from A*u=f (call solver::initialize(A) for example)
     * This method is abstract. It has to be implemented by a class which inherits
     * from this class.
     *
     * @param[in,out] solution  The solution from A * solution = rhs. Also used as starting
     *                          solution for an IterativeSolver.
     * @param[in,out] rhs       The right hand side 
     */
    void solve( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs );

    /**
     * @brief Initializes the solver with rhs and solution.
     *
     * @param[in]  rhs      The right hand side of the system of equations
     * @param[out] solution The allocated memory and starting solution for the system
     */
    virtual void solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs );

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
     * @brief Contingently calculates the current residual based on the
     *        coefficients, rhs and solution currently associated with this.
     *
     * Should be used internally only, because the three vectors
     * mentioned above have to be initialized.
     *
     * @return The current residual
     */
    const lama::Vector<ValueType>& getResidual() const;

    /**
     * @brief Gets the matrix A from A*u=f.
     *
     * @return The coefficient matrix A.
     * @throws Exception if solve
     */
    const lama::Matrix<ValueType>& getCoefficients() const;

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual Solver<ValueType>* copy() = 0;

    /**
     * @brief Status independent solver informations
     *
     * The runtime will contain all that data that is never copied when
     * a solver is copied.
     */
    struct SolverRuntime: public common::NonCopyable
    {
        SolverRuntime();
        virtual ~SolverRuntime();

        /**
         * @brief The coefficient matrix A, here a pointer reference
         *
         * Note: A solver assumes that the matrix is not changed during the lifetime of the solver
         */
        const lama::Matrix<ValueType>* mCoefficients;

        /**
         * @brief The right-hand-side f.
         */
        const lama::Vector<ValueType>* mRhs;

        /**
         * @brief The solution u (using the SolutionProxy).
         */
        mutable SolutionProxy<ValueType> mSolution;

        /**
         * @brief The residual.
         */
        std::unique_ptr<lama::Vector<ValueType>> mResidual;

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
    virtual SolverRuntime& getRuntime() = 0;

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const SolverRuntime& getRuntime() const = 0;

protected:

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

};

/** 
 * Definiton of corresponding shared pointer type for the class Solver<ValueType> by a type alias.
 *
 *  \code
 *      SolverPtr<ValueType> x( Solver<ValueType>::getSolver( "CG" ) );
 *      std::shared_ptr<Solver<ValueType> > x( Solver<ValueType>::getSolver( "CG" ) );
 *  \endcode
*/
template<typename ValueType>
using SolverPtr = std::shared_ptr<Solver<ValueType> >;

} /* end namespace solver */

} /* end namespace scai */
