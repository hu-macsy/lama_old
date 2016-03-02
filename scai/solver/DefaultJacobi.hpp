/**
 * @file DefaultJacobi.hpp
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
 * @brief DefaultJacobi.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/OmegaSolver.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{

class COMMON_DLL_IMPORTEXPORT DefaultJacobi:
		public OmegaSolver,
		public Solver::Register<DefaultJacobi>
{
public:
    DefaultJacobi( const std::string& id );

    DefaultJacobi( const std::string& id, LoggerPtr logger );

    DefaultJacobi( const std::string& id, const lama::Scalar omega ); //2nd param Matrix.Scalar

    DefaultJacobi( const std::string& id, const lama::Scalar omega, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    DefaultJacobi( const DefaultJacobi& other );

    virtual ~DefaultJacobi();

    /**
     * @brief Initializes the solver by calculating D^(-1)*C from A = D + C.
     *
     * This specialized initialization for the Jacobi solver does the
     * one-time calculation of D^(-1)*C where D is the diagonal Matrix
     * of A and C is the rest (A = D + C).
     *
     * @param coefficients The matrix A from A*u=f
     */
    virtual void initialize( const lama::Matrix& coefficients );

    virtual void solve( lama::Vector& solution, const lama::Vector& rhs );

    virtual void solveInit( lama::Vector& solution, const lama::Vector& rhs );

    virtual void solveFinalize();

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct DefaultJacobiRuntime: OmegaSolverRuntime
    {
        DefaultJacobiRuntime();
        virtual ~DefaultJacobiRuntime();

        common::shared_ptr<lama::Matrix> mDiagonalTimesLU;
        common::shared_ptr<lama::Matrix> mDiagonalInverted;
        common::shared_ptr<lama::Vector> mDiagonalTimesRhs;
        common::shared_ptr<lama::Vector> mOldSolution;
        SolutionProxy mProxyOldSolution;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual DefaultJacobiRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const DefaultJacobiRuntime& getConstRuntime() const;

    static std::string createValue();
    static Solver* create( const std::string name );

protected:
    DefaultJacobiRuntime mDefaultJacobiRuntime;

    /**
     * @brief Performs one Jacobi iteration based on Matrix/Vector operations
     *
     * In addition to the Jacobi iteration the term D^(-1)*f from A=u*f and
     * A = D + C is evaluated before the first iteration.
     */
    virtual void iterate();

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    template<typename ValueType>
    void initialize(const lama::Matrix& coefficients);

    template<typename ValueType>
    void iterate();
};

} /* end namespace solver */

} /* end namespace scai */
