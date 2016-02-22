/**
 * @file Richardson.hpp
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
 * @brief Richardson.hpp
 * @author David Schissler
 * @date 17.04.2015
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

class COMMON_DLL_IMPORTEXPORT Richardson:
		public OmegaSolver,
		public Solver::Register<Richardson>
{
public:
    Richardson( const std::string& id );

    Richardson( const std::string& id, LoggerPtr logger );

    Richardson( const std::string& id, const lama::Scalar omega ); //2nd param Matrix.Scalar

    Richardson( const std::string& id, const lama::Scalar omega, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    Richardson( const Richardson& other );

    virtual ~Richardson();

    /**
     * @param coefficients  The matrix A from A*u=f
     *
     *This method converges if      omega < 2 / ||A||_2     holds.
     *To assure this we choose omega s.t. the lower bound holds
     *        omega < 2 / ||A||_F <= 2 / ||A||_2
     *with
     *||.||_2 ...... spectral norm,
     *||.||_F ...... Frobenius norm.
     *
     *In particular, we take omega = (2/3) * 2 / ||A||_F  in initialize()
     *if there is no initialized omega.
    */

    virtual void initialize( const lama::Matrix& coefficients );

    virtual void solveInit( lama::Vector& solution, const lama::Vector& rhs );

    virtual void solveFinalize();

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct RichardsonRuntime: OmegaSolverRuntime
    {
        RichardsonRuntime();
        virtual ~RichardsonRuntime();

        common::unique_ptr<lama::Vector> mOldSolution;
        SolutionProxy mProxyOldSolution;
    };
    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual RichardsonRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const RichardsonRuntime& getConstRuntime() const;

    static std::string createValue();
    static Solver* create( const std::string name );

protected:
    RichardsonRuntime mRichardsonRuntime;

    using OmegaSolver::mOmega;

    /**
     * @brief Performs one Richardson iteration based on Matrix/Vector operations
     *
     *
     */
    virtual void iterate();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

private:

    template<typename T>
    void initialize( const lama::Matrix& coefficients );

    template<typename T>
    void iterate();
};

} /* end namespace solver */

} /* end namespace scai */
