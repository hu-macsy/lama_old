/**
 * @file InverseSolver.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Contains the class InversteSolver.
 * @author Matthias Makulla
 * @date 06.04.2011
 * $Id$
 */
#ifndef LAMA_INVERSESOLVER_HPP_
#define LAMA_INVERSESOLVER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/solver/Solver.hpp>

// others
#include <lama/matrix/DenseMatrix.hpp>

// logging
#include <logging/Logger.hpp>

namespace lama
{

/**
 * @brief Uses matrix inverse to solve the equation system.
 */
class LAMA_DLL_IMPORTEXPORT InverseSolver: public Solver
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
    virtual void initialize( const Matrix& coefficients );

    /**
     * @brief Solves the equation system with the given rhs and stores the
     *        result in the given vector.
     *
     * Solves the equation system with the given rhs. Must be initialized first.
     */
    virtual void solveImpl();

    virtual void setContext( ContextPtr context );

    void computeInverse( Matrix& matrix ) const;

    struct InverseSolverRuntime: SolverRuntime
    {
        InverseSolverRuntime();
        virtual ~InverseSolverRuntime();

        boost::shared_ptr<Matrix> mInverse;
    };

    virtual SolverPtr copy();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual InverseSolverRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const InverseSolverRuntime& getConstRuntime() const;

protected:
    InverseSolverRuntime mInverseSolverRuntime;

private:

    template<typename T>
    void invert( DenseMatrix<T>& matrix, IndexType* const permutation ) const;

    template<typename T>
    void decompose( DenseMatrix<T>& matrix, IndexType* const permutation ) const;

    void logStartSolve();
    void logEndSolve();

    LAMA_LOG_DECL_STATIC_LOGGER(logger);
};

}

#endif // LAMA_INVERSESOLVER_HPP_
