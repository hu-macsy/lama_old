/**
 * @file OmegaSolver.hpp
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
 * @brief OmegaSolver.hpp
 * @author Kai Buschulte
 * @date 10.08.2011
 * @since 1.0.0
 */

#ifndef LAMA_OMEGASOLVER_HPP_
#define LAMA_OMEGASOLVER_HPP_

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/IterativeSolver.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

class OmegaSolver;
typedef common::shared_ptr<OmegaSolver> OldSolutionHandlerPtr;

/**
 * @brief The OldSolutionHandler class only manages the omega parameter
 * For solvers like Jacobi and SOR
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
    OmegaSolver( const std::string& id, const Scalar omega );

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
    OmegaSolver( const std::string& id, const Scalar omega, LoggerPtr logger );

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
    virtual void initialize( const Matrix& coefficients );

    /**
     * @brief Sets the omega parameter of this.
     *
     * @param[in] omega The omega parameter of the omega solver.
     */
    void setOmega( const Scalar omega );

    /**
     * @brief Returns omega.
     *
     * @return Omega.
     */
    Scalar getOmega() const;

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

    Scalar mOmega;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */

#endif // LAMA_OMEGASOLVER_HPP_
