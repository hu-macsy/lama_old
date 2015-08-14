/**
 * @file CG.hpp
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
 * @brief CG.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */
#ifndef LAMA_CG_HPP_
#define LAMA_CG_HPP_

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/IterativeSolver.hpp>
namespace scai
{

namespace lama
{

/**
 * @brief The class CG represents a IterativeSolver which uses the krylov subspace CG method
 *        to solve a system of linear equations iteratively.
 */
class COMMON_DLL_IMPORTEXPORT CG: public IterativeSolver
{
public:
    /**
     * @brief Creates a CG solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    CG( const std::string& id );

    /**
     * @brief Create a CG solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    CG( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    CG( const CG& other );

    virtual ~CG();

    virtual void initialize( const Matrix& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct CGRuntime: IterativeSolverRuntime
    {
        CGRuntime();
        virtual ~CGRuntime();

        common::shared_ptr<Vector> mP;
        common::shared_ptr<Vector> mQ;
        common::shared_ptr<Vector> mZ;
        Scalar mPScalar;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual CGRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const CGRuntime& getConstRuntime() const;

    double getAverageIterationTime() const;
    double getAveragePreconditionerTime() const;

protected:

    virtual void iterate();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    CGRuntime    mCGRuntime;

private:

    double totalIterationTime;
    double totalPreconditionerTime;
};

} /* end namespace lama */

} /* end namespace scai */

#endif // LAMA_CG_HPP_
