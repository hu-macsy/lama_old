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

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/IterativeSolver.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief The class CG represents a IterativeSolver which uses the krylov subspace CG method
 *        to solve a system of linear equations iteratively.
 *
 * Remark: 
 * The scalars in the algorithm are set to zero if they are smaller than machine precision
 * (3*eps) to avoid devision by zero. In this case the solution doesn't change anymore.
 */
class COMMON_DLL_IMPORTEXPORT CG:
	public IterativeSolver,
	public Solver::Register<CG>
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

    virtual void initialize( const lama::Matrix& coefficients );

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

        lama::VectorPtr mP;
        lama::VectorPtr mQ;
        lama::VectorPtr mZ;
        lama::Scalar mPScalar;
        lama::Scalar mEps;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual CGRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const CGRuntime& getConstRuntime() const;

    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    virtual void iterate();

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    CGRuntime    mCGRuntime;
};

} /* end namespace solver */

} /* end namespace scai */
