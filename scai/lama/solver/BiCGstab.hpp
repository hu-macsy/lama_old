/**
 * @file BiCGstab.hpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief BiCGstab.hpp
 * @author lschubert
 * @date 06.08.2013
 * $Id$
 */

#pragma once

// for dll import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/Solver.hpp>
#include <scai/lama/solver/IterativeSolver.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief The class BiCGstab represents a IterativeSolver which uses the krylov subspace stabilized BiCG method
 *        to solve a system of linear equations iteratively.
 *
 * Remark: 
 * The scalars in the algorithm are set to zero if they are smaller than machine precision
 * (3*eps) to avoid devision by zero. In this case the solution doesn't change anymore.
 */
class COMMON_DLL_IMPORTEXPORT BiCGstab:
        public IterativeSolver,
        public Solver::Register<BiCGstab>
{
public:
    /**
     * @brief Creates a BiCG solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    BiCGstab( const std::string& id );

    /**
     * @brief Create a BiCG solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    BiCGstab( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    BiCGstab( const BiCGstab& other );

    virtual ~BiCGstab();

    virtual void initialize( const Matrix& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct BiCGstabRuntime: IterativeSolverRuntime
    {
        BiCGstabRuntime();
        virtual ~BiCGstabRuntime();

        common::shared_ptr<Vector> mRes0;
        common::shared_ptr<Vector> mVecV;
        common::shared_ptr<Vector> mVecP;
        common::shared_ptr<Vector> mVecS;
        common::shared_ptr<Vector> mVecT;

        Scalar mEps;
        Scalar mResNorm;
        Scalar mOmega;
        Scalar mAlpha;
        Scalar mBeta;
        Scalar mRhoOld;
        Scalar mRhoNew;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual BiCGstabRuntime& getRuntime();
    /** 
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( Vector& solution, const Vector& rhs );
    
    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const BiCGstabRuntime& getConstRuntime() const;

    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    virtual void iterate();

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    BiCGstabRuntime    mBiCGstabRuntime;
};

} /* end namespace lama */

} /* end namespace scai */
