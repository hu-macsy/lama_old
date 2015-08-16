/**
 * @file CGNR.hpp
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
 * @brief CGNR.hpp
 * @author David Schissler
 * @date 27.05.2015
 * @since 
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/solver/IterativeSolver.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace lama
{
/**
 * @brief The class CGNR represents a IterativeSolver which uses the krylov subspace Conjugate
 * Gradients for Normal Residuals (CGNR) method to solve a system of linear equations iteratively.
 * Remarks:
 * 1. assignTranspose is not supported yet for DenseMatrix.
 * 2. The scalars in the algorithm are set to zero if they are smaller then machine
 * precision (3*eps) to avoid devision by zero. In this case the solution doesn't change anymore.
 * 3. In this case it makes less sense to take the residual regarding to some norm itself since 
 * it has to be additionally computed in each iterate() (contrary to e.g. BiCGstab solver;
 * no higher costs). 
 */
class COMMON_DLL_IMPORTEXPORT CGNR: public IterativeSolver
{
public:
    /**
    * @brief Creates a CGNR solver with a given ID.
    *
    * @param id The ID for the solver.
    */
    CGNR( const std::string& id );
    /**
    * @brief Create a CGNR solver with a given ID and a given logger.
    *
    * @param id        The ID of the solver.
    * @param logger    The logger which shall be used by the solver
    */
    CGNR( const std::string& id, LoggerPtr logger );

    /**
    * @brief Copy constructor that copies the status independent solver information
    */
    CGNR( const CGNR& other );

    virtual ~CGNR();

    virtual void initialize( const Matrix& coefficients );

    /**
    * @brief Copies the status independent solver informations to create a new instance of the same
    * type
    *
    * @return shared pointer of the copied solver
    */
    virtual SolverPtr copy();

    struct CGNRRuntime: IterativeSolverRuntime
    {
        CGNRRuntime();
        virtual ~CGNRRuntime();

    common::shared_ptr<Matrix> mTransposedMat;
	common::shared_ptr<Vector> mVecD;
    common::shared_ptr<Vector> mVecW;
    common::shared_ptr<Vector> mVecZ;

    Scalar mEps;
    Scalar mNormVecZ;
    };
    /**
    * @brief Returns the complete configuration of the derived class
    */
    virtual CGNRRuntime& getRuntime();
    /** 
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( Vector& solution, const Vector& rhs );

    /**
    * @brief Returns the complete const configuration of the derived class
    */
    virtual const CGNRRuntime& getConstRuntime() const;
    

protected:

    CGNRRuntime mCGNRRuntime;
    /**
     * @brief Performs one CGNR iteration based on Matrix/Vector operations. 
     */
    virtual void iterate();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */
