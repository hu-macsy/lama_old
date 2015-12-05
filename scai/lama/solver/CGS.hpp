/**
 * @file CGS.hpp
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
 * @brief CGS.hpp
 * @author David Schissler
 * @date 18.05.2015
 * @since 
 */

#pragma once

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
 * @brief The class CGS represents a IterativeSolver which uses the krylov subspace Conjugate Gradient Squared
 * method to solve a system of linear equations iteratively. Keep in mind that this method is not stable.
 *
 * Remarks:
 * 1. This method is not numerically stable. This effect gets a bit annulled by 2.
 * 2. The scalars in the algorithm are set to zero if the norm of the residual is smaller than 
 * machine precision (3*eps) to avoid devision by zero. In this case the solution doesn't
 * change anymore.
 * 3. In this case it makes sense to take the residual since we have to update the residual in each
 * iterate() anyways (contrary to e.g. TFQMR solver).
 */
class COMMON_DLL_IMPORTEXPORT CGS:
		public IterativeSolver,
		public Solver::Register<CGS>
{
public:
    /**
     * @brief Creates a CGS solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    CGS( const std::string& id );

    /**
     * @brief Create a CGS solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    CGS( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    CGS( const CGS& other );

    virtual ~CGS();

    virtual void initialize( const Matrix& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct CGSRuntime: IterativeSolverRuntime
    {
        CGSRuntime();
        virtual ~CGSRuntime();

        common::shared_ptr<Vector> mRes0;
        common::shared_ptr<Vector> mVecP;
        common::shared_ptr<Vector> mVecQ;
        common::shared_ptr<Vector> mVecU;
        common::shared_ptr<Vector> mVecT;

        Scalar mEps;
        Scalar mNormRes;
        Scalar mInnerProdRes;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual CGSRuntime& getRuntime();
    /** 
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( Vector& solution, const Vector& rhs );
    
    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const CGSRuntime& getConstRuntime() const;

    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    virtual void iterate();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    CGSRuntime    mCGSRuntime;
};

} /* end namespace lama */

} /* end namespace scai */
