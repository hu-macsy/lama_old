/**
 * @file MINRES.hpp
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
 * @brief MINRES.hpp
 * @author David Schissler
 * @date 29.05.2015
 * @since 
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/IterativeSolver.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{
/**
 * @brief The class MINRES represents a IterativeSolver which uses the krylov subspace Minimum Residual (MINRES)
 * method to solve a system of linear equations iteratively.
 */
class COMMON_DLL_IMPORTEXPORT MINRES:
	public IterativeSolver,
	public Solver::Register<MINRES>
{
public:
    /**
    * @brief Creates a MINRES solver with a given ID.
    *
    * @param id The ID for the solver.
    */
    MINRES( const std::string& id );
    /**
    * @brief Create a MINRES solver with a given ID and a given logger.
    *
    * @param id        The ID of the solver.
    * @param logger    The logger which shall be used by the solver
    */
    MINRES( const std::string& id, LoggerPtr logger );

    /**
    * @brief Copy constructor that copies the status independent solver information
    */
    MINRES( const MINRES& other );

    virtual ~MINRES();

    virtual void initialize( const lama::Matrix& coefficients );

    /**
    * @brief Copies the status independent solver informations to create a new instance of the same
    * type
    *
    * @return shared pointer of the copied solver
    */
    virtual SolverPtr copy();

    struct MINRESRuntime: IterativeSolverRuntime
    {
        MINRESRuntime();
        virtual ~MINRESRuntime();        

    common::shared_ptr<lama::Vector> mVecV;
	common::shared_ptr<lama::Vector> mVecVOld;
    common::shared_ptr<lama::Vector> mVecVNew;
    common::shared_ptr<lama::Vector> mVecP;
    common::shared_ptr<lama::Vector> mVecPOld;
    common::shared_ptr<lama::Vector> mVecPNew;

    lama::Scalar mAlpha;
    lama::Scalar mBetaNew;
    lama::Scalar mBeta;
    lama::Scalar mC;
    lama::Scalar mCOld;
    lama::Scalar mCNew;
    lama::Scalar mS;
    lama::Scalar mSOld;
    lama::Scalar mSNew;
    lama::Scalar mZeta;

    lama::Scalar mEps;
    };
    /**
    * @brief Returns the complete configuration of the derived class
    */
    virtual MINRESRuntime& getRuntime();
    /** 
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( lama::Vector& solution, const lama::Vector& rhs );

    /**
    * @brief Returns the complete const configuration of the derived class
    */
    virtual const MINRESRuntime& getConstRuntime() const;
    
    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    MINRESRuntime mMINRESRuntime;
    /**
     * @brief Performs one MINRES iteration based on Matrix/Vector operations
     */
    virtual void iterate();
    void Lanczos();
    void applyGivensRotation();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;
};

} /* end namespace solver */

} /* end namespace scai */
