/**
 * @file CGS.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief CGS.hpp
 * @author David Schissler
 * @date 18.05.2015
 */

#pragma once

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
 * @brief The class CGS represents a IterativeSolver which uses the krylov subspace Conjugate Gradient Squared
 * method to solve a system of linear equations iteratively. Keep in mind that this method is not stable.
 *
 * Remarks:
 * 1. This method is numerically unstable. This effect gets a bit annulled by (2.).
 * 2. The scalars in the algorithm are set to zero if they are smaller than machine precision
 * (3*eps) to avoid devision by zero. In this case the solution doesn't change anymore.
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

    virtual void initialize( const lama::Matrix& coefficients );

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

        common::shared_ptr<lama::Vector> mRes0;
        common::shared_ptr<lama::Vector> mVecP;
        common::shared_ptr<lama::Vector> mVecQ;
        common::shared_ptr<lama::Vector> mVecU;
        common::shared_ptr<lama::Vector> mVecT;
        common::shared_ptr<lama::Vector> mVecPT;
        common::shared_ptr<lama::Vector> mVecUT;
        common::shared_ptr<lama::Vector> mVecTemp;

        lama::Scalar mEps;
        lama::Scalar mNormRes;
        lama::Scalar mInnerProdRes;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual CGSRuntime& getRuntime();
    /** 
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( lama::Vector& solution, const lama::Vector& rhs );
    
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

} /* end namespace solver */

} /* end namespace scai */
