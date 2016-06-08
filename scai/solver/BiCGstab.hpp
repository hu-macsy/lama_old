/**
 * @file BiCGstab.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 * @brief BiCGstab.hpp
 * @author lschubert
 * @date 06.08.2013
 */

#pragma once

// for dll import
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

    virtual void initialize( const lama::Matrix& coefficients );

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

        common::shared_ptr<lama::Vector> mRes0;
        common::shared_ptr<lama::Vector> mVecV;
        common::shared_ptr<lama::Vector> mVecP;
        common::shared_ptr<lama::Vector> mVecS;
        common::shared_ptr<lama::Vector> mVecT;
        common::shared_ptr<lama::Vector> mVecPT;
        common::shared_ptr<lama::Vector> mVecST;
        common::shared_ptr<lama::Vector> mVecTT;

        lama::Scalar mEps;
        lama::Scalar mResNorm;
        lama::Scalar mOmega;
        lama::Scalar mAlpha;
        lama::Scalar mBeta;
        lama::Scalar mRhoOld;
        lama::Scalar mRhoNew;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual BiCGstabRuntime& getRuntime();
    /** 
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( lama::Vector& solution, const lama::Vector& rhs );
    
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

} /* end namespace solver */

} /* end namespace scai */
