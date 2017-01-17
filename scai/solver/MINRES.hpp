/**
 * @file MINRES.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief MINRES.hpp
 * @author David Schissler
 * @date 29.05.2015
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
