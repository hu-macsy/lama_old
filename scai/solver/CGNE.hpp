/**
 * @file CGNE.hpp
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
 * @brief CGNE.hpp
 * @author David Schissler
 * @date 30.10.2015
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
 * @brief The class CGNE represents an IterativeSolver which uses the krylov subspace Conjugate
 * Gradients for Normal Equations (CGNE) method to solve a system of linear equations of type  
 * 
 *                             (A * A^t) * x = b
 * iteratively 
 * where 
 * A  .............. is some matrix (not necessary square) 
 * A^t.............. is the transposed of A
 * x  .............. solution vector
 * b  .............. rhs vector
 * Remark:
 * The scalars in the algorithm are set to zero if they are smaller then machine
 * precision to avoid division by zero. In this case the solution doesn't change anymore.
 *
 */
class COMMON_DLL_IMPORTEXPORT CGNE:
		public IterativeSolver,
		public Solver::Register<CGNE>
{
public:
    /**
    * @brief Creates a CGNE solver with a given ID.
    *
    * @param id The ID for the solver.
    */
    CGNE( const std::string& id );
    /**
    * @brief Create a CGNE solver with a given ID and a given logger.
    *
    * @param id        The ID of the solver.
    * @param logger    The logger which shall be used by the solver
    */
    CGNE( const std::string& id, LoggerPtr logger );

    /**
    * @brief Copy constructor that copies the status independent solver information
    */
    CGNE( const CGNE& other );

    virtual ~CGNE();

    virtual void initialize( const lama::Matrix& coefficients );

    /**
    * @brief Copies the status independent solver informations to create a new instance of the same
    * type
    *
    * @return shared pointer of the copied solver
    */
    virtual SolverPtr copy();

    struct CGNERuntime: IterativeSolverRuntime
    {
        CGNERuntime();
        virtual ~CGNERuntime();

    common::shared_ptr<lama::Matrix> mTransposedMat;
	common::shared_ptr<lama::Vector> mVecP;
    common::shared_ptr<lama::Vector> mVecZ;

    lama::Scalar mEps;
    };
    /**
    * @brief Returns the complete configuration of the derived class
    */
    virtual CGNERuntime& getRuntime();
    /** 
    * @brief Initializes vectors and values of the runtime
    */
    virtual void solveInit( lama::Vector& solution, const lama::Vector& rhs );

    /**
    * @brief Returns the complete const configuration of the derived class
    */
    virtual const CGNERuntime& getConstRuntime() const;
    
    static std::string createValue();
    static Solver* create( const std::string name );

protected:

    CGNERuntime mCGNERuntime;
    /**
     * @brief Performs one CGNE iteration based on Matrix/Vector operations. 
     */
    virtual void iterate();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace solver */

} /* end namespace scai */
