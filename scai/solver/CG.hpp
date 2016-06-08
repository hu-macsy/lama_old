/**
 * @file CG.hpp
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief CG.hpp
 * @author Matthias Makulla
 * @date 06.04.2011
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
