/**
 * @file BiCG.hpp
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
 * @brief BiCG.hpp
 * @author Lauretta Schubert
 * @date 03.07.2013
 */

#pragma once

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/CG.hpp>

// internal scai libraries
#include <scai/lama/matrix/SparseMatrix.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief The class BiCG represents a IterativeSolver which uses the krylov subspace BiCG method
 *        to solve a system of linear equations iteratively.
 */
class COMMON_DLL_IMPORTEXPORT BiCG:

    public CG,
    public Solver::Register<BiCG>

{
public:

    /**
     * @brief Creates a BiCG solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    BiCG( const std::string& id );

    /**
     * @brief Create a BiCG solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    BiCG( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    BiCG( const BiCG& other );

    virtual ~BiCG();

    virtual void initialize( const lama::_Matrix& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual SolverPtr copy();

    struct BiCGRuntime: CGRuntime
    {
        BiCGRuntime();
        virtual ~BiCGRuntime();

        common::shared_ptr<lama::_Matrix> mTransposeA;
        common::shared_ptr<lama::_Vector> mP2;
        common::shared_ptr<lama::_Vector> mQ2;
        common::shared_ptr<lama::_Vector> mZ2;
        lama::Scalar mPScalar2;
        mutable common::shared_ptr<lama::_Vector> mResidual2;
    };

    const lama::_Vector& getResidual2() const;

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual BiCGRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const BiCGRuntime& getConstRuntime() const;

    /**
     * @brief returns value used for registration of this solver
     */
    static std::string createValue();

    /**
     * @brief create a new BiCG solver with the corresponding name
     *
     * This method is used as create routine in the Solver factory.
     */
    static Solver* create( const std::string name );

protected:

    virtual void iterate();

    void print( lama::_Vector& vec, size_t n );

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    BiCGRuntime    mBiCGRuntime;
};

} /* end namespace solver */

} /* end namespace scai */
