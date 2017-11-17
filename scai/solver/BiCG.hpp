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
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT BiCG:

    public CG<ValueType>,
    public _Solver::Register<BiCG<ValueType> >

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

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual BiCG<ValueType>* copy();

    struct BiCGRuntime: CG<ValueType>::CGRuntime
    {
        lama::MatrixPtr<ValueType> mConjTransposeA;

        lama::DenseVector<ValueType> mP2;
        lama::DenseVector<ValueType> mQ2;
        lama::DenseVector<ValueType> mZ2;
        ValueType mPScalar2;
        mutable lama::DenseVector<ValueType> mResidual2;
    };

    const lama::Vector<ValueType>& getResidual2() const;

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual BiCGRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const BiCGRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

protected:

    using IterativeSolver<ValueType>::mPreconditioner;

    virtual void iterate();

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    BiCGRuntime    mBiCGRuntime;
};

} /* end namespace solver */

} /* end namespace scai */
