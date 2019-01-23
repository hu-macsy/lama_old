/**
 * @file CG.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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

namespace lama
{
    template<typename ValueType> class Matrix;
    template<typename ValueType> class Vector;
}

namespace solver
{

/**
 * @brief The class CG represents a IterativeSolver which uses the krylov subspace CG method
 *        to solve a system of linear equations iteratively.
 */
template<class ValueType>
class COMMON_DLL_IMPORTEXPORT CG:

    public IterativeSolver<ValueType>,
    public _Solver::Register<CG<ValueType> >    // register at solver factory

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

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same type
     */
    virtual CG<ValueType>* copy();

    using typename IterativeSolver<ValueType>::IterativeSolverRuntime;

    struct CGRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        /** Initialize the runtime with, allocates the temporary vectors. */

        void initialize();

        std::unique_ptr<lama::Vector<ValueType>> mP;
        std::unique_ptr<lama::Vector<ValueType>> mQ;
        std::unique_ptr<lama::Vector<ValueType>> mZ;

        ValueType mPScalar;
  
        using Solver<ValueType>::SolverRuntime::mCoefficients;
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual CGRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const CGRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

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
