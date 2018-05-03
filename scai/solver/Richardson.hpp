/**
 * @file solver/Richardson.hpp
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
 * @brief Richardson.hpp
 * @author David Schissler
 * @date 17.04.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/Solver.hpp>
#include <scai/solver/OmegaSolver.hpp>

// logging
#include <scai/logging/Logger.hpp>

namespace scai
{

namespace solver
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Richardson:

    public OmegaSolver<ValueType>,
    public _Solver::Register<Richardson<ValueType> >

{
public:

    Richardson( const std::string& id );

    Richardson( const std::string& id, LoggerPtr logger );

    Richardson( const std::string& id, const ValueType omega );

    Richardson( const std::string& id, const ValueType omega, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    Richardson( const Richardson& other );

    virtual ~Richardson();

    /**
     * @param coefficients  The matrix A from A*u=f
     *
     *This method converges if      omega < 2 / ||A||_2     holds.
     *To assure this we choose omega s.t. the lower bound holds
     *        omega < 2 / ||A||_F <= 2 / ||A||_2
     *with
     *||.||_2 ...... spectral norm,
     *||.||_F ...... Frobenius norm.
     *
     *In particular, we take omega = (2/3) * 2 / ||A||_F  in initialize()
     *if there is no initialized omega.
    */

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    virtual void solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual Richardson<ValueType>* copy();

    struct RichardsonRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        std::unique_ptr<lama::Vector<ValueType>> mOldSolution;
        std::unique_ptr<lama::Vector<ValueType>> mX;
    };
    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual RichardsonRuntime& getRuntime();

    /**
     * @brief Returns the complete const configuration of the derived class
     */
    virtual const RichardsonRuntime& getRuntime() const;

    // static method that delivers the key for registration in solver factor

    static SolverCreateKeyType createValue();

    // static method for create by factory

    static _Solver* create();

protected:

    RichardsonRuntime mRichardsonRuntime;

    /**
     * @brief Performs one Richardson iteration based on _Matrix/Vector operations
     */
    virtual void iterate();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /**
     *  @brief own implementation of Printable::writeAt
     */
    virtual void writeAt( std::ostream& stream ) const;

};

} /* end namespace solver */

} /* end namespace scai */
