/**
 * @file Kaczmarz.hpp
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
 * @brief Implementation of solver using the Kaczmarz method
 * @author Thomas Brandes
 * @date 18.11.2016
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
 * @brief The class Kaczmarz represents a IterativeSolver which uses the Kaczmarz method
 *        to solve a system of linear equations iteratively.
 *
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Kaczmarz:

    public IterativeSolver<ValueType>,
    public _Solver::Register<Kaczmarz<ValueType> >

{
public:
    /**
     * @brief Creates a Kaczmarz solver with a given ID.
     *
     * @param id The ID for the solver.
     */
    Kaczmarz( const std::string& id );

    /**
     * @brief Create a Kaczmarz solver with a given ID and a given logger.
     *
     * @param id        The ID of the solver.
     * @param logger    The logger which shall be used by the solver
     */
    Kaczmarz( const std::string& id, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    Kaczmarz( const Kaczmarz& other );

    virtual ~Kaczmarz();

    virtual void initialize( const lama::Matrix<ValueType>& coefficients );

    /**
     * @brief Copies the status independent solver informations to create a new instance of the same
     * type
     *
     * @return shared pointer of the copied solver
     */
    virtual Kaczmarz<ValueType>* copy();

    struct KaczmarzRuntime: IterativeSolver<ValueType>::IterativeSolverRuntime
    {
        std::unique_ptr<lama::Vector<ValueType>> mRow;   // temporary for row, might be sparse or dense
    };

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual KaczmarzRuntime& getRuntime();

    /**
     * @brief Returns the complete configuration of the derived class
     */
    virtual const KaczmarzRuntime& getRuntime() const;

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

    KaczmarzRuntime mKaczmarzRuntime;
};

} /* end namespace solver */

} /* end namespace scai */
