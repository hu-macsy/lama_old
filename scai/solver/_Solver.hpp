/**
 * @file _Solver.hpp
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
 * @brief Solver base class for all solvers.
 * @author Thomas Brandes
 * @date 06.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

#include <scai/solver/logger/SolverLogger.hpp>

// logging
#include <scai/logging.hpp>

// std
#include <string>
#include <memory>

namespace scai
{

/** @brief Namespace for all solver classes, used for all stuff of project solver */

namespace solver
{

class _Solver;

typedef std::pair<common::ScalarType, std::string> SolverCreateKeyType;

/* ========================================================================= */
/*     _Solver : Helper class as common base class for all typed solvers     */
/* ========================================================================= */

/** Common base class for all typed Solver classes. 
 *
 *  This class is helpful for dynamic creation of solver objects where the value type
 *  is only available at runtime. It also contains the type-independent member variables
 *  of the Solver class.
 */
class COMMON_DLL_IMPORTEXPORT _Solver : 

    public common::Printable,
    public common::Factory<SolverCreateKeyType, _Solver*>
{
public:

    /**
     *  Provide a more convenient interface to the create method of the factory.
     */
    static _Solver* getSolver( const common::ScalarType scalarType, const std::string& solverType );

    /**
     * Get logger create values only for a certain type 
     *
     * @param[out] values is a vector of strings containing the types of available solvers
     * @param[in]  stype specifies the scalar type for which solver types are queried.
     *
     * * only registered solvers can be queried
     * * some solvers might only be supported for certain value types
     */
    static void getTypedCreateValues( std::vector<std::string>& values, const common::ScalarType stype );

    /**
     * @brief Queries the value type of the solver 
     *
     * This method allows a safe reinterpret_cast of untyped solvers to typed solvers.
     *
     * \code
     *    _Solver& s = ...
     *    if ( s.getValueType() == ScalarType::DOUBLE )
     *    {
     *        Solver<double>& sd = reinterpret_cast<Solver<double>&>( s );
     *        ...
     *    }
     * \endcode
     */
    virtual common::ScalarType getValueType() const = 0;

     /**
      * @brief Create a solver with a given ID and a default logger.
      */
    _Solver( const std::string& id );

     /**
      * @brief Create a solver with a given ID and a given logger.
      *
      * @param id        The ID of the solver.
      * @param logger    The logger which shall be used by the solver
      */
    _Solver( const std::string& id, LoggerPtr logger );
 
    /** 
     *  @brief copy constructor
     *
     *  Be careful: copied solver shares the logger with the original one
     */
    _Solver( const _Solver& other );

    /**
     * @brief Set or redefine mLogger (convergence history, ... )
     */
    void setLogger( LoggerPtr logger );

    /**
     * @brief Switches the loglevel of mLogger
     *
     * Be careful: as solvers can share the logger, this might have side
     * effects for other solvers. On the other side, it might be useful
     * for solver hierarchies as used in multigrid methods.
     */
    void setLogLevel( LogLevel level );

    /**
     * @brief Get the id of this solver
     */
    inline const std::string& getId() const;

    /**
     * @brief reset the id of this solver
     */
    void setId( const std::string& id );

private:

    /**
     * @brief The ID of this solver.
     */
    std::string mId;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /**
     * @brief The solver logger, should never be NULL pointer
     */
    LoggerPtr mLogger;
};

/** Define the shared pointer class for Solver */

typedef std::shared_ptr<_Solver> _SolverPtr;

/* ---------------------------------------------------------------------- */
/*   Implementation of inline methods                                     */
/* ---------------------------------------------------------------------- */

const std::string& _Solver::getId() const
{
    return mId;
}

} /* end namespace solver */

} /* end namespace scai */
