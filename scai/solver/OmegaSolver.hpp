/**
 * @file OmegaSolver.hpp
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
 * @brief OmegaSolver.hpp
 * @author Kai Buschulte
 * @date 10.08.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/solver/IterativeSolver.hpp>

// logging
#include <scai/logging.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief The OldSolutionHandler class only manages the omega parameter
 * For solvers like a Jacobi.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT OmegaSolver: public IterativeSolver<ValueType>
{
public:

    /**
     * @brief Creates a solver with the given id.
     *
     * @param[in] id The id of the solver.
     */
    OmegaSolver( const std::string& id );

    /**
     * @brief Creates a solver with the given id and omega.
     *
     * @param[in] id    The id of the solver.
     * @param[in] omega The omega parameter which is used by the jacobi solver.
     */
    OmegaSolver( const std::string& id, const ValueType omega );

    /**
     * @brief Creates a solver with the given id and logger.
     *
     * @param[in] id     The id of the solver.
     * @param[in] logger The logger which is used by this solver.
     */
    OmegaSolver( const std::string& id, LoggerPtr logger );

    /**
     * @brief Creates a solver with the given id, omega and logger.
     *
     * @param[in] id     The id of the solver.
     * @param[in] omega  The omega parameter which is used by the jacobi solver.
     * @param[in] logger The logger used by the solver.
     */
    OmegaSolver( const std::string& id, const ValueType omega, LoggerPtr logger );

    /**
     * @brief Copy constructor that copies the status independent solver information
     */
    OmegaSolver( const OmegaSolver& other );

    /**
     * @brief Destructor.
     */
    virtual ~OmegaSolver();

    /**
     * @brief Sets the omega parameter of this.
     *
     * @param[in] omega The omega parameter of the omega solver.
     */
    void setOmega( const ValueType omega );

    /**
     * @brief Returns omega.
     *
     * @return Omega.
     */
    ValueType getOmega() const;

    /**
     * @brief Override the pure copy method to get a covariant return type
     */
    virtual OmegaSolver* copy() = 0;

protected:

    ValueType mOmega;
};

template<typename ValueType>
OmegaSolver<ValueType>::OmegaSolver( const std::string& id ) : 

    IterativeSolver<ValueType>( id ), 
    mOmega( 0.5 )
{
}

template<typename ValueType>
OmegaSolver<ValueType>::OmegaSolver( const std::string& id, const ValueType omega ) : 

    IterativeSolver<ValueType>( id ), 
    mOmega( omega )
{
}

template<typename ValueType>
OmegaSolver<ValueType>::OmegaSolver( const std::string& id, LoggerPtr logger ) : 

     IterativeSolver<ValueType>( id, logger ), 
     mOmega( 0.5 )
{
}

template<typename ValueType>
OmegaSolver<ValueType>::OmegaSolver( const std::string& id, const ValueType omega, LoggerPtr logger ) : 

    IterativeSolver<ValueType>( id, logger ), 
    mOmega( omega )
{
}

template<typename ValueType>
OmegaSolver<ValueType>::OmegaSolver( const OmegaSolver<ValueType>& other ) : 

    IterativeSolver<ValueType>( other ), 
    mOmega( other.mOmega )
{
}

template<typename ValueType>
OmegaSolver<ValueType>::~OmegaSolver()
{
}

template<typename ValueType>
void OmegaSolver<ValueType>::setOmega( const ValueType omega )
{
    mOmega = omega;
}

template<typename ValueType>
ValueType OmegaSolver<ValueType>::getOmega() const
{
    return mOmega;
}

} /* end namespace solver */

} /* end namespace scai */
