/**
 * @file SolutionProxy.hpp
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
 * @brief Contains the Interface of the Class SolutionProxy
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/lama/Vector.hpp>

namespace scai
{

namespace solver
{

/**
 * @brief The SolutionProxy is used to avoid multiple residual calculations
 *
 * The SolutionProxy is used to avoid needless residual re-calculation.
 * Basically it is a wraps around pointer to a vector. Each not read-only
 * access to the underlying vector marks the proxy as dirty and by that signals
 * the solver to recalculate the residual.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT SolutionProxy
{
public:

    /**
     * @brief Creates a SolutionProxy with no associated vector.
     */
    SolutionProxy();

    /**
     * @brief Creates a SolutionProxy for the given pointer
     *
     * SolutionProxy does not take owner ship of the passed pointer. Therefore the
     * caller has to make sure that the associated vector exists as long as this
     * SolutionProxy exists.
     *
     * @param[in] solution   The pointer to the vector which the proxy will wrap.
     */
    SolutionProxy( lama::Vector<ValueType>* solution );

    /**
     * @brief SolutionProxy destructor.
     *
     * The destructor ~SolutionProxy() does not destroy the associated vector.
     */
    ~SolutionProxy();

    /**
     * @brief Returns a constant reference to the underlying vector, dirty flag remains unchanged
     *
     * @return constant reference to the underlying vector.
     */
    const lama::Vector<ValueType>& getConstReference() const;

    /**
     * @brief Associates the given Vector Pointer with this SolutionProxy.
     *
     * The assignment drops the probably existing association to anther Vector.
     *
     * @param[in] newVector the Vector to which the SolutionProxy shall point to.
     */
    void operator=( lama::Vector<ValueType>* newVector );

    /**
     * @brief Determines if the proxy is dirty and the residual needs to be
     *        recomputed.
     *
     * @return Value to determine if the proxy is dirty.
     */
    bool isDirty() const;

    /**
     * @brief Sets/unsets the isDirtY flag of the SolutionProxy
     *
     * @param[in] isDirty   value determining whether the proxy is dirty or not.
     */
    void setDirty( bool isDirty );

    /**
     * @brief Returns a reference to the underlying vector, marks it as dirty
     *
     * @return Reference to the underlying vector.
     */
    lama::Vector<ValueType>& getReference();

private:

    /**
     * @brief The underlying solution vector.
     */
    lama::Vector<ValueType>* mSolution;

    /**
     * @brief Flag which determines, if the Proxy is dirty or not.
     */
    bool mIsDirty;
};

} /* end namespace solver */

} /* end namespace scai */
