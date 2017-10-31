/**
 * @file Vector.hpp
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
 * @brief Definition of an abstract class for distributed vectors of a given type
 * @author Thomas Brandes
 * @date 30.10.2017
 */
#pragma once

#include <scai/lama/_Vector.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Definition of an abstract class that represents a distributed one-dimensional vector 
 *        of a certain value type.
 *
 * @tparam ValueType stands for the type of the entries in the distributed vector.
 *
 * Sparse and dense vector class derive from this abstract base class.
 */
template <typename ValueType>
class COMMON_DLL_IMPORTEXPORT Vector: public _Vector

{
public:

    /** Desctructor. */

    virtual ~Vector();

    /** Implementation of _Vector::getValueType */

    virtual common::scalar::ScalarType getValueType() const;

    using _Vector::getContext;
    using _Vector::getDistribution;

protected:

    /**
     * Constructor of replicated vector for derived classes by size and/or context
     *
     * @param[in] size    number of entries for the vector
     * @param[in] context is optional, will be Host context.
     *
     * Note: this constructor overrides also the default constructor.
     */
    explicit Vector( const IndexType size = 0, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Constructor of Vector for derived classes by distribution
     *
     * @param[in] distribution  the distribution to use for the new Vector.
     * @param[in] context       is optional, will be Host context.
     */
    explicit Vector( dmemo::DistributionPtr distribution, hmemo::ContextPtr context = hmemo::ContextPtr() );

    /**
     * @brief Creates a copy of the passed _Vector.
     *
     * @param[in] other   the Vector to take a copy from.
     *
     * Inherits size/distribution, context and content of the passed vector
     */
    Vector( const _Vector& other );

    /** Override the default copy constructor */

    Vector( const Vector<ValueType>& other );
};
  
} /* end namespace lama */

} /* end namespace scai */
