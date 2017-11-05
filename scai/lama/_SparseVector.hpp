/**
 * @file _SparseVector.hpp
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
 * @brief Definition of template class that stands for a sparse vector of a certain type.
 * @author Thomas Brandes
 * @date 16.01.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/Vector.hpp>

#include <scal/hmemo/HArray.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Common base class for all sparse vectors to deal with untyped sparse vectors.
 *
 * A sparse vector is a distributed one-dimensional array where only non-zero values are explicitly stored.
 * This base class provides common methods for typed sparse vectors and allows access to the
 * data via untyped heterogeneous arrays.
 */

class COMMON_DLL_IMPORTEXPORT _SparseVector :

    public Vector

{

public:

    /** Implementation of pure method Vector::getVectorKind() */

    inline virtual VectorKind getVectorKind() const;

    /**
     * @brief get a constant reference to local non-zero values of this Sparse Vector.
     *
     * @return  a constant reference to the local non-zero values of this dense vector.
     */
    virtual const hmemo::_HArray& getNonZeroValues() const = 0;

    virtual const hmemo::HArray<IndexType>& getNonZeroIndexes() const = 0;

    /** Query the zero value, i.e. default value at positions not in nonZeroIndexes. */

    virtual Scalar getZero() const = 0;

    /**
     * @brief Implementation of pure method Vector::isConsistent 
     */
    virtual bool isConsistent() const;

    /**
     * @brief Create a new sparse vector of a certain type 
     *
     * @param type is the value type of the vector
     * @return new allocated SparseVector of the given type 
     * @throw  common::Exception if no dense vector of this type is registered in factory.
     */

    static _SparseVector* create( common::ScalarType );

    // make operators and methods of Vector visible for _SparseVector

    // make operators and methods of Vector visible for _SparseVector

    using Vector::operator=;

protected:

    // All constructors here just call the corresponing constructors of Vector 

    _SparseVector( const IndexType n );

    _SparseVector( const IndexType n, hmemo::ContextPtr context );

    _SparseVector( const dmemo::DistributionPtr dist );

    _SparseVector( const dmemo::DistributionPtr dist, hmemo::ContextPtr context );

    _SparseVector( const _SparseVector& other );

    _SparseVector( const Vector& other );

    /** Common logger for all types of sparse vectors. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};


Vector::VectorKind _SparseVector::getVectorKind() const
{
    return Vector::SPARSE;
}

} /* end namespace lama */

} /* end namespace scai */
