/**
 * @file _DenseVector.hpp
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
 * @brief Definition of common base clase for dense vectors.
 * @author Thomas Brandes
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/Vector.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Common base class for all dense vectors to deal with untyped dense vectors.
 *
 * A dense vector is a distributed one-dimensional array where each value is explicitly available.
 * This base class provides common methods for typed dense vectors and allows access to the
 * data via untyped heterogeneous arrays.
 */

class COMMON_DLL_IMPORTEXPORT _DenseVector :

    public Vector

{

public:

    /** Implementation of pure method Vector::getVectorKind() */

    inline virtual VectorKind getVectorKind() const;

    /**
     * @brief get a reference to local values of this dense vector.
     *
     * @return  a reference to the local values of this dense vector.
     *
     * This method allows to modifiy the values of the array but it 
     * should not be used to change the size of the array.
     *
     * Derived classes can overwrite this method and might use a covariant return type,
     * i.e. a typed version of a heterogeneous array.
     */

    virtual hmemo::_HArray& getLocalValues() = 0;

    /**
     * @brief get a constant reference to local values of this Dense Vector.
     *
     * @return  a constant reference to the local values of this dense vector.
     */
    virtual const hmemo::_HArray& getLocalValues() const = 0;

    /**
     * @brief This method fills/initializes an allocated vector 
     *
     * @param[in] startValue value for the first element
     * @param[in] inc increment between the elements
     *
     * This is a pure routine that is implemented individually by derived classes.
     */
    virtual void fillRange( const Scalar startValue, const Scalar inc ) = 0;

    /**
     * @brief allocate a replicated vector and fill it with range values
     *
     * @param[in] n becomes the size of the vector
     * @param[in] startValue value for the first element, $vector[0] = startValue$
     * @param[in] inc increment between two elements, i.e. $vector[i+1] - vector[i] = inc$
     */
    void setRange( const IndexType n, const Scalar startValue, const Scalar inc )
    {
        allocate( n );
        fillRange( startValue, inc );
    }

    /**
     * This method initializes a (distributed) dense vector with a sequence of values
     *
     * @param[in] distribution determines global/local size of the vector
     * @param[in] startValue value for the first elemen
     * @param[in] inc increment between the element
     */
    void setRange( dmemo::DistributionPtr distribution, const Scalar startValue, const Scalar inc )
    {
        allocate( distribution );
        fillRange( startValue, inc );
    }

    /**
     * @brief Create a new dense vector of a certain type 
     *
     * @param type is the value type of the vector
     * @return new allocated DenseVector of the given type 
     * @throw  common::Exception if no dense vector of this type is registered in factory.
     */
     
    static _DenseVector* create( common::ScalarType );

    // make operators and methods of Vector visible for _DenseVector

    using Vector::operator=;

    /**
     * @brief Implementation of pure method Vector::isConsistent 
     */
    virtual bool isConsistent() const;

protected:

    // All constructors here just call the corresponing constructors of Vector 

    _DenseVector( const IndexType n );

    _DenseVector( const IndexType n, hmemo::ContextPtr context );

    _DenseVector( const dmemo::DistributionPtr dist );

    _DenseVector( const dmemo::DistributionPtr dist, hmemo::ContextPtr context );

    _DenseVector( const _DenseVector& other );

    _DenseVector( const Vector& other );

    /** Common logger for all types of dense vectors. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* ------------------------------------------------------------------------- */

Vector::VectorKind _DenseVector::getVectorKind() const
{
    return Vector::DENSE;
}

} /* end namespace lama */

} /* end namespace scai */
