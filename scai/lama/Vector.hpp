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

#include <scai/common/shared_ptr.hpp>

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

    /** Definiton of corresponding shared pointer type for this class 
     *
     *  \code
     *      Vector<ValueType>::Ptr x( Vector<ValueType>::getVector( VectorKind::SPARSE ) );
     *      std::shared_ptr<Vector<ValueType> > x( Vector<ValueType>::getVector( VectorKind::DENSE ) );
     *  \endcode
     */
    typedef common::shared_ptr<Vector<ValueType> > Ptr;

    /** Create a new vector of a certain kind but with same value type */

    static Vector<ValueType>* getVector( VectorKind kind );

    /** Desctructor. */

    virtual ~Vector();

    /** Implementation of _Vector::getValueType */

    virtual common::scalar::ScalarType getValueType() const;

    using _Vector::getContext;
    using _Vector::getDistribution;
    using _Vector::operator=;

    /** Help class to observe the further use of operator[] for Vector */

    class VectorElemProxy
    {
    public:

        /** Proxy constructed by ref to the array and the index value. */

        inline VectorElemProxy( Vector<ValueType>& vector, const IndexType i );

        /** Proxy for a vector element can be used to get its value, type conversion to ValueType
         *
         *  @returns current value of the vector element as a single value 
         */
        inline operator ValueType() const;

        /** indexed value proxy can be assigned a scalar */

        inline VectorElemProxy& operator= ( ValueType val );

        /** Override the default assignment operator to avoid ambiguous interpretation of a[i] = b[i] */

        inline VectorElemProxy& operator= ( const VectorElemProxy& other );

    private:

        Vector<ValueType>& mVector;
        IndexType mIndex;
    };

    /**
     *  Indexing of a distributed vector returns a proxy so that this operator can be used
     *  on lhs and rhs of an assignment.
     */
    VectorElemProxy operator[]( const IndexType i )
    {
        return VectorElemProxy( *this, i );
    }

    /**
     *  Indexing of a const distributed vector returns directly the corresponding element.
     */
    ValueType operator[]( const IndexType i ) const
    {
        Scalar s = getValue( i );
        return s.getValue<ValueType>();
    }

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
  
/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>::VectorElemProxy::VectorElemProxy( Vector<ValueType>& vector, const IndexType i ) :

    mVector( vector ),
    mIndex( i )

{
}

template<typename ValueType>
Vector<ValueType>::VectorElemProxy::operator ValueType() const
{
    Scalar s = mVector.getValue( mIndex );
    return s.getValue<ValueType>();
}

template<typename ValueType>
typename Vector<ValueType>::VectorElemProxy& Vector<ValueType>::VectorElemProxy::operator= ( ValueType val )
{
    mVector.setValue( mIndex, Scalar( val ) );
    return *this;
}

template<typename ValueType>
typename Vector<ValueType>::VectorElemProxy& Vector<ValueType>::VectorElemProxy::operator= ( const Vector<ValueType>::VectorElemProxy& other )
{
    Scalar tmp = other.mVector.getValue( other.mIndex );
    mVector.setValue( mIndex, tmp );
    return *this;
}

} /* end namespace lama */

} /* end namespace scai */
