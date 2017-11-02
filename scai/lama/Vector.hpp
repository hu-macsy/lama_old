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

#include <memory>
#include <scai/common/TypeTraits.hpp>

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

    /** Create a new vector of a certain kind but with same value type */

    static Vector<ValueType>* getVector( VectorKind kind );

    /** Desctructor. */

    virtual ~Vector();

    /** Overwrite _Vector::newVector to get the covariant return type */

    virtual Vector<ValueType>* newVector( void ) const = 0;

    /** Overwrite _Vector::copy to get the covariant return type */

    virtual Vector<ValueType>* copy( void ) const = 0;

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

    /**
     * @brief Returns the L1 norm of this.
     *
     * @return the L1 norm of this.
     *
     * l1Norm computes the sum of the absolute values of this.
     */
    virtual NormType<ValueType> l1Norm() const = 0;

    /**
     * @brief Returns the L2 norm of this.
     *
     * @return the L2 norm of this.
     *
     * l2Norm computes the sum of the absolute values of this.
     */
    virtual NormType<ValueType> l2Norm() const = 0;

    /**
     * @brief Returns the max norm of this.
     *
     * @return the max norm of this.
     *
     * maxNorm computes the value of this with the largest magnitude.
     */
    virtual NormType<ValueType> maxNorm() const = 0;

    /**
     * @brief Returns the max norm of the difference with another vector
     *
     *  v1.maxDiffNorm( v2 ) is equivalent to:
     *
     *  \code
     *      Vector<ValueType> tmp = v1 - v2;
     *      maxNorm( tmp )
     *  \endcode
     *
     *  But it avoids the temporary vector wherever possible
     */
    virtual NormType<ValueType> maxDiffNorm( const _Vector& other ) const = 0;

    /**
     * @brief Returns the global minimum value of this.
     *
     * @return   the global minimum value of this vector.
     */
    virtual ValueType min() const = 0;

    /**
     * @brief Returns the global maximum value of this.
     *
     * @return the global maximum value of this vector.
     */
    virtual ValueType max() const = 0;

    /**
     * @brief Returns the sum of all vector elements.
     *
     * @return the sum of all vector elements.
     *
     * As the summation of the values depends on the mapping of the values to
     * the processors, this routine might return slightly different results
     * for different parallel environments.
     */
    virtual ValueType sum() const = 0;

    /**
     * @brief Returns the dot product of this and other.
     *
     * @param[in] other   the vector to calculate the dot product with.
     * @return            the dot product of this and other
     */
    virtual ValueType dotProduct( const _Vector& other ) const = 0;

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

    // Implementations of pure _Vector methods to guarantee upward compatibilty

    Scalar _l1Norm() const;

    Scalar _l2Norm() const;

    Scalar _maxNorm() const;

    Scalar _maxDiffNorm( const _Vector& other ) const;

    Scalar _sum() const;
    Scalar _min() const;
    Scalar _max() const;

    Scalar _dotProduct( const _Vector& other ) const;
};
  
/** 
 * Definiton of corresponding shared pointer type for the class Vector<ValueType> by a type alias.
 *
 *  \code
 *      VectorPtr<ValueType> x( Vector<ValueType>::getVector( VectorKind::SPARSE ) );
 *      std::shared_ptr<Vector<ValueType> > x( Vector<ValueType>::getVector( VectorKind::DENSE ) );
 *  \endcode
*/
template<typename ValueType>
using VectorPtr = std::shared_ptr<Vector<ValueType> >;

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
