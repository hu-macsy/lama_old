/**
 * @file LArray.hpp
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
 * @brief Definition of LArray as HArray with additional kernel functionality.
 * @author Thomas Brandes
 * @date 18.11.2015
 */
#pragma once

#include <scai/utilskernel/HArrayUtils.hpp>

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/Math.hpp>

namespace scai
{

namespace utilskernel
{

/** The LArray is derived form HArray but offers some more functionality.
 *
 *  - constructor with array of values
 *  - copy operator with arbitrary array (implicit conversions)
 *  - assignment operator with scalar (initialization on any device)
 *  - assignment operator with arbitrary array (implicit conversions)
 *
 *  In contrary to the HArray the LArray can only be instantiated for
 *  arithmetic types and uses already some kernels of the kernel
 *  library (initialization, conversion, scale, ... )
 */

template<typename ValueType>
class LArray : public hmemo::HArray<ValueType>
{
public:

    /** Help class to observe the further use of operator[] in LArray */
    class IndexProxy
    {
    public:

        /** Proxy constructed by ref to the array and the index value. */

        IndexProxy( LArray<ValueType>& array, const IndexType i ) :

            mArray( array ),
            mIndex( i )
        {
        }

        /** indexed value proxy can be used to get its value */

        operator ValueType() const
        {
            return HArrayUtils::getVal<ValueType>( mArray, mIndex );
        }

        /** indexed value proxy can be assigned a value */

        IndexProxy& operator= ( ValueType val )
        {
            HArrayUtils::setVal( mArray, mIndex, val );
            return *this;
        }

        /** Override the default assignment operator to avoid ambiguous interpretation of a[i] = b[i] */

        IndexProxy& operator= ( const IndexProxy& other )
        {
            ValueType tmp = HArrayUtils::getVal<ValueType>( other.mArray, other.mIndex );
            HArrayUtils::setVal( mArray, mIndex, tmp );
            return *this;
        }

    private:

        LArray<ValueType>& mArray;
        IndexType mIndex;

    };

    /* -------------------------------------------------------------------------------------- */
    /* -   Constructors                                                                       */
    /* -------------------------------------------------------------------------------------- */

    /** Default constructor with no arguments. */

    LArray()
    {
    }

    /** Construct an array and give it a first touch */

    explicit LArray( hmemo::ContextPtr context ) :

        hmemo::HArray<ValueType>( context )

    {
    }

    /** @brief Construcor with context and size.
     *
     *  This constructor also allocates memory on the specified context.
     */
    explicit LArray( const IndexType n,
                     hmemo::ContextPtr context = hmemo::ContextPtr() )
    {
        this->resize( n );    // only sets the size, no allocate

        if ( context.get() != NULL )
        {
            this->reserve( context, n );  // includes also the first touch
        }
    }

    /** @brief Construcor with context and size and initial value.
     *
     *  @param n is the size of the array
     *  @param value is the initial value
     *  @param context is location where initialization is done, if not specified its the host
     */

    explicit LArray( const IndexType n,
                     const ValueType value,
                     hmemo::ContextPtr context = hmemo::ContextPtr() )
    {
        // SCAI_ASSERT( context.get(), "NULL context" )
        this->resize( n );  // size of the array must be known before a value can be assigned

        // context == NULL might happen by DenseVector

        if ( context.get() )
        {
            HArrayUtils::setScalar( *this, value, binary::COPY, context );
        }
        else
        {
            HArrayUtils::setScalar( *this, value, binary::COPY, hmemo::Context::getHostPtr() );
        }
    }

    /** @brief Construcor with context and size and initial startValue and increment.
     *
     *  @param n is the size of the array
     *  @param startValue is the initial value
     *  @param inc is the increment for the sequence of values
     *  @param context is location where initialization is done, if not specified its the host
     */

    explicit LArray( const IndexType n,
                     const ValueType startValue,
                     const ValueType inc,
                     hmemo::ContextPtr context = hmemo::ContextPtr() )
    {
        // SCAI_ASSERT( context.get(), "NULL context" )
        this->resize( n );  // size of the array must be known before a value can be assigned

        // context == NULL might happen by DenseVector

        if ( context.get() )
        {
            HArrayUtils::setSequence( *this, startValue, inc, n, context );
        }
        else
        {
            HArrayUtils::setSequence( *this, startValue, inc, n, hmemo::Context::getHostPtr() );
        }
    }

    /** @brief Construcor with size and initial values
     *
     *  @param n is the size of the array
     *  @param values is array with values, must be at least n
     *  @param context is location for which the array is valid afterwards
     */

    template<typename OtherValueType>
    explicit LArray( const IndexType n,
                     const OtherValueType* values,
                     hmemo::ContextPtr context = hmemo::Context::getHostPtr() )
    {
        hmemo::HArrayRef<OtherValueType> tmp( n, values );
        HArrayUtils::assign( *this, tmp, context );
    }

    /** Override the default copy construtor */

    LArray( const LArray<ValueType>& other ) : hmemo::HArray<ValueType>()
    {
        HArrayUtils::assign( *this, other );
    }

    /** Copy constructor that works with HArray of any type. */

    LArray( const hmemo::_HArray& other ) : hmemo::HArray<ValueType>()
    {
        HArrayUtils::assign( *this, other );
    }

    /** Copy constructor that works with HArray of any type and specifies context */

    LArray( const hmemo::_HArray& other, hmemo::ContextPtr context ) : hmemo::HArray<ValueType>()
    {
        HArrayUtils::assign( *this, other, context );
    }

    LArray& operator= ( const LArray<ValueType>& other )
    {
        HArrayUtils::assign( *this, other );
        return *this;
    }

    LArray& operator= ( const hmemo::_HArray& other )
    {
        HArrayUtils::assign( *this, other );
        return *this;
    }

    LArray& operator*= ( const hmemo::_HArray& other )
    {
        HArrayUtils::setArray( *this, other, binary::MULT );
        return *this;
    }

    LArray& operator*= ( const ValueType val )
    {
        HArrayUtils::setScalar( *this, val, binary::MULT );
        return *this;
    }

    LArray& operator/= ( const hmemo::_HArray& other )
    {
        HArrayUtils::setArray( *this, other, binary::DIVIDE );
        return *this;
    }

    LArray& operator/= ( const ValueType val )
    {
        HArrayUtils::setScalar( *this, val, binary::DIVIDE );
        return *this;
    }

    LArray& operator+= ( const hmemo::_HArray& other )
    {
        HArrayUtils::setArray( *this, other, binary::ADD );
        return *this;
    }

    LArray& operator+= ( const ValueType val )
    {
        HArrayUtils::setScalar( *this, val, binary::ADD );
        return *this;
    }

    LArray& operator-= ( const hmemo::_HArray& other )
    {
        HArrayUtils::setArray( *this, other, binary::SUB );
        return *this;
    }

    LArray& operator-= ( const ValueType val )
    {
        HArrayUtils::setScalar( *this, val, binary::SUB );
        return *this;
    }

    /** Assignment operator to initialize an array with a certain value.
     *
     *  @brief val is value to be assigned.
     *
     *  The size of the array remains unchanged and should have been set before.
     *  Initialization is done on the first context of the array.
     */
    LArray& operator= ( const ValueType val )
    {
        //  assignment is done on the first touch memory/context
        hmemo::ContextPtr context = this->getFirstTouchContextPtr();
        SCAI_ASSERT( context.get(), "No first touch context" )
        HArrayUtils::setScalar( *this, val, binary::COPY, this->getFirstTouchContextPtr() );
        return *this;
    }

    IndexProxy operator[] ( const IndexType i )
    {
        return IndexProxy( *this, i );
    }

    ValueType operator[] ( const IndexType i ) const
    {
        return HArrayUtils::getVal<ValueType>( *this, i );
    }

    void setRandom( IndexType n, float fillRate = 1.0f, hmemo::ContextPtr context = hmemo::ContextPtr() )
    {
        HArrayUtils::setRandomImpl( *this, n, fillRate, context );
    }

    /** Get the minimal value of an array */

    ValueType min() const
    {
        return HArrayUtils::reduce( *this, binary::MIN );
    }

    /** Get the maximal value of an array */

    ValueType max() const
    {
        return HArrayUtils::reduce( *this, binary::MAX );
    }

    /** Get the maximal value of an array */

    ValueType maxNorm() const
    {
        return HArrayUtils::reduce( *this, binary::ABS_MAX );
    }

    /** Get the sum of all array elements */

    ValueType sum() const
    {
        return HArrayUtils::reduce( *this, binary::ADD );
    }

    /** Compute the sum of magnitudes, for complex numbers it is the sum of real and imag part */

    ValueType l1Norm() const
    {
        return HArrayUtils::asum( *this );
    }

    ValueType l2Norm() const
    {
        return HArrayUtils::nrm2( *this );
    }

    /** Build the max diff norm with another LAMA array */

    ValueType maxDiffNorm( const hmemo::HArray<ValueType>& other ) const
    {
        return HArrayUtils::absMaxDiffVal( *this, other );
    }

    /** Compute the dotproduct with another LAMA array. */

    ValueType dotProduct( const hmemo::HArray<ValueType>& other ) const
    {
        return HArrayUtils::dotProduct( *this, other );
    }

    /** Compute the inverse/reciprocal in-place, all elements should be non-zero  */

    void invert()
    {
        // use the binary op DIVIDE that is also supported
        HArrayUtils::binaryOpScalar1( *this, ValueType( 1 ), *this, binary::DIVIDE );
    }

    /** Compute the conj in-place */

    void conj()
    {
        HArrayUtils::unaryOp( *this, *this, unary::CONJ );
    }

    void exp()
    {
        HArrayUtils::unaryOp( *this, *this, unary::EXP );
    }

    void log()
    {
        HArrayUtils::unaryOp( *this, *this, unary::LOG );
    }

    void floor()
    {
        HArrayUtils::unaryOp( *this, *this, unary::FLOOR );
    }

    void ceil()
    {
        HArrayUtils::unaryOp( *this, *this, unary::CEIL );
    }

    void sqrt()
    {
        HArrayUtils::unaryOp( *this, *this, unary::SQRT );
    }

    void sin()
    {
        HArrayUtils::unaryOp( *this, *this, unary::SIN );
    }

    void cos()
    {
        HArrayUtils::unaryOp( *this, *this, unary::COS );
    }

    void tan()
    {
        HArrayUtils::unaryOp( *this, *this, unary::TAN );
    }

    void atan()
    {
        HArrayUtils::unaryOp( *this, *this, unary::ATAN );
    }

    void powBase( ValueType base )
    {
        HArrayUtils::binaryOpScalar1( *this, base, *this, binary::POW );
    }

    void powExp( ValueType exp )
    {
        HArrayUtils::binaryOpScalar2( *this, *this, exp, binary::POW );
    }

    void powBase( const hmemo::HArray<ValueType>& base )
    {
        return HArrayUtils::binaryOp( *this, base, *this, binary::POW );
    }

    void powExp( const hmemo::HArray<ValueType>& exp )
    {
        return HArrayUtils::binaryOp( *this, *this, exp, binary::POW );
    }
};

} /* end namespace utilskernel */

} /* end namespace scai */
