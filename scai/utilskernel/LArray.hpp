/**
 * @file LArray.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
            ValueType tmp = HArrayUtils::getVal<ValueType>( mArray, mIndex );
            HArrayUtils::setVal( mArray, mIndex, tmp );
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
            HArrayUtils::setScalar( *this, value, common::reduction::COPY, context );
        }
        else
        {
            HArrayUtils::setScalar( *this, value, common::reduction::COPY, hmemo::Context::getHostPtr() );
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
        HArrayUtils::assignOp( *this, other, common::reduction::MULT );
        return *this;
    }

    LArray& operator*= ( const ValueType val )
    {
        HArrayUtils::setScalar( *this, val, common::reduction::MULT );
        return *this;
    }

    LArray& operator/= ( const hmemo::_HArray& other )
    {
        HArrayUtils::assignOp( *this, other, common::reduction::DIVIDE );
        return *this;
    }

    LArray& operator/= ( const ValueType val )
    {
        HArrayUtils::setScalar( *this, val, common::reduction::DIVIDE );
        return *this;
    }

    LArray& operator+= ( const hmemo::_HArray& other )
    {
        HArrayUtils::assignOp( *this, other, common::reduction::ADD );
        return *this;
    }

    LArray& operator+= ( const ValueType val )
    {
        HArrayUtils::setScalar( *this, val, common::reduction::ADD );
        return *this;
    }

    LArray& operator-= ( const hmemo::_HArray& other )
    {
        HArrayUtils::assignOp( *this, other, common::reduction::SUB );
        return *this;
    }

    LArray& operator-= ( const ValueType val )
    {
        HArrayUtils::setScalar( *this, val, common::reduction::SUB );
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
        HArrayUtils::setScalar( *this, val, common::reduction::COPY, this->getFirstTouchContextPtr() );
        return *this;
    }

    IndexProxy operator[] ( const IndexType i )
    {
        return IndexProxy( *this, i );
    }

    ValueType operator[] ( const IndexType i ) const
    {
        return HArrayUtils::getVal( *this, i );
    }

    /** Get the minimal value of an array */

    ValueType min() const
    {
        return HArrayUtils::reduce( *this, common::reduction::MIN );
    }

    /** Get the maximal value of an array */

    ValueType max() const
    {
        return HArrayUtils::reduce( *this, common::reduction::MAX );
    }

    /** Get the maximal value of an array */

    ValueType maxNorm() const
    {
        return HArrayUtils::reduce( *this, common::reduction::ABS_MAX );
    }

    /** Get the sum of all array elements */

    ValueType sum() const
    {
        return HArrayUtils::reduce( *this, common::reduction::ADD );
    }

    /** Compute the sum of magnitudes, for complex numbers it is the sum of real and imag part */

    ValueType l1Norm() const
    {
        return HArrayUtils::asum( *this );
    }

    ValueType l2Norm() const
    {
        return common::Math::sqrt( HArrayUtils::dotProduct( *this, *this ) );
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
        HArrayUtils::invert( *this );
    }

    /** Compute the conj in-place */

    void conj()
    {
        HArrayUtils::conj( *this );
    }

};

} /* end namespace utilskernel */

} /* end namespace scai */
