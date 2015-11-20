/**
 * @file LAMAArray.hpp
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
 * @brief Definition of LAMAArray as HArray with more functionality
 * @author Thomas Brandes
 * @date 18.11.2015
 */
#pragma once

#include <scai/lama/HArrayUtils.hpp>

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace lama
{

/** The LAMAArray is derived form HArray but offers some more functionality. 
 *
 *  - copy operator with arbitrary array (implicit conversions)
 *  - assignment operator with scalar (initialization on any device)
 *  - assignment operator with arbitrary array (implicit conversions)
 *
 *  In contrary to the HArray the LAMAArray can only be instantiated for
 *  arithmetic types and uses already some kernels of the kernel
 *  library (initialization, conversion, scale, ... )
 */

template<typename ValueType> 
class LAMAArray : public hmemo::HArray<ValueType>
{
public:

    LAMAArray()
    {
    }

    explicit LAMAArray( hmemo::ContextPtr context ) :

        hmemo::HArray<ValueType>( context )
   
    {
    }

    /** @brief Construcor with context and size.
     *
     *  This constructor only allocates memory on the given device.
     */
    explicit LAMAArray( const IndexType n,
                        hmemo::ContextPtr context = hmemo::ContextPtr() )
    {
        this->resize( n );            // reserve does not include the resize.

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

    explicit LAMAArray( const IndexType n, 
                        const ValueType value, 
                        hmemo::ContextPtr context = hmemo::Context::getHostPtr() ) 
    {
        // SCAI_ASSERT( context.get(), "NULL context" )

        this->resize( n );  // size of the array must be known before a value can be assigned

        // context == NULL might happen by DenseVector

        if ( context.get() )
        {
            HArrayUtils::assignScalar( *this, value, context );
        }
        else
        {
            HArrayUtils::assignScalar( *this, value, hmemo::Context::getHostPtr() );
        }
    }

    /** @brief Construcor with size and initial values
     *
     *  @param n is the size of the array
     *  @param values is array with values, must be at least n
     *  @param context is location for which the array is valid afterwards
     */

    template<typename OtherValueType>
    explicit LAMAArray( const IndexType n, 
                        const OtherValueType* values,
                        hmemo::ContextPtr context = hmemo::Context::getHostPtr() ) 
    {
        hmemo::HArrayRef<OtherValueType> tmp( n, values );
        this->reserve( context, n );  // includes also the first touch
        HArrayUtils::assign( *this, tmp );
    }

    LAMAArray( const LAMAArray<ValueType>& other ) : hmemo::HArray<ValueType>()
    {
        HArrayUtils::assign( *this, other );
    }

    LAMAArray( const hmemo::ContextArray& other ) : hmemo::HArray<ValueType>()
    {
        HArrayUtils::assign( *this, other );
    }

    LAMAArray& operator= ( const LAMAArray<ValueType>& other )
    {
        HArrayUtils::assign( *this, other );
        return *this;
    }

    LAMAArray& operator= ( const hmemo::ContextArray& other )
    {
        HArrayUtils::assign( *this, other );
        return *this;
    }

    LAMAArray& operator= ( const ValueType val )
    {
        //  assignment is done on the first touch memory/context

        hmemo::ContextPtr context = this->getFirstTouchContextPtr();
        SCAI_ASSERT( context.get(), "No first touch context" )
        HArrayUtils::assignScalar( *this, val, this->getFirstTouchContextPtr() );
        return *this;
    }
};

} /* end namespace lama */

} /* end namespace scai */
