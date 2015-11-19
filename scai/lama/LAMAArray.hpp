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
 *  arithmetic types and uses already some kernels of the Kernel
 *  library (initialization, conversion, scale, ... )
 */

template<typename ValueType> 
class LAMAArray : public hmemo::HArray<ValueType>
{
public:

    LAMAArray() : hmemo::HArray<ValueType>()
    {
    }

    explicit LAMAArray( hmemo::ContextPtr context ) : hmemo::HArray<ValueType>( context )
    {
    }

    explicit LAMAArray( hmemo::ContextPtr context, const IndexType n ) : hmemo::HArray<ValueType>( context )
    {
        this->resize( n );  // Note: this does not yet allocate memory
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
    }

    LAMAArray& operator= ( const hmemo::ContextArray& other )
    {
        HArrayUtils::assign( *this, other );
    }

    LAMAArray& operator= ( const ValueType val )
    {
        //  assignment is done on the first touch memory/context

        HArrayUtils::assignScalar( *this, val, this->getFirstTouchContextPtr() );
    }
};

} /* end namespace lama */

} /* end namespace scai */
