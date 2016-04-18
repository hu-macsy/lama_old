/**
 * @file HArrays.hpp
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
 * @brief Vector with all HArray one for each supported type
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
 */

#include <scai/hmemo/_HArray.hpp>

#include <vector>

typedef scai::common::shared_ptr<scai::hmemo::_HArray> ArrayPtr;

/** Class for a list of matrix storage pointers, one for each supported 
 *  matrix storage format and each supported arithmetic type.
 */

class HArrays : public std::vector<ArrayPtr> 
{

public:

    /** Constructor creates already the list with all storage pointers. */

    HArrays( scai::hmemo::ContextPtr ctx = scai::hmemo::ContextPtr() ) 
    {
        using namespace scai::common;
        using namespace scai::hmemo;

        std::vector<scalar::ScalarType> values;  //  all create values

        _HArray::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            ArrayPtr arrayPtr( _HArray::create( values[i] ) );

            if ( ctx )
            {
                arrayPtr->prefetch( ctx );
            }

            push_back( arrayPtr );
        }
    }

    // Destructor will free all arrays due to use of shared pointers
};
