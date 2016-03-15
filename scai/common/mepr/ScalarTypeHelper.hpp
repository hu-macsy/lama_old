/**
 * @file common/mepr/ScalarTypeHelper.hpp
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
 * @brief Conversion from ScalarType to ValueType and other way round
 * @author Eric Schricker
 * @date 10.03.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai {

namespace common {

namespace mepr {

/*
 * Forward declartion
 */

template<typename TList>
struct ScalarTypeHelper;

/*
 * Termination
 */
template<>
struct ScalarTypeHelper<NullType>
{
    static long sizeOf( const scalar::ScalarType& )
    {
        return 0;
    }

    static scalar::ScalarType getBySize( const long )
    {
        return scalar::UNKNOWN;
    }
};

/*
 * Step n
 */
template<typename H, typename T>
struct ScalarTypeHelper< TypeList<H,T> >
{
    static long sizeOf( const scalar::ScalarType& s )
    {
        if( s == common::getScalarType<H>() )
        {
            return sizeof( H );
        }
        else
        {
            return ScalarTypeHelper<T>::sizeOf( s );
        }
    }

    static scalar::ScalarType getBySize( const long size )
    {
        if( sizeof( H ) == size )
        {
            return TypeTraits<H>::stype;
        }
        else
        {
            return ScalarTypeHelper< T >::getBySize( size );
        }
    }
};

} /* end namespace mepr */

} /* end namespace common */

} /* end namespace scai */
