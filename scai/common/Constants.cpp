/**
 * @file Constants.cpp
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
 * @brief Class for defining constants
 * @author Eric Schricker
 * @date 28.04.2015
 * @since 1.1.0
 */

#include <scai/common/Constants.hpp>
#include <scai/common/SCAITypes.hpp>

#include <scai/common/Complex.hpp>

#include <limits>

#include <boost/preprocessor.hpp>
#include <boost/scoped_array.hpp>

namespace scai
{

namespace common
{

template<typename ValueType>
const ValueType Constants<ValueType>::generateEps()
{
//    ValueType eps = 1.0;
//
//    while( eps + (ValueType) 1.0 > 1.0 )
//        eps *= 0.5;
//
//    eps *= 64.0;
//
//    return eps;
	return std::numeric_limits<ValueType>::epsilon();
}

template<>
const ComplexFloat Constants<ComplexFloat>::generateEps()
{
	return ComplexFloat( std::numeric_limits<float>::epsilon() );
}

template<>
const ComplexDouble Constants<ComplexDouble>::generateEps()
{
	return ComplexDouble( std::numeric_limits<double>::epsilon() );
}

template<typename ValueType>
const ValueType Constants<ValueType>::generateSfmin()
{
    ValueType sfmin = std::numeric_limits<ValueType>::min();
    ValueType small = ValueType(1) / std::numeric_limits<ValueType>::max();

    if( small >= sfmin )
    {
    	sfmin = small * ( ValueType( 1.0 ) + std::numeric_limits<ValueType>::epsilon() );
    }

    return sfmin;
}

template<>
const ComplexFloat Constants<ComplexFloat>::generateSfmin()
{
    float sfmin = std::numeric_limits<float>::min();
    float small = float(1) / std::numeric_limits<float>::max();

    if( small >= sfmin )
    {
    	sfmin = small * ( float( 1.0 ) + std::numeric_limits<float>::epsilon() );
    }

    return ComplexFloat( sfmin );
}

template<>
const ComplexDouble Constants<ComplexDouble>::generateSfmin()
{
	double sfmin = std::numeric_limits<double>::min();
	double small = double(1) / std::numeric_limits<double>::max();

    if( small >= sfmin )
    {
    	sfmin = small * ( double( 1.0 ) + std::numeric_limits<double>::epsilon() );
    }

    return ComplexDouble( sfmin );
}

template<typename ValueType> const ValueType Constants<ValueType>::eps = Constants<ValueType>::generateEps();
template<typename ValueType> const ValueType Constants<ValueType>::sfmin = Constants<ValueType>::generateSfmin();

template<typename ValueType> const ValueType Constants<ValueType>::one = ValueType(1.0);
template<typename ValueType> const ValueType Constants<ValueType>::zero = ValueType(0.0);
template<typename ValueType> const ValueType Constants<ValueType>::minusone = ValueType(-1.0);

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

#define COMMON_CONSTANTS_INSTANTIATE(z, I, _)                                      \
        template class COMMON_DLL_IMPORTEXPORT Constants<ARRAY_TYPE##I>;

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, COMMON_CONSTANTS_INSTANTIATE, _ )

#undef COMMON_CONSTANTS_INSTANTIATE

} /* end namespace common */

} /* end namespace scai */
