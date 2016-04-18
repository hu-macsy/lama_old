/**
 * @file lama/mepr/CRTPMatrixStorageWrapper.hpp
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
 * @brief Wrapper for templated calls in Solver
 * @author Eric Schricker
 * @date 14.04.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/ScalarType.hpp>

#include <scai/lama/Scalar.hpp>

namespace scai 
{
	
namespace solver
{

namespace mepr 
{

template<typename TList>
struct SolverEps;

template<>
struct SolverEps<common::mepr::NullType>
{
	static lama::Scalar get( const common::scalar::ScalarType& )
	{
		return lama::Scalar( 0.0 );
	}
};

template<typename H, typename T>
struct SolverEps<common::mepr::TypeList<H, T> >
{
	static lama::Scalar get( const common::scalar::ScalarType& type )
	{
		if( common::TypeTraits<H>::stype == type )
		{
			return lama::Scalar( common::TypeTraits<H>::getEps() );
		}
		else
		{
			return SolverEps<T>::get( type );
		}
	}
};


} /* end namespace mepr */
} /* end namespace solver */
} /* end namepsace scai */