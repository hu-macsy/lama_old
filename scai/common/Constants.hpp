/**
 * @file Constants.hpp
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

#pragma once

// for dll_import
#include <scai/common/config.hpp>

namespace scai
{

namespace common
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Constants {

public:
	static const ValueType eps;
	static const ValueType sfmin;
	static const ValueType zero;
	static const ValueType one;
	static const ValueType minusone;

private:
	Constants(){}
	virtual ~Constants(){}

	static const ValueType generateEps();
	static const ValueType generateSfmin();
};

} /* end namespace common */

} /* end namespace scai */
