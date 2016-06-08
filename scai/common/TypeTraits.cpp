/**
 * @file TypeTraits.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Type specific implementations for arithmetic value types.
 * @author Thomas Brandes
 * @date 17.11.2015
 */


#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace common
{

const scalar::ScalarType TypeTraits<IndexType>::stype;
const scalar::ScalarType TypeTraits<float>::stype;
const scalar::ScalarType TypeTraits<double>::stype;
const scalar::ScalarType TypeTraits<long double>::stype;

#ifdef SCAI_COMPLEX_SUPPORTED
const scalar::ScalarType TypeTraits<ComplexFloat>::stype;
const scalar::ScalarType TypeTraits<ComplexDouble>::stype;
const scalar::ScalarType TypeTraits<ComplexLongDouble>::stype;
#endif

}  // namespace common

}  // namespace scai

