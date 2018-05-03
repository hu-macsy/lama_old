/**
 * @file examples/DemoMath.cpp
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
 * @brief Examples for using TypeTraits
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#ifdef SCAI_COMPLEX_SUPPORTED
#include <scai/common/Complex.hpp>
#endif

#include <scai/common/Math.hpp>

#include <iostream>

using scai::common::Math;

template<typename ValueType>
void testRoutine( ValueType x )
{
    ValueType absX = Math::abs( x );
    std::cout << "abs ( " << x << " ) = " << absX << std::endl;
    ValueType sqrtX = Math::sqrt( x );
    std::cout << "sqrt ( " << x << " ) = " << sqrtX << std::endl;
    ValueType conjX = Math::conj( x );
    std::cout << "conj ( " << x << " ) = " << conjX << std::endl;
}

int main()
{
    testRoutine<float>( -1 );
    testRoutine<double>( -1 );
#ifdef SCAI_COMPLEX_SUPPORTED
    testRoutine<scai::ComplexFloat>( scai::ComplexFloat( 2, -1 ) );
#endif
}
