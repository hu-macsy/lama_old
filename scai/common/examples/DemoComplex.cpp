/**
 * @file common/examples/DemoComplex.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Using complex numbers of scai/common
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#include <scai/common/Complex.hpp>
#include <scai/common/Math.hpp>

#include <iostream>

using scai::common::Complex;

int main()
{
    Complex<float> a( 1.0f, 1.0f );
    Complex<double> b( 2.0, -1.0 );
    std::cout << "a = " << a << ", b = " << b << std::endl;
    double b_re = b.real();
    double b_im = b.imag();
    std::cout << "b.real() = " << b_re << ", b.imag() = " << b_im << std::endl;
    std::cout << "conj( b ) = " << scai::common::Math::conj( b ) << std::endl;
    // not allowed for std::complex
    std::cout << "static_cast<double>( b ) = " << static_cast<double>( b ) << std::endl;
    std::cout << "abs( b ) = " << scai::common::Math::abs( b ) << std::endl;

    if ( a == b )
    {
        std::cout << a << " == " << b << std::endl;
    }
    else if ( a != b )
    {
        std::cout << a << " == " << b << std::endl;
    }
    else
    {
        std::cout << "a = " << a << " and b = " << b << " are neihter equal nor not equal" << std::endl;
    }
}
