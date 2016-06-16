/**
 * @file L1Norm.cpp
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
 * @brief L1Norm.cpp
 * @author Jiri Kraus
 * @date 14.06.2011
 */

// hpp
#include <scai/lama/norm/L1Norm.hpp>

template scai::lama::Norm::Register<scai::lama::L1Norm>::RegisterGuard 
scai::lama::Norm::Register<scai::lama::L1Norm>::registerGuard;

namespace scai
{

namespace lama
{

L1Norm::L1Norm()
{
}

L1Norm::~L1Norm()
{
}

std::string L1Norm::createValue()
{
   return "L1";
}

Norm* L1Norm::create()
{
   return new L1Norm();
}

void L1Norm::writeAt( std::ostream& stream ) const
{
    stream << "L1Norm";
}

Scalar L1Norm::apply( const Scalar& scalar ) const
{
    return l1Norm( scalar );
}

Scalar L1Norm::apply( const Vector& vector ) const
{
    return l1Norm( vector );
}

Scalar L1Norm::apply( const Matrix& matrix ) const
{
    return l1Norm( matrix );
}

Scalar l1Norm( const Scalar& scalar )
{
    return abs( scalar );
}

Scalar l1Norm( const Vector& vector )
{
    return vector.l1Norm();
}

Scalar l1Norm( const Matrix& matrix )
{
    return matrix.l1Norm();
}

} /* end namespace lama */

} /* end namespace scai */
