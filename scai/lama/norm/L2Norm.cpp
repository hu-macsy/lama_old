/**
 * @file L2Norm.cpp
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
 * @brief L2Norm.cpp
 * @author Jiri Kraus
 * @date 01.06.2011
 */

// hpp
#include <scai/lama/norm/L2Norm.hpp>

template scai::lama::Norm::Register<scai::lama::L2Norm>::RegisterGuard
scai::lama::Norm::Register<scai::lama::L2Norm>::registerGuard;

namespace scai
{

namespace lama
{

L2Norm::L2Norm()
{
}

L2Norm::~L2Norm()
{
}

std::string L2Norm::createValue()
{
   return "L2";
}

Norm* L2Norm::create()
{
   return new L2Norm();
}

void L2Norm::writeAt( std::ostream& stream ) const
{   
    stream << "L2Norm";
}

Scalar L2Norm::apply( const Scalar& scalar ) const
{
    return l2Norm( scalar );
}

Scalar L2Norm::apply( const Vector& vector ) const
{
    return l2Norm( vector );
}

Scalar L2Norm::apply( const Matrix& matrix ) const
{
    return l2Norm( matrix );
}

Scalar l2Norm( const Scalar& scalar )
{
    return abs( scalar );
}

Scalar l2Norm( const Vector& vector )
{
    return vector.l2Norm();
}

Scalar l2Norm( const Matrix& matrix )
{
    return matrix.l2Norm();
}

} /* end namespace lama */

} /* end namespace scai */
