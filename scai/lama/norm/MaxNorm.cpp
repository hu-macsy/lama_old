/**
 * @file MaxNorm.cpp
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
 * @brief Implementations for derived norm class MaxNorm.
 * @author Thomas Brandes, Jiri Kraus
 * @date 14.06.2011
 */

// hpp
#include <scai/lama/norm/MaxNorm.hpp>


namespace scai
{

template lama::Norm::Register<lama::MaxNorm>::RegisterGuard
lama::Norm::Register<lama::MaxNorm>::registerGuard;

namespace lama
{

MaxNorm::MaxNorm()
{
}

MaxNorm::~MaxNorm()
{
}

std::string MaxNorm::createValue()
{
   return "Max";
}

Norm* MaxNorm::create()
{
   return new MaxNorm();
}

void MaxNorm::writeAt( std::ostream& stream ) const
{   
    stream << "MaxNorm";
}

Scalar MaxNorm::apply( const Scalar& scalar ) const
{
    return maxNorm( scalar );
}

Scalar MaxNorm::apply( const Vector& vector ) const
{
    return maxNorm( vector );
}

Scalar MaxNorm::apply( const Matrix& matrix ) const
{
    return maxNorm( matrix );
}

Scalar maxNorm( const Scalar& scalar )
{
    return abs( scalar );
}

Scalar maxNorm( const Vector& vector )
{
    return vector.maxNorm();
}

Scalar maxNorm( const Matrix& matrix )
{
    return matrix.maxNorm();
}

} /* end namespace lama */

} /* end namespace scai */
