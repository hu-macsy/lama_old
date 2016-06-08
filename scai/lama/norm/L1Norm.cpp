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
 * @endlicense
 *
 * @brief L1Norm.cpp
 * @author Jiri Kraus
 * @date 14.06.2011
 */

// hpp
#include <scai/lama/norm/L1Norm.hpp>

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
