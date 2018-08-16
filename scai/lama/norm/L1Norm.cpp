/**
 * @file L1Norm.cpp
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
 * @brief L1Norm.cpp
 * @author Jiri Kraus
 * @date 14.06.2011
 */

// hpp
#include <scai/lama/norm/L1Norm.hpp>

#include <scai/common/macros/instantiate.hpp>


namespace scai
{

namespace lama
{

template<typename ValueType>
L1Norm<ValueType>::L1Norm()
{
}

template<typename ValueType>
L1Norm<ValueType>::~L1Norm()
{
}

template<typename ValueType>
std::string L1Norm<ValueType>::createValue()
{
    return "L1";
}

template<typename ValueType>
Norm<ValueType>* L1Norm<ValueType>::create()
{
    return new L1Norm();
}

template<typename ValueType>
void L1Norm<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "L1Norm";
}

template<typename ValueType>
RealType<ValueType> L1Norm<ValueType>::apply( const ValueType& scalar ) const
{
    return l1Norm( scalar );
}

template<typename ValueType>
RealType<ValueType> L1Norm<ValueType>::apply( const Vector<ValueType>& vector ) const
{
    return l1Norm( vector );
}

template<typename ValueType>
RealType<ValueType> L1Norm<ValueType>::apply( const Matrix<ValueType>& matrix ) const
{
    return l1Norm( matrix );
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( L1Norm, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
