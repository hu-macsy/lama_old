/**
 * @file L2Norm.cpp
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
 * @brief L2Norm.cpp
 * @author Jiri Kraus
 * @date 01.06.2011
 */

// hpp
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
L2Norm<ValueType>::L2Norm()
{
}

template<typename ValueType>
L2Norm<ValueType>::~L2Norm()
{
}

template<typename ValueType>
std::string L2Norm<ValueType>::createValue()
{
    return "L2";
}

template<typename ValueType>
Norm<ValueType>* L2Norm<ValueType>::create()
{
    return new L2Norm();
}

template<typename ValueType>
void L2Norm<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "L2Norm";
}

template<typename ValueType>
RealType<ValueType> L2Norm<ValueType>::apply( const ValueType& scalar ) const
{
    return l2Norm( scalar );
}

template<typename ValueType>
RealType<ValueType> L2Norm<ValueType>::apply( const Vector<ValueType>& vector ) const
{
    return l2Norm( vector );
}

template<typename ValueType>
RealType<ValueType> L2Norm<ValueType>::apply( const Matrix<ValueType>& matrix ) const
{
    return l2Norm( matrix );
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( L2Norm, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
