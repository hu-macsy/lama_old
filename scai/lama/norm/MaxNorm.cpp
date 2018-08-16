/**
 * @file MaxNorm.cpp
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
 * @brief Implementations for derived norm class MaxNorm.
 * @author Thomas Brandes, Jiri Kraus
 * @date 14.06.2011
 */

// hpp
#include <scai/lama/norm/MaxNorm.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
MaxNorm<ValueType>::MaxNorm()
{
}

template<typename ValueType>
MaxNorm<ValueType>::~MaxNorm()
{
}

template<typename ValueType>
std::string MaxNorm<ValueType>::createValue()
{
    return "Max";
}

template<typename ValueType>
Norm<ValueType>* MaxNorm<ValueType>::create()
{
    return new MaxNorm();
}

template<typename ValueType>
void MaxNorm<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "MaxNorm";
}

template<typename ValueType>
RealType<ValueType> MaxNorm<ValueType>::apply( const ValueType& scalar ) const
{
    return maxNorm( scalar );
}

template<typename ValueType>
RealType<ValueType> MaxNorm<ValueType>::apply( const Vector<ValueType>& vector ) const
{
    return maxNorm( vector );
}

template<typename ValueType>
RealType<ValueType> MaxNorm<ValueType>::apply( const Matrix<ValueType>& matrix ) const
{
    return maxNorm( matrix );
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( MaxNorm, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
