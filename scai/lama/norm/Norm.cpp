/**
 * @file Norm.cpp
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
 * @brief Norm.cpp
 * @author Jiri Kraus
 * @date 01.06.2011
 */

// hpp
#include <scai/lama/norm/Norm.hpp>

#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
Norm<ValueType>::Norm()
{
}

template<typename ValueType>
Norm<ValueType>::~Norm()
{
}

template<typename ValueType>
RealType<ValueType> Norm<ValueType>::operator()( const ValueType& scalar ) const
{
    return apply( scalar );
}

template<typename ValueType>
RealType<ValueType> Norm<ValueType>::operator()( const Vector<ValueType>& vector ) const
{
    return apply( vector );
}

template<typename ValueType>
RealType<ValueType> Norm<ValueType>::operator()( const Matrix<ValueType>& matrix ) const
{
    return apply( matrix );
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Norm, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
