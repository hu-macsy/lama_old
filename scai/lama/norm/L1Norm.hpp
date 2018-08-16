/**
 * @file L1Norm.hpp
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
 * @brief L1Norm.hpp
 * @author Jiri Kraus
 * @date 14.06.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/norm/Norm.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT L1Norm:

    public Norm<ValueType>,
    public Norm<ValueType>::template Register<L1Norm<ValueType> >

{
public:

    L1Norm();

    virtual ~L1Norm();

    virtual RealType<ValueType> apply( const ValueType& scalar ) const;

    virtual RealType<ValueType> apply( const Vector<ValueType>& vector ) const;

    virtual RealType<ValueType> apply( const Matrix<ValueType>& matrix ) const;

    static std::string createValue();

    static Norm<ValueType>* create();

    /** Override Printable::writeAt with version for this class. */

    virtual void writeAt( std::ostream& stream ) const;
};

template<typename ValueType>
inline RealType<ValueType> l1Norm( const ValueType& scalar )
{
    return common::Math::abs( scalar );
}

template<typename ValueType>
inline RealType<ValueType> l1Norm( const Vector<ValueType>& vector )
{
    return vector.l1Norm();
}

template<typename ValueType>
inline RealType<ValueType> l1Norm( const Matrix<ValueType>& matrix )
{
    return matrix.l1Norm();
}


} /* end namespace lama */

} /* end namespace scai */
