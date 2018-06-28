/**
 * @file MaxNorm.hpp
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
 * @brief MaxNorm.hpp
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
class COMMON_DLL_IMPORTEXPORT MaxNorm:

    public Norm<ValueType>,
    public Norm<ValueType>::template Register<MaxNorm<ValueType> >

{
public:

    MaxNorm();

    virtual ~MaxNorm();

    virtual RealType<ValueType> apply( const ValueType& scalar ) const;

    virtual RealType<ValueType> apply( const Vector<ValueType>& vector ) const;

    virtual RealType<ValueType> apply( const Matrix<ValueType>& matrix ) const;

    /**
     *  Getter routine for key of this derived class used in Norm factory
     */
    static std::string createValue();

    /**
     *  Create method is just function version of constructor.
     */
    static Norm<ValueType>* create();

    /** Override Printable::writeAt with version for this class. */

    virtual void writeAt( std::ostream& stream ) const;
};

template<typename ValueType>
inline RealType<ValueType> maxNorm( const ValueType& scalar )
{
    return common::Math::abs( scalar );
}


template<typename ValueType>
inline RealType<ValueType> maxNorm( const Vector<ValueType>& vector )
{
    return vector.maxNorm();
}

template<typename ValueType>
inline RealType<ValueType> maxNorm( const Matrix<ValueType>& matrix )
{
    return matrix.maxNorm();
}

} /* end namespace lama */

} /* end namespace scai */
