/**
 * @file L2Norm.hpp
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
 * @brief L2Norm.hpp
 * @author Jiri Kraus
 * @date 01.06.2011
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

/**
 * @brief L2Norm is a functor specialization of Norm that calculates the l2 norm
 *        for the passed values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT L2Norm:

    public Norm<ValueType>,
    public Norm<ValueType>::template Register<L2Norm<ValueType> >

{
public:
    /**
     * @brief Constructs a L2Norm functor.
     */
    L2Norm();

    /**
     * @brief Destroys this L2Norm functor.
     */
    virtual ~L2Norm();

    /**
     * @brief calculates the l2 norm of the passed Scalar.
     *
     * @param[in] scalar    the scalar to calculate the l2 norm of.
     *
     * @return              the l2 norm of the Scalar scalar.
     */
    virtual RealType<ValueType> apply( const ValueType& scalar ) const;

    /**
     * @brief calculates the l2 norm of the passed Vector.
     *
     * @param[in] vector    the Vector to calculate the l2 norm of.
     *
     * @return              the l2 norm of the Vector vector.
     */
    virtual RealType<ValueType> apply( const Vector<ValueType>& vector ) const;

    /**
     * @brief calculates the l2 norm of the passed matrix.
     *
     * @param[in] matrix    the _Matrix to calculate the l2 norm of.
     *
     * @return              the l2 norm of matrix.
     */
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

/**
 * @brief calculates the l2 norm of the passed Scalar.
 *
 * @param[in] scalar    the Scalar to calculate the l2 norm of.
 *
 * @return              the l2 norm of the Scalar scalar.
 */
template<typename ValueType>
inline RealType<ValueType> l2Norm( const ValueType& scalar )
{
    return common::Math::abs( scalar );
}

/**
 * @brief calculates the l2 norm of the passed Vector.
 *
 * @param[in] vector    the Vector to calculate the l2 norm of.
 *
 * @return              the l2 norm of the Vector vector.
 */
template<typename ValueType>
inline RealType<ValueType> l2Norm( const Vector<ValueType>& vector )
{
    return vector.l2Norm();
}

/**
 * @brief calculates the l2 norm of the passed matrix.
 *
 * @param[in] matrix    the _Matrix to calculate the l2 norm of.
 *
 * @return              the l2 norm of the _Matrix matrix.
 */
template<typename ValueType>
inline RealType<ValueType> l2Norm( const Matrix<ValueType>& matrix )
{
    return matrix.l2Norm();
}

} /* end namespace lama */

} /* end namespace scai */
