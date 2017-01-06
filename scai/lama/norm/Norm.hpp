/**
 * @file Norm.hpp
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
 * @brief Norm.hpp
 * @author Jiri Kraus
 * @date 01.06.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/Vector.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/common/Factory.hpp>

namespace scai
{

namespace lama
{

class Norm;

typedef common::shared_ptr<Norm> NormPtr;

/**
 * @brief Norm is a abstract functor to calculate the norms for the passed values.
 *
 * The Functor Norm is mainly used by the stopping criteria ResidualThreshold and
 * ResidualStagnation, to allow a to customize the Norms for these stopping criteria.
 */
class COMMON_DLL_IMPORTEXPORT Norm : 

   public common::Factory<std::string, Norm*>,
   public common::Printable

{
public:
    /**
     * @brief Constructs a empty norm object.
     */
    Norm();

    /**
     * @brief Destroys this.
     */
    virtual ~Norm();

    /**
     * @brief Calculates the norm of the passed Scalar.
     *
     * To call this is equivalent to call apply with the same argument.
     *
     * @param[in] scalar    the Scalar to caluclate the norm for.
     * @return              the norm of scalar.
     */
    Scalar operator()( const Scalar& scalar ) const;

    /**
     * @brief Calculates the norm of the passed Vector.
     *
     * To call this is equivalent to call apply with the same argument.
     *
     * @param[in] vector    the Vector to caluclate the norm for.
     * @return              the norm of vector.
     */
    Scalar operator()( const Vector& vector ) const;

    /**
     * @brief Calculates the norm of the passed Matrix.
     *
     * To call this is equivalent to call apply with the same argument.
     *
     * @param[in] matrix    the Matrix to caluclate the norm for.
     * @return              the norm of matrix.
     */
    Scalar operator()( const Matrix& matrix ) const;

    /**
     * @brief Calculates the norm of the passed Scalar.
     *
     * @param[in] scalar    the Scalar to caluclate the norm for.
     * @return              the norm of scalar.
     */
    virtual Scalar apply( const Scalar& scalar ) const = 0;

    /**
     * @brief Calculates the norm of the passed Vector.
     *
     * @param[in] vector    the Vector to caluclate the norm for.
     * @return              the norm of vector.
     */
    virtual Scalar apply( const Vector& vector ) const = 0;

    /**
     * @brief Calculates the norm of the passed Matrix.
     *
     * @param[in] matrix    the Matrix to caluclate the norm for.
     * @return              the norm of matrix.
     */
    virtual Scalar apply( const Matrix& matrix ) const = 0;
};

} /* end namespace lama */

} /* end namespace scai */
