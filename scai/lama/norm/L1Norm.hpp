/**
 * @file L1Norm.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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

class COMMON_DLL_IMPORTEXPORT L1Norm: public scai::lama::Norm
{
public:
    L1Norm();
    virtual ~L1Norm();

    virtual Scalar apply( const Scalar& scalar ) const;

    virtual Scalar apply( const Vector& vector ) const;

	virtual Scalar apply( const Matrix& matrix ) const;
};

Scalar l1Norm( const Scalar& scalar );

Scalar l1Norm( const Vector& vector );

Scalar l1Norm( const Matrix& matrix );

} /* end namespace lama */

} /* end namespace scai */
