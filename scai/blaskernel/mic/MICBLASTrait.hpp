/**
 * @file MICBLASTrait.hpp
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
 * @brief Definitions for MKL interface on mic
 * @author Eric Stricker
 * @date 21.01.2016
 */

#pragma once

// macros
#define MIC_BLAS_NAME( name, prefix ) prefix##name

#define MIC_BLAS_CALL( name, prefix, ... )	\
		SCAI_MIC_CALL( CUBLAS_BLAS_NAME( name, prefix ), __VAR_ARGS__ )

// external
#include <mkl_blas.h>

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT MICBLASTrait
{
public:
	typedef int BLASIndexType;
	typedef char BLASTrans;
};

} /* end namespace blaskernel */

} /* end namespace scai */
