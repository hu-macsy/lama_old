/**
 * @file utilskernel/external/FFTW3Types.hpp
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
 * @brief Configuration of found fftw3 library types
 * @author Lauretta Schubert
 * @date 03.11.2016
 */
#pragma once

#define SCAI_NUMERIC_TYPES_FFTW3 float, double

#define SCAI_NUMERIC_TYPES_FFTW3_LIST SCAI_TYPELIST( SCAI_NUMERIC_TYPES_FFTW3 )
