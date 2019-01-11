/**
 * @file scai/common.hpp
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
 * @brief Single header file that includes all header files of the common subproject.
 * @author Thomas Brandes
 * @date 16.06.2015
 */

#pragma once

#pragma message("It is not recommended to include whole common")

#ifdef SCAI_COMPLEX_SUPPORTED
#include <scai/common/Complex.hpp>
#endif

#include <scai/common/config.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/ContextType.hpp>
#include <scai/common/Factory1.hpp>
#include <scai/common/Factory.hpp>
#include <scai/common/LibModule.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/NonCopyable.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/Printable.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/thread.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Walltime.hpp>

// CUDA
#include <scai/common/cuda/CUDACallable.hpp>

// Exceptions
#include <scai/common/exception/AssertException.hpp>
#include <scai/common/exception/Exception.hpp>
#include <scai/common/exception/InvalidArgumentException.hpp>
#include <scai/common/exception/UnsupportedException.hpp>

// Macros
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/inline.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/loop.hpp>
#include <scai/common/macros/unique_name.hpp>
#include <scai/common/macros/unused.hpp>
