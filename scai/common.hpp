/**
 * @file common.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief common.hpp
 * @author Thomas Brandes
 * @date 16.06.2015
 */

#pragma once

#pragma message("It is not recommended to include whole common")

#include <scai/common/bind.hpp>
#include <scai/common/Complex.hpp>
#include <scai/common/config.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/ContextType.hpp>
#include <scai/common/Factory1.hpp>
#include <scai/common/Factory.hpp>
#include <scai/common/function.hpp>
#include <scai/common/LibModule.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/NonCopyable.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/preprocessor.hpp>
#include <scai/common/Printable.hpp>
#include <scai/common/ReductionOp.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/shared_ptr.hpp>
#include <scai/common/Thread.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/weak_ptr.hpp>

// CUDA
#include <scai/common/cuda/CUDACallable.hpp>

// Exceptions
#include <scai/common/exception/AssertException.hpp>
#include <scai/common/exception/Exception.hpp>
#include <scai/common/exception/NotSupportedValueTypeException.hpp>
#include <scai/common/exception/UnsupportedException.hpp>

// Macros
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/inline.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/typeloop.hpp>
#include <scai/common/macros/unique_name.hpp>
#include <scai/common/macros/unused.hpp>

// MIC
#include <scai/common/mic/MICCallable.hpp>
