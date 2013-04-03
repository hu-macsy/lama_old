/**
 * @file CUDAInterface.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief CUDAInterface.cpp
 * @author Jiri Kraus
 * @date 26.05.2011
 * $Id$
 */

// hpp
#include <lama/cuda/CUDAInterface.hpp>

// others
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAUtils.hpp>
#include <lama/cuda/CUDACSRUtils.hpp>
#include <lama/cuda/CUDAELLUtils.hpp>
#include <lama/cuda/CUDAJDSUtils.hpp>
#include <lama/cuda/CUDADIAUtils.hpp>
#include <lama/cuda/CUDACOOUtils.hpp>

#include <lama/cuda/CUDABLAS1.hpp>
#include <lama/cuda/CUDABLAS2.hpp>
#include <lama/cuda/CUDABLAS3.hpp>
#include <lama/cuda/CUDALAPACK.hpp>

using std::auto_ptr;

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDAInterface::logger, "LAMAInterface.CUDAInterface" )

CUDAInterface::CUDAInterface() : LAMAInterface()
{
    // set typed function pointers with the corresponding CUDA implementations

    CUDABLAS1::setInterface( BLAS );
    CUDABLAS2::setInterface( BLAS );
    CUDABLAS3::setInterface( BLAS );
    CUDALAPACK::setInterface( BLAS );

    // Use of SCALAPACK is not supported in CUDA, no entries required

    CUDAUtils::setInterface( Utils );

    CUDACSRUtils::setInterface( CSRUtils );
    CUDAELLUtils::setInterface( ELLUtils );
    CUDAJDSUtils::setInterface( JDSUtils );
    CUDADIAUtils::setInterface( DIAUtils );
    CUDACOOUtils::setInterface( COOUtils );

    CUDABLAS1::setInterface( BLAS );
}

CUDAInterface::~CUDAInterface()
{
}

// register an incarnation of CUDAInterface at the LAMA Interface registry that is used for the CUDA Context

LAMA_LAMAINTERFACE_REGISTRATION( Context::CUDA, CUDAInterface )

} // namespace

