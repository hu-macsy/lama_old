/**
 * @file OpenMPInterface.cpp
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
 * @brief Implementation of methods for class OpenMPInterface
 * @author Thomas Brandes
 * @date 18.05.2011
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPInterface.hpp>

#include <lama/openmp/OpenMPUtils.hpp>
#include <lama/openmp/OpenMPCSRUtils.hpp>
#include <lama/openmp/OpenMPELLUtils.hpp>
#include <lama/openmp/OpenMPJDSUtils.hpp>
#include <lama/openmp/OpenMPDIAUtils.hpp>
#include <lama/openmp/OpenMPCOOUtils.hpp>
#include <lama/openmp/OpenMPDenseUtils.hpp>

#include <lama/openmp/OpenMPBLAS1.hpp>
#include <lama/openmp/OpenMPBLAS2.hpp>
#include <lama/openmp/OpenMPBLAS3.hpp>
#include <lama/openmp/OpenMPLAPACK.hpp>
#include <lama/openmp/OpenMPSCALAPACK.hpp>

#include <lama/LAMAInterfaceRegistry.hpp>

// boost
#include <boost/bind.hpp>

#include <cmath>

using std::auto_ptr;
using boost::bind;
using boost::ref;
using boost::cref;
using std::sqrt;

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPInterface::logger, "LAMAInterface.OpenMPInterface" )

OpenMPInterface::OpenMPInterface() : LAMAInterface()
{
    // set typed function pointers with the corresponding OpenMP implementations

    OpenMPBLAS1::setInterface( BLAS );
    OpenMPBLAS2::setInterface( BLAS );
    OpenMPBLAS3::setInterface( BLAS );
    OpenMPLAPACK::setInterface( BLAS );

    // Use of SCALAPACK is optional in LAMA

#ifdef LAMA_MKL_SCALAPACK
    OpenMPSCALAPACK::setInterface( BLAS );
#endif

    OpenMPUtils::setInterface( Utils );

    OpenMPCSRUtils::setInterface( CSRUtils );
    OpenMPELLUtils::setInterface( ELLUtils );
    OpenMPJDSUtils::setInterface( JDSUtils );
    OpenMPDIAUtils::setInterface( DIAUtils );
    OpenMPCOOUtils::setInterface( COOUtils );
    OpenMPDenseUtils::setInterface( DenseUtils );
}

OpenMPInterface::~OpenMPInterface()
{
}

LAMA_LAMAINTERFACE_REGISTRATION( Context::Host, OpenMPInterface )

} // namespace
