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
 * @brief OpenMPInterface.cpp
 * @author Jiri Kraus
 * @date 18.05.2011
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPInterface.hpp>

// others
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
#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>

#include <lama/task/TaskSyncToken.hpp>

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

OpenMPInterface::OpenMPInterface()
    : LAMAInterface()
{
    /** BLAS */

    /** Level 1 */

    /** float */
    mFloatBLAS1Interface.scal = &OpenMPBLAS1::scal<float>;
    mFloatBLAS1Interface.nrm2 = &OpenMPBLAS1::nrm2<float>;
    mFloatBLAS1Interface.asum = &OpenMPBLAS1::asum<float>;
    mFloatBLAS1Interface.iamax = &OpenMPBLAS1::iamax<float>;
    mFloatBLAS1Interface.viamax = &OpenMPBLAS1::viamax<float>;
    mFloatBLAS1Interface.swap = &OpenMPBLAS1::swap<float>;
    mFloatBLAS1Interface.copy = &OpenMPBLAS1::copy<float>;
    mFloatBLAS1Interface.axpy = &OpenMPBLAS1::axpy<float>;
    mFloatBLAS1Interface.dot = &OpenMPBLAS1::dot<float>;
    mFloatBLAS1Interface.sum = &OpenMPBLAS1::sum<float>;
    mFloatBLAS1Interface.rot = &OpenMPBLAS1::rot<float>;
    mFloatBLAS1Interface.rotm = &OpenMPBLAS1::rotm<float>;
    mFloatBLAS1Interface.rotg = &OpenMPBLAS1::rotg<float>;
    mFloatBLAS1Interface.rotmg = &OpenMPBLAS1::rotmg<float>;
    mFloatBLAS1Interface.ass = &OpenMPBLAS1::ass<float>;

    /** double */
    mDoubleBLAS1Interface.scal = &OpenMPBLAS1::scal<double>;
    mDoubleBLAS1Interface.nrm2 = &OpenMPBLAS1::nrm2<double>;
    mDoubleBLAS1Interface.asum = &OpenMPBLAS1::asum<double>;
    mDoubleBLAS1Interface.iamax = &OpenMPBLAS1::iamax<double>;
    mDoubleBLAS1Interface.viamax = &OpenMPBLAS1::viamax<double>;
    mDoubleBLAS1Interface.swap = &OpenMPBLAS1::swap<double>;
    mDoubleBLAS1Interface.copy = &OpenMPBLAS1::copy<double>;
    mDoubleBLAS1Interface.axpy = &OpenMPBLAS1::axpy<double>;
    mDoubleBLAS1Interface.dot = &OpenMPBLAS1::dot<double>;
    mDoubleBLAS1Interface.sum = &OpenMPBLAS1::sum<double>;
    mDoubleBLAS1Interface.rot = &OpenMPBLAS1::rot<double>;
    mDoubleBLAS1Interface.rotm = &OpenMPBLAS1::rotm<double>;
    mDoubleBLAS1Interface.rotg = &OpenMPBLAS1::rotg<double>;
    mDoubleBLAS1Interface.rotmg = &OpenMPBLAS1::rotmg<double>;
    mDoubleBLAS1Interface.ass = &OpenMPBLAS1::ass<double>;

    /** Level 2 */

    /** float */
    mFloatBLAS2Interface.gemv = &OpenMPBLAS2::gemv<float>;
    mFloatBLAS2Interface.symv = &OpenMPBLAS2::symv<float>;
    mFloatBLAS2Interface.trmv = &OpenMPBLAS2::trmv<float>;
    mFloatBLAS2Interface.trsv = &OpenMPBLAS2::trsv<float>;
    mFloatBLAS2Interface.gbmv = &OpenMPBLAS2::gbmv<float>;
    mFloatBLAS2Interface.sbmv = &OpenMPBLAS2::sbmv<float>;
    mFloatBLAS2Interface.tbmv = &OpenMPBLAS2::tbmv<float>;
    mFloatBLAS2Interface.tbsv = &OpenMPBLAS2::tbsv<float>;
    mFloatBLAS2Interface.ger = &OpenMPBLAS2::ger<float>;
    mFloatBLAS2Interface.syr = &OpenMPBLAS2::syr<float>;
    mFloatBLAS2Interface.syr2 = &OpenMPBLAS2::syr2<float>;
    mFloatBLAS2Interface.spmv = &OpenMPBLAS2::spmv<float>;
    mFloatBLAS2Interface.spr = &OpenMPBLAS2::spr<float>;
    mFloatBLAS2Interface.spr2 = &OpenMPBLAS2::spr2<float>;
    mFloatBLAS2Interface.tpmv = &OpenMPBLAS2::tpmv<float>;
    mFloatBLAS2Interface.tpsv = &OpenMPBLAS2::tpsv<float>;

    /** double */
    mDoubleBLAS2Interface.gemv = &OpenMPBLAS2::gemv<double>;
    mDoubleBLAS2Interface.symv = &OpenMPBLAS2::symv<double>;
    mDoubleBLAS2Interface.trmv = &OpenMPBLAS2::trmv<double>;
    mDoubleBLAS2Interface.trsv = &OpenMPBLAS2::trsv<double>;
    mDoubleBLAS2Interface.gbmv = &OpenMPBLAS2::gbmv<double>;
    mDoubleBLAS2Interface.sbmv = &OpenMPBLAS2::sbmv<double>;
    mDoubleBLAS2Interface.tbmv = &OpenMPBLAS2::tbmv<double>;
    mDoubleBLAS2Interface.tbsv = &OpenMPBLAS2::tbsv<double>;
    mDoubleBLAS2Interface.ger = &OpenMPBLAS2::ger<double>;
    mDoubleBLAS2Interface.syr = &OpenMPBLAS2::syr<double>;
    mDoubleBLAS2Interface.syr2 = &OpenMPBLAS2::syr2<double>;
    mDoubleBLAS2Interface.spmv = &OpenMPBLAS2::spmv<double>;
    mDoubleBLAS2Interface.spr = &OpenMPBLAS2::spr<double>;
    mDoubleBLAS2Interface.spr2 = &OpenMPBLAS2::spr2<double>;
    mDoubleBLAS2Interface.tpmv = &OpenMPBLAS2::tpmv<double>;
    mDoubleBLAS2Interface.tpsv = &OpenMPBLAS2::tpsv<double>;

    /** Level 3 */

    /** float */
    mFloatBLAS3Interface.gemm = &OpenMPBLAS3::gemm<float>;
    mFloatBLAS3Interface.symm = &OpenMPBLAS3::symm<float>;
    mFloatBLAS3Interface.trmm = &OpenMPBLAS3::trmm<float>;
    mFloatBLAS3Interface.trsm = &OpenMPBLAS3::trsm<float>;
    mFloatBLAS3Interface.syrk = &OpenMPBLAS3::syrk<float>;
    mFloatBLAS3Interface.syrk2 = &OpenMPBLAS3::syrk2<float>;

    /** double */
    mDoubleBLAS3Interface.gemm = &OpenMPBLAS3::gemm<double>;
    mDoubleBLAS3Interface.symm = &OpenMPBLAS3::symm<double>;
    mDoubleBLAS3Interface.trmm = &OpenMPBLAS3::trmm<double>;
    mDoubleBLAS3Interface.trsm = &OpenMPBLAS3::trsm<double>;
    mDoubleBLAS3Interface.syrk = &OpenMPBLAS3::syrk<double>;
    mDoubleBLAS3Interface.syrk2 = &OpenMPBLAS3::syrk2<double>;

    /** LAPACK */

    /** float */
    mFloatLAPACKInterface.getrf = &OpenMPLAPACK::getrf<float>;
    mFloatLAPACKInterface.getri = &OpenMPLAPACK::getri<float>;
    mFloatLAPACKInterface.getinv = &OpenMPLAPACK::getinv<float>;
    mFloatLAPACKInterface.trtrs = &OpenMPLAPACK::trtrs<float>;
    mFloatLAPACKInterface.tptrs = &OpenMPLAPACK::tptrs<float>;
    mFloatLAPACKInterface.laswp = &OpenMPLAPACK::laswp<float>;

    /** double */
    mDoubleLAPACKInterface.getrf = &OpenMPLAPACK::getrf<double>;
    mDoubleLAPACKInterface.getri = &OpenMPLAPACK::getri<double>;
    mDoubleLAPACKInterface.getinv = &OpenMPLAPACK::getinv<double>;
    mDoubleLAPACKInterface.trtrs = &OpenMPLAPACK::trtrs<double>;
    mDoubleLAPACKInterface.tptrs = &OpenMPLAPACK::tptrs<double>;
    mDoubleLAPACKInterface.laswp = &OpenMPLAPACK::laswp<double>;

#ifdef LAMA_MKL_SCALAPACK

    /** float */

    mFloatSCALAPACKInterface.pgetrf = &OpenMPSCALAPACK::pgetrf<float>;
    mFloatSCALAPACKInterface.pgetri = &OpenMPSCALAPACK::pgetri<float>;
    mFloatSCALAPACKInterface.inverse = &OpenMPSCALAPACK::inverse<float>;

    /** double */

    mDoubleSCALAPACKInterface.pgetrf = &OpenMPSCALAPACK::pgetrf<double>;
    mDoubleSCALAPACKInterface.pgetri = &OpenMPSCALAPACK::pgetri<double>;
    mDoubleSCALAPACKInterface.inverse = &OpenMPSCALAPACK::inverse<double>;

#else

    // nothing to do here, as we just do not register routines

#endif

    /** utils */

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
