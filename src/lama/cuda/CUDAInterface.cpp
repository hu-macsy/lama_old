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

LAMA_LOG_DEF_LOGGER( CUDAInterface::logger, "LAMAInterface.CUDAInterface" );

CUDAInterface::CUDAInterface()
    : LAMAInterface()
{
    /** BLAS */

    /** Level 1 */

    /** float */
    mFloatBLAS1Interface.scal = &CUDABLAS1::scal<float>;
    mFloatBLAS1Interface.nrm2 = &CUDABLAS1::nrm2<float>;
    mFloatBLAS1Interface.asum = &CUDABLAS1::asum<float>;
    mFloatBLAS1Interface.iamax = &CUDABLAS1::iamax<float>;
    mFloatBLAS1Interface.viamax = &CUDABLAS1::viamax<float>;
    mFloatBLAS1Interface.swap = &CUDABLAS1::swap<float>;
    mFloatBLAS1Interface.copy = &CUDABLAS1::copy<float>;
    mFloatBLAS1Interface.axpy = &CUDABLAS1::axpy<float>;
    mFloatBLAS1Interface.dot = &CUDABLAS1::dot<float>;
    mFloatBLAS1Interface.sum = &CUDABLAS1::sum<float>;
    mFloatBLAS1Interface.rot = &CUDABLAS1::rot<float>;
    mFloatBLAS1Interface.rotm = &CUDABLAS1::rotm<float>;
//    mFloatBLAS1Interface.rotg    = &CUDABLAS1::rotg<float>;
//    mFloatBLAS1Interface.rotmg   = &CUDABLAS1::rotmg<float>;
    mFloatBLAS1Interface.ass = &CUDABLAS1::ass<float>;

    /** double */
    mDoubleBLAS1Interface.scal = &CUDABLAS1::scal<double>;
    mDoubleBLAS1Interface.nrm2 = &CUDABLAS1::nrm2<double>;
    mDoubleBLAS1Interface.asum = &CUDABLAS1::asum<double>;
    mDoubleBLAS1Interface.iamax = &CUDABLAS1::iamax<double>;
    mDoubleBLAS1Interface.viamax = &CUDABLAS1::viamax<double>;
    mDoubleBLAS1Interface.swap = &CUDABLAS1::swap<double>;
    mDoubleBLAS1Interface.copy = &CUDABLAS1::copy<double>;
    mDoubleBLAS1Interface.axpy = &CUDABLAS1::axpy<double>;
    mDoubleBLAS1Interface.dot = &CUDABLAS1::dot<double>;
    mDoubleBLAS1Interface.sum = &CUDABLAS1::sum<double>;
    mDoubleBLAS1Interface.rot = &CUDABLAS1::rot<double>;
    mDoubleBLAS1Interface.rotm = &CUDABLAS1::rotm<double>;
//    mDoubleBLAS1Interface.rotg    = &CUDABLAS1::rotg<double>;
//    mDoubleBLAS1Interface.rotmg   = &CUDABLAS1::rotmg<double>;
    mDoubleBLAS1Interface.ass = &CUDABLAS1::ass<double>;

    /** Level 2 */

    /** float */
    mFloatBLAS2Interface.gemv = &CUDABLAS2::gemv<float>;
    mFloatBLAS2Interface.symv = &CUDABLAS2::symv<float>;
    mFloatBLAS2Interface.trmv = &CUDABLAS2::trmv<float>;
    mFloatBLAS2Interface.trsv = &CUDABLAS2::trsv<float>;
    mFloatBLAS2Interface.gbmv = &CUDABLAS2::gbmv<float>;
    mFloatBLAS2Interface.sbmv = &CUDABLAS2::sbmv<float>;
    mFloatBLAS2Interface.tbmv = &CUDABLAS2::tbmv<float>;
    mFloatBLAS2Interface.tbsv = &CUDABLAS2::tbsv<float>;
    mFloatBLAS2Interface.ger = &CUDABLAS2::ger<float>;
    mFloatBLAS2Interface.syr = &CUDABLAS2::syr<float>;
    mFloatBLAS2Interface.syr2 = &CUDABLAS2::syr2<float>;
    mFloatBLAS2Interface.spmv = &CUDABLAS2::spmv<float>;
    mFloatBLAS2Interface.spr = &CUDABLAS2::spr<float>;
    mFloatBLAS2Interface.spr2 = &CUDABLAS2::spr2<float>;
    mFloatBLAS2Interface.tpmv = &CUDABLAS2::tpmv<float>;
    mFloatBLAS2Interface.tpsv = &CUDABLAS2::tpsv<float>;

    /** double */
    mDoubleBLAS2Interface.gemv = &CUDABLAS2::gemv<double>;
    mDoubleBLAS2Interface.symv = &CUDABLAS2::symv<double>;
    mDoubleBLAS2Interface.trmv = &CUDABLAS2::trmv<double>;
    mDoubleBLAS2Interface.trsv = &CUDABLAS2::trsv<double>;
    mDoubleBLAS2Interface.gbmv = &CUDABLAS2::gbmv<double>;
    mDoubleBLAS2Interface.sbmv = &CUDABLAS2::sbmv<double>;
    mDoubleBLAS2Interface.tbmv = &CUDABLAS2::tbmv<double>;
    mDoubleBLAS2Interface.tbsv = &CUDABLAS2::tbsv<double>;
    mDoubleBLAS2Interface.ger = &CUDABLAS2::ger<double>;
    mDoubleBLAS2Interface.syr = &CUDABLAS2::syr<double>;
    mDoubleBLAS2Interface.syr2 = &CUDABLAS2::syr2<double>;
    mDoubleBLAS2Interface.spmv = &CUDABLAS2::spmv<double>;
    mDoubleBLAS2Interface.spr = &CUDABLAS2::spr<double>;
    mDoubleBLAS2Interface.spr2 = &CUDABLAS2::spr2<double>;
    mDoubleBLAS2Interface.tpmv = &CUDABLAS2::tpmv<double>;
    mDoubleBLAS2Interface.tpsv = &CUDABLAS2::tpsv<double>;

    /** Level 3 */

    /** float */
    mFloatBLAS3Interface.gemm = &CUDABLAS3::gemm<float>;
//    mFloatBLAS3Interface.symm  = &CUDABLAS3::symm<float>;
//    mFloatBLAS3Interface.trmm  = &CUDABLAS3::trmm<float>;
    mFloatBLAS3Interface.trsm = &CUDABLAS3::trsm<float>;
//    mFloatBLAS3Interface.syrk  = &CUDABLAS3::syrk<float>;
//    mFloatBLAS3Interface.syrk2 = &CUDABLAS3::syrk2<float>;

    /** double */
    mDoubleBLAS3Interface.gemm = &CUDABLAS3::gemm<double>;
//    mDoubleBLAS3Interface.symm  = &CUDABLAS3::symm<double>;
//    mDoubleBLAS3Interface.trmm  = &CUDABLAS3::trmm<double>;
    mDoubleBLAS3Interface.trsm = &CUDABLAS3::trsm<double>;
//    mDoubleBLAS3Interface.syrk  = &CUDABLAS3::syrk<double>;
//    mDoubleBLAS3Interface.syrk2 = &CUDABLAS3::syrk2<double>;

    /** LAPACK */

    /** float */
//    mFloatLAPACKInterface.lamch = &CUDALAPACK::lamch<float>;
//    mFloatLAPACKInterface.getrf = &CUDALAPACK::getrf<float>;
//    mFloatLAPACKInterface.getri = &CUDALAPACK::getri<float>;
//    mFloatLAPACKInterface.trtrs = &CUDALAPACK::trtrs<float>;
//    mFloatLAPACKInterface.tptrs = &CUDALAPACK::tptrs<float>;
    mFloatLAPACKInterface.laswp = &CUDALAPACK::laswp<float>;

    /** double */
//    mDoubleLAPACKInterface.lamch = &CUDALAPACK::lamch<double>;
//    mDoubleLAPACKInterface.getrf = &CUDALAPACK::getrf<double>;
//    mDoubleLAPACKInterface.getri = &CUDALAPACK::getri<double>;
//    mDoubleLAPACKInterface.trtrs = &CUDALAPACK::trtrs<double>;
//    mDoubleLAPACKInterface.tptrs = &CUDALAPACK::tptrs<double>;
    mDoubleLAPACKInterface.laswp = &CUDALAPACK::laswp<double>;

//    /** SCALAPACK */
//
//    /** float */
//
//    mFloatSCALAPACKInterface.pgetrf = &CUDASCALAPACK::pgetrf<float>;
//    mFloatSCALAPACKInterface.pgetri = &CUDASCALAPACK::pgetri<float>;
//
//    /** double */
//
//    mDoubleSCALAPACKInterface.pgetrf = &CUDASCALAPACK::pgetrf<double>;
//    mDoubleSCALAPACKInterface.pgetri = &CUDASCALAPACK::pgetri<double>;

// Each interface part provides an own implementation to
// add function pointers in the tables.

    CUDAUtils::setInterface( Utils );

    CUDACSRUtils::setInterface( CSRUtils );
    CUDAELLUtils::setInterface( ELLUtils );
    CUDAJDSUtils::setInterface( JDSUtils );
    CUDADIAUtils::setInterface( DIAUtils );
    CUDACOOUtils::setInterface( COOUtils );
}

CUDAInterface::~CUDAInterface()
{
}

// static lama::LAMAInterfaceRegistration<CUDAInterface> cudaRegisterObj( Context::CUDA );

LAMA_LAMAINTERFACE_REGISTRATION( Context::CUDA, CUDAInterface )

} // namespace

