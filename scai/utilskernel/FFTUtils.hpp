/**
 * @file utilskernel/FFTUtils.hpp
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
 * @brief Fast fourier transformation for heterorgeneous arrays
 * @author Thomas Brandes
 * @date 21.03.2018
 */
#pragma once

#include <scai/hmemo/HArray.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace utilskernel
{

/** Static class of fft routines on HArray data
 *
 */
class COMMON_DLL_IMPORTEXPORT FFTUtils
{
public:

#ifdef SCAI_COMPLEX_SUPPORTED

/** Compute the discrete fourier transform of an array using the FFT algorithm
 *
 *  @param[out] result is the result array
 *  @param[in]  x  is the input array
 *  @param[in]  n  padding length (n > x.size()) or truncate length , optional
 *  @param[in]  direction must be either 1 (forward) or -1 (backward, inverse)
 *  @param[in]  ctx preferred context for execution   
 */
template<typename ValueType>
static void fft1D( 
    hmemo::HArray<common::Complex<RealType<ValueType>>>& out, 
    const hmemo::HArray<ValueType>& in,
    const IndexType n,
    const hmemo::ContextPtr ctx = hmemo::ContextPtr() );

template<typename ValueType>
static void ifft1D( 
    hmemo::HArray<ValueType>& out, 
    const hmemo::HArray<common::Complex<RealType<ValueType>>>& in, 
    const IndexType n,
    const hmemo::ContextPtr ctx = hmemo::ContextPtr() );

/** Compute the discrete fourier transform multiple vectors using the FFT algorithm
 *
 *  @param[out] result is the result array, size will be many * n
 *  @param[in]  x is the input array, size must be multiple of many
 *  @param[in]  many is the number of rows, valid for input and result array
 *  @param[in]  n is the padding length for each row vector
 *  @param[in]  direction must be either 1 (forward) or -1 (backward, inverse)
 *  @param[in]  ctx preferred context for execution   
 */
template<typename ValueType>
static void fft_many( 
    hmemo::HArray<common::Complex<RealType<ValueType>>>& result, 
    const hmemo::HArray<ValueType>& x, 
    const IndexType many,
    const IndexType n, 
    const int direction,
    const hmemo::ContextPtr ctx = hmemo::ContextPtr() );

template<typename ValueType>
static void fftcall(
    hmemo::HArray<common::Complex<RealType<ValueType>>>& data,
    const IndexType k,
    const IndexType n,
    const IndexType m,
    int direction,
    const hmemo::ContextPtr ctx = hmemo::ContextPtr() );

#endif

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    FFTUtils();  // static class, no objects outside

    FFTUtils( const FFTUtils& ) = delete;

};

} /* end namespace utilskernel */

} /* end namespace scai */

