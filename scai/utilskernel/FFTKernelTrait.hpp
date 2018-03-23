/**
 * @file FFTKernelTrait.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Struct with traits for abstraction of different FFT-libraries
 * @author Eric Schricker
 * @date 29.09.2016
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/Complex.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

/** Namespace for utilities on heterogeneous arrays (HArray) and derived class LArray */

namespace utilskernel
{

/** Structure with traits for all Utils kernels. */

struct FFTKernelTrait
{
    /** @brief Trait for register kernel function paddedForward1 that applies a forward fft on padded data
     */
    template <typename ValueType>
    struct paddedForward1D
    {
        /** @brief one dimensional fft
         *
         *  @param[in] n is the size of the in-array
         *  @param[in] npad is the size of the out-array
         *  @param[in] in is an array of values
         *  @param[out] out is an array which contains the transformed values
         */

        typedef void ( *FuncType ) ( 
            const IndexType n, 
            const IndexType npad, 
            const ValueType in[], 
            common::Complex<RealType<ValueType>> out[] );

        static const char* getId()
        {
            return "FFTKernel.paddedForward1D";
        }
    };

    /** @brief Trait for register kernel function paddedForward1 that applies a backward fft on padded data
     */
    template <typename ValueType>
    struct paddedBackward1D
    {
        /** @brief one dimensional fft
         *
         *  @param[in] n is the size of the in-array
         *  @param[in] npad is the size of the out-array
         *  @param[in] in is an array of values
         *  @param[out] out is an array which contains the transformed values
         */
        typedef void ( *FuncType ) ( 
            const IndexType n, 
            const IndexType npad, 
            const ValueType in[], 
            common::Complex<RealType<ValueType>> out[] );

        static const char* getId()
        {
            return "FFTKernel.paddedBackward1D";
        }
    };

    template <typename ValueType>
    struct fft
    {
        /** @brief one dimensional fft in-place
         *
         *  @param[in,out] array used for input and output
         *  @param[in] n is the size of the array, must be power of 2
         *  @param[in] m is the log of n so that n == 2**m
         *  @param[in] direction is either 1 (forward) or -1 (backward)
         */
        typedef void ( *FuncType ) ( 
            common::Complex<ValueType> array[],
            const IndexType n,
            const IndexType m,
            const int direction );

        static const char* getId()
        {
            return "FFTKernel.fft";
        }
    };

    template <typename ValueType>
    struct fftK
    {
        /** @brief one dimensional fft in-place for multiple (row) vectors
         *
         *  @param[in,out] array used for input and output, size is k x n
         *  @param[in] k is the number of rows
         *  @param[in] n is the size of the array, must be power of 2
         *  @param[in] m is the log of n so that n == 2**m
         *  @param[in] direction is either 1 (forward) or -1 (backward)
         */
        typedef void ( *FuncType ) ( 
            common::Complex<ValueType> array[],
            const IndexType k,
            const IndexType n,
            const IndexType m,
            const int direction );

        static const char* getId()
        {
            return "FFTKernel.fftK";
        }
    };

};

} /* end namespace utilskernel */

} /* end namespace scai */
