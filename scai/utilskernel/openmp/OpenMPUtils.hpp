/**
 * @file OpenMPUtils.hpp
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
 * @brief Implementation of general utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/CompareOp.hpp>
#include <scai/common/UnaryOp.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace utilskernel
{

/** General utilities of the LAMA Interface implemented in OpenMP  */

class COMMON_DLL_IMPORTEXPORT OpenMPUtils
{
public:

    /** OpenMP implementation of UtilKernelTrait::validIndexes */

    static bool validIndexes( const IndexType array[], const IndexType n, const IndexType size );

    /** OpenMP implementation for UtilKernelTrait::reduce */

    template<typename ValueType>
    static ValueType reduce(
        const ValueType array[],
        const IndexType n,
        const ValueType zero,
        const common::BinaryOp op );

    /** OpenMP implementation for UtilKernelTrait::reduce2 */

    template<typename ValueType>
    static ValueType reduce2(
        const ValueType array1[],
        const ValueType array2[],
        const IndexType n,
        const common::BinaryOp binOp,
        const ValueType zero,
        const common::BinaryOp redOp );

    /** OpenMP implementation for UtilKernelTrait::allCompare */

    template<typename ValueType>
    static bool allCompare(
        const ValueType array1[],
        const ValueType array2[],
        const IndexType n,
        const common::CompareOp op );

    /** OpenMP implementation for UtilKernelTrait::allCompareScalar */

    template<typename ValueType>
    static bool allCompareScalar(
        const ValueType array[],
        const ValueType scalar,
        const IndexType n,
        const common::CompareOp op );

    /** OpenMP implementation for UtilKernelTrait::setVal */

    template<typename ValueType>
    static void setVal( ValueType array[], const IndexType n, const ValueType val, const common::BinaryOp op );

    /** OpenMP implementation for UtilKernelTrait::scaleVectorAddScalar */

    template<typename ValueType>
    static void scaleVectorAddScalar( ValueType array1[], const ValueType array2[], const IndexType n, const ValueType alpha, const ValueType beta );

    /** OpenMP implementation for UtilKernelTrait::setOrder */

    template<typename ValueType>
    static void setOrder( ValueType array[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::setSequence */

    template<typename ValueType>
    static void setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::getValue */

    template<typename ValueType>
    static ValueType getValue( const ValueType* array, const IndexType i );

    /** OpenMP implementation for UtilKernelTrait::isSorted */

    template<typename ValueType>
    static bool isSorted( const ValueType array[], const IndexType n, const common::CompareOp op );

    /** OpenMP implementation for UtilKernelTrait::set */

    template<typename ValueType1, typename ValueType2>
    static void set( ValueType1 out[], const ValueType2 in[], const IndexType n, const common::BinaryOp op );

    /** OpenMP implementation for UtilKernelTrait::setSection */

    template<typename ValueType1, typename ValueType2>
    static void setSection(
        ValueType1 out[],
        const IndexType inc1,
        const ValueType2 in[],
        const IndexType inc2,
        const IndexType n,
        const common::BinaryOp op );

    /** OpenMP implementation for UtilKernelTrait::fillSection */

    template<typename ValueType>
    static void fillSection(
        ValueType out[],
        const IndexType inc,
        const ValueType val,
        const IndexType n,
        const common::BinaryOp op );

    /** OpenMP implementation for UtilKernelTrait::unaryOp */

    template<typename ValueType>
    static void unaryOp( ValueType out[], const ValueType in[], const IndexType n, const common::UnaryOp op );

    /** OpenMP implementation for UtilKernelTrait::binaryOp */

    template<typename ValueType>
    static void binaryOp( ValueType out[], const ValueType in1[], const ValueType in2[], const IndexType n, const common::BinaryOp op );

    /** OpenMP implementation for UtilKernelTrait::binaryOpScalar */

    template<typename ValueType>
    static void binaryOpScalar( 
        ValueType out[], 
        const ValueType in[], 
        const ValueType value, 
        const IndexType n, 
        const common::BinaryOp op,
        const bool swapScalar );

    /** OpenMP implementation for UtilKernelTrait::setGather */

    template<typename ValueType1, typename ValueType2>
    static void setGather(
        ValueType1 out[],
        const ValueType2 in[],
        const IndexType indexes[],
        const common::BinaryOp op,
        const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::setGatherSparse */

    template<typename ValueType1, typename ValueType2>
    static void setGatherSparse(
        ValueType1 target[],
        const ValueType2 sourceZeroVal,
        const ValueType2 sourceNonZeroValues[],
        const IndexType sourceNonZeroIndexes[],
        const IndexType sourceNNZ,
        const IndexType indexes[],
        const common::BinaryOp op,
        const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::scatterVal */

    template<typename ValueType>
    static void scatterVal( ValueType out[], const IndexType indexes[], const ValueType value, const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::setScatter */

    template<typename ValueType1, typename ValueType2>
    static void setScatter( ValueType1 out[],
                            const IndexType indexes[],
                            const bool unique,
                            const ValueType2 in[],
                            const common::BinaryOp op,
                            const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::scan */

    template<typename ValueType>
    static ValueType scan( 
        ValueType array[], 
        const IndexType n, 
        const ValueType first,
        const bool exclusive,
        const bool append );

    /** OpenMP implementation for UtilKernelTrait::unscan */

    template<typename ValueType>
    static ValueType unscan( ValueType array[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::binarySearch */

    static void binarySearch( IndexType outPos[],
                              const IndexType indexes[], const IndexType m,
                              const IndexType inPos[], const IndexType n );

    /** OpenMP implementation for UtilKernelTrait::sort */

    template<typename ValueType>
    static void sort(
        IndexType perm[],
        ValueType outValues[],
        const ValueType inValues[],
        const IndexType n,
        const bool ascending );

    /** OpenMP implementation for UtilKernelTrait::sortInPlace */

    template<typename ValueType>
    static void sortInPlace(
        IndexType indexes[],
        ValueType values[],
        const IndexType n,
        const bool ascending );

    /** Compute the inverse permutation as specified in UtilKernelTrait::setInversePerm */

    static void setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n );

    /** Count bucket sizes for values mapped to buckets, see UtilKernelTrait::countBuckets */

    template<typename BucketType>
    static void countBuckets( IndexType bucketSizes[], const BucketType nBuckets, const BucketType bucketMap[], const IndexType n );

    /** Resort indexes 0, ..., n-1 according to their mapping to buckets, see UtilKernelTrait::sortInBuckets */

    template<typename BucketType>
    static void sortInBuckets( IndexType sortedIndexes[],
                               IndexType offsets[],
                               const BucketType nBuckets,
                               const BucketType bucketMap[],
                               const IndexType n );

    /** OpenMP implementation of SparseKernelTrait::countAddSparse */

    static IndexType countAddSparse(
        const IndexType indexes1[],
        const IndexType n1,
        const IndexType indexes2[],
        const IndexType n2 );

    /** OpenMP implementation of SparseKernelTrait::addSparse */

    template<typename ValueType>
    static IndexType addSparse(
        IndexType indexes[],
        ValueType values[],
        const IndexType indexes1[],
        const ValueType values1[],
        const ValueType zero1,
        const IndexType n1,
        const ValueType alpha,
        const IndexType indexes2[],
        const ValueType values2[],
        const ValueType zero2,
        const IndexType n2,
        const ValueType beta );

    /** OpenMP implementation of SparseKernelTrait::binopSparse */

    template<typename ValueType>
    static IndexType binopSparse(
        IndexType indexes[],
        ValueType values[],
        const IndexType indexes1[],
        const ValueType values1[],
        const ValueType zero1,
        const IndexType n1,
        const IndexType indexes2[],
        const ValueType values2[],
        const ValueType zero2,
        const IndexType n2,
        const common::BinaryOp op );

    /** OpenMP implementation of SparseKernelTrait::joinSparse */

    template<typename ValueType>
    static IndexType joinSparse(
        IndexType indexes[],
        ValueType values[],
        const IndexType indexes1[],
        const ValueType values1[],
        const IndexType n1,
        const IndexType indexes2[],
        const ValueType values2[],
        const IndexType n2 );

private:

    /** Optimized reduce for common::BinaryOp::ADD as reduction operator. */

    template<typename ValueType>
    static ValueType reduceSum( const ValueType array[], const IndexType n, const ValueType zero );

    template<typename ValueType>
    static ValueType reduceMaxVal( const ValueType array[], const IndexType n, const ValueType zero );

    template<typename ValueType>
    static ValueType reduceMinVal( const ValueType array[], const IndexType n, const ValueType zero );

    template<typename ValueType>
    static ValueType reduceAbsMaxVal( const ValueType array[], const IndexType n, const ValueType zero );

    /** The following method is the same as reduce but will not switch for optimized routines any more. */

    template<typename ValueType>
    static ValueType reduceBinOp(
        const ValueType array[],
        const IndexType n,
        const ValueType zero,
        const common::BinaryOp op );

    template<typename ValueType>
    static ValueType absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n );

    template<typename ValueType>
    static ValueType scanSerial( ValueType array[], const IndexType n, const ValueType first, const bool exclusive );

    template<typename ValueType>
    static ValueType scanParallel( PartitionId numThreads, ValueType array[], const IndexType n, const ValueType first, const bool exclusive );

    /** OpenMP implementation of UtilsKernelTrait::countNonZeros */

    template<typename ValueType>
    static IndexType countNonZeros( const ValueType denseArray[], const IndexType n, const ValueType zero, const ValueType eps );

    /** OpenMP implementation of UtilsKernelTrait::compress */

    template<typename TargetType, typename SourceType>
    static IndexType compress(
        TargetType sparseArray[],
        IndexType sparseIndexes[],
        const SourceType denseArray[],
        const IndexType n,
        const SourceType zero,
        const SourceType eps );

    /** OpenMP implementation of SparseKernelTrait::allCompareSparse */

    template<typename ValueType>
    static IndexType allCompareSparse(
        bool& allFlag,
        const IndexType indexes1[],
        const ValueType values1[],
        const ValueType zero1,
        const IndexType n1,
        const IndexType indexes2[],
        const ValueType values2[],
        const ValueType zero2,
        const IndexType n2,
        const common::CompareOp op );

    /** OpenMP implementation of SparseKernelTrait::mergeSparse */

    template<typename ValueType>
    static IndexType mergeSparse(
        IndexType indexes[],
        ValueType values[],
        const IndexType indexes1[],
        const ValueType values1[],
        const IndexType n1,
        const IndexType indexes2[],
        const ValueType values2[],
        const IndexType n2,
        const common::BinaryOp op );

    template<typename ValueType>
    static void sortValues( ValueType outValues[], const ValueType inValues[], const IndexType n, const bool ascending );

    template<typename KeyType, typename ValueType>
    static void qsort( KeyType keys[], ValueType values[], IndexType left, IndexType right, bool ascending );

    /** Compute the inverse permutation as specified in UtilKernelTrait::setInversePerm */

    /** Routine that registers all methods at the kernel registry. */

    struct BaseKernels
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    template<typename ValueType>
    struct NumericKernels
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    template<typename ValueType>
    struct ArrayKernels
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    template<typename ValueType, typename OtherValueType>
    struct BinOpKernels
    {
        static void registerKernels( const scai::kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    OpenMPUtils();

    /** Destructor for unregistration. */

    ~OpenMPUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPUtils guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace utilskernel */

} /* end namespace scai */
