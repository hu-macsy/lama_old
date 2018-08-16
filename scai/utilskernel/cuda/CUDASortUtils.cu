/**
 * @file utilskernel/cuda/CUDASortUtils.cu
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
 * @brief Implementation of CSR utilities with CUDA
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/utilskernel/cuda/CUDASortUtils.hpp>
#include <scai/utilskernel/cuda/CUDASparseUtils.hpp>



// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/Math.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/SparseKernelTrait.hpp>

// thrust
#include <thrust/iterator/constant_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>


#include <complex.h>


using namespace scai::common;

namespace scai
{

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( CUDASortUtils::logger, "CUDA.Utils" )

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDASortUtils::sortBoth( ValueType array[], IndexType perm[], const IndexType n, bool ascending )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;

    SCAI_REGION( "CUDA.Utils.sortBoth" )

    thrust::device_ptr<RealType> array_d( reinterpret_cast<RealType*>( array ) );
    thrust::device_ptr<IndexType> perm_d( perm );

    // stable sort, descending order, so override default comparison

    if ( ascending )
    {
        thrust::stable_sort_by_key( array_d, array_d + n, perm_d, thrust::less<RealType>() );
    }
    else
    {
        thrust::stable_sort_by_key( array_d, array_d + n, perm_d, thrust::greater<RealType>() );
    }
}

template<typename ValueType>
void CUDASortUtils::sortValues( ValueType array[], const IndexType n, bool ascending )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;

    SCAI_REGION( "CUDA.Utils.sortValues" )

    thrust::device_ptr<RealType> array_d( reinterpret_cast<RealType*>( array ) );

    if ( ascending )
    {
        thrust::sort( array_d, array_d + n, thrust::less<RealType>() );
    }
    else
    {
        thrust::sort( array_d, array_d + n, thrust::greater<RealType>() );
    }
}

template<typename ValueType>
struct myLess
{
    const ValueType* mArray;

    myLess( const ValueType array[] ) : mArray( array )
    {
    }

    __host__ __device__
    bool operator()( const IndexType& i1, const IndexType& i2 )
    {
        return mArray[i1] < mArray[i2];
    }
};


template<typename ValueType>
struct myGreater
{
    const ValueType* mArray;

    myGreater( const ValueType array[] ) : mArray( array )
    {
    }

    __host__ __device__
    bool operator()( const IndexType& i1, const IndexType& i2 )
    {
        return mArray[i1] > mArray[i2];
    }
};

template<typename ValueType>
void CUDASortUtils::sortPerm( IndexType perm[], const ValueType array[], const IndexType n, bool ascending )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;

    SCAI_REGION( "CUDA.Utils.sortValues" )

    thrust::device_ptr<IndexType> perm_d( perm );

    if ( TypeTraits<ValueType>::stype != TypeTraits<RealType>::stype )
    {
        COMMON_THROWEXCEPTION( "sort unsupported for complex values" )
    }
    else
    {
        const RealType* absArray = reinterpret_cast<const RealType*>( array );

        // stable sort, descending order, so override default comparison

        if ( ascending )
        {
            thrust::stable_sort( perm_d, perm_d + n, myLess<RealType>( absArray ) );
        }
        else
        {
            thrust::stable_sort( perm_d, perm_d + n, myGreater<RealType>( absArray ) );
        }
    }
}

template<typename ValueType>
void CUDASortUtils::sort(
    IndexType perm[],
    ValueType outValues[],
    const ValueType inValues[],
    const IndexType n,
    bool ascending )
{
    SCAI_REGION( "CUDA.Utils.sort" )

    SCAI_LOG_INFO( logger, "sort " << n << " values, ascending = " << ascending )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS

    if ( perm == NULL )
    {
        // just sort the values, but can only be done in place

        if ( inValues != outValues )
        {
            CUDASparseUtils::set( outValues, inValues, n, BinaryOp::COPY );
        }

        if ( n > 1 )
        {
            sortValues( outValues, n, ascending );
        }

        return;
    }

    // Initialize permutation array

    thrust::device_ptr<IndexType> perm_d( perm );
    thrust::sequence( perm_d, perm_d + n );

    if ( outValues == NULL )
    {
        if ( n > 1 )
        {
            sortPerm( perm, inValues, n, ascending );
        }

        return;
    }

    if ( inValues != outValues )
    {
        CUDASparseUtils::set( outValues, inValues, n, common::BinaryOp::COPY );
    }

    if ( n > 1 )
    {
        sortBoth( outValues, perm, n, ascending );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDASortUtils::sortInPlace(
    IndexType indexes[],
    ValueType values[],
    const IndexType n,
    bool ascending )
{
    SCAI_REGION( "CUDA.Utils.sort" )

    SCAI_LOG_INFO( logger, "sort " << n << " values, ascending = " << ascending )

    if ( n <= 1 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<IndexType> array_d( indexes );
    thrust::device_ptr<ValueType> values_d( values );

    if ( ascending )
    {
        thrust::sort_by_key( array_d, array_d + n, values_d, thrust::less<IndexType>() );
    }
    else
    {
        thrust::sort_by_key( array_d, array_d + n, values_d, thrust::greater<IndexType>() );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASortUtils::RegArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;

    SCAI_LOG_DEBUG( logger, "registerV array UtilsKernel CUDA [" << flag
                    << "] --> ValueType = " << common::getScalarType<ValueType>() )

    // Note: these kernels will be instantiated for numeric types + IndexType

    KernelRegistry::set<UtilKernelTrait::sort<ValueType> >( sort, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::sortInPlace<ValueType> >( sortInPlace, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDASortUtils::CUDASortUtils()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDASortUtils::~CUDASortUtils()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDASortUtils CUDASortUtils::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
