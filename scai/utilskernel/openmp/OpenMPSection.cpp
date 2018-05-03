/**
 * @file OpenMPSection.cpp
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
 * @brief Implementation of kernel operations for multidimension sections with OpenMP 
 * @author Thomas Brandes
 * @date 15.05.2017
 */

// hpp
#include <scai/utilskernel/openmp/OpenMPSection.hpp>

// local library
#include <scai/utilskernel/SectionKernelTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

namespace scai
{

using common::BinaryOp;
using common::UnaryOp;
using common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPSection::logger, "OpenMP.Section" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPSection::copy2( 
    ValueType targetSection[],
    const IndexType sizes[],
    const IndexType targetDistances[],
    const ValueType sourceSection[],
    const IndexType sourceDistances[] )
{
    SCAI_REGION( "OpenMP.Section.copy2" )

    #pragma omp parallel for
    for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
    {
        for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
        {
            targetSection[i0 * targetDistances[0] + i1 * targetDistances[1]] = 
                sourceSection[i0 * sourceDistances[0] + i1 * sourceDistances[1]];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPSection::assign( 
    ValueType targetSection[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType targetDistances[],
    const ValueType sourceSection[],
    const IndexType sourceDistances[],
    const BinaryOp op,
    const bool swapOperands )
{
    if ( nDims == 0 )
    {
        SCAI_REGION( "OpenMP.Section.assign0" )

        if ( swapOperands )
        {
            *targetSection = applyBinary( *sourceSection, op, *targetSection );
        }
        else
        {
            *targetSection = applyBinary( *targetSection, op, *sourceSection );
        }
    }
    else if ( nDims == 1 )
    {
        SCAI_REGION( "OpenMP.Section.assign1" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            ValueType& target = targetSection[i0 * targetDistances[0]];

            if ( swapOperands )
            {
                target = applyBinary( sourceSection[i0 * sourceDistances[0]], op, target );
            }
            else
            {
                target = applyBinary( target, op, sourceSection[i0 * sourceDistances[0]] );
            }
        }
    }
    else if ( nDims == 2 )
    {
        if ( op == BinaryOp::COPY && !swapOperands )
        {
            copy2( targetSection, sizes, targetDistances, sourceSection, sourceDistances );
            return;
        }

        SCAI_REGION( "OpenMP.Section.assign2" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                ValueType& target = targetSection[i0 * targetDistances[0] + i1 * targetDistances[1]];
                const ValueType& source = sourceSection[i0 * sourceDistances[0] + i1 * sourceDistances[1]];

                if ( swapOperands )
                {
                    target = applyBinary( source, op, target );
                }
                else
                {
                    target = applyBinary( target, op, source );
                }
            }
        }
    }
    else if ( nDims == 3 )
    {
        SCAI_REGION( "OpenMP.Section.assign3" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizes[2]; ++i2 )
                {
                    ValueType& target = targetSection[i0 * targetDistances[0] + i1 * targetDistances[1]
                                                      + i2 * targetDistances[2]];
                    const ValueType& source = sourceSection[i0 * sourceDistances[0] + i1 * sourceDistances[1]
                                                      + i2 * sourceDistances[2]];

                    if ( swapOperands )
                    {
                        target = applyBinary( source, op, target );
                    }
                    else
                    {
                        target = applyBinary( target, op, source );
                    }
                }
            }
        }
    }
    else if ( nDims == 4 )
    {
        SCAI_REGION( "OpenMP.Section.assign4" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizes[2]; ++i2 )
                {
                    for ( IndexType i3 = 0; i3 < sizes[3]; ++i3 )
                    {
                        ValueType& target = targetSection[i0 * targetDistances[0] + i1 * targetDistances[1]
                                                        + i2 * targetDistances[2] + i3 * targetDistances[3]];
                        const ValueType& source = sourceSection[i0 * sourceDistances[0] + i1 * sourceDistances[1]
                                                        + i2 * sourceDistances[2] + i3 * sourceDistances[3]];

                        if ( swapOperands )
                        {
                            target = applyBinary( source, op, target );
                        }
                        else
                        {
                            target = applyBinary( target, op, source );
                        }
                    }
                }
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "dim = " << nDims << " unsupported number of dims for sections" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void OpenMPSection::unaryOp( 
    TargetValueType targetSection[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType targetDistances[],
    const SourceValueType sourceSection[],
    const IndexType sourceDistances[],
    const common::UnaryOp op )
{
    SCAI_LOG_INFO( logger, "unaryOp<" << common::TypeTraits<TargetValueType>::id() 
                            << ", " << common::TypeTraits<SourceValueType>::id() 
                            << ">, #dims = " << nDims << ", op = " << op )
    if ( nDims == 0 )
    {
        SCAI_REGION( "OpenMP.Section.UnaryOp0" )

        *targetSection = static_cast<TargetValueType>( applyUnary( op, *sourceSection ) );
    }
    else if ( nDims == 1 )
    {
        SCAI_REGION( "OpenMP.Section.UnaryOp1" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            TargetValueType& target = targetSection[i0 * targetDistances[0]];
            const SourceValueType& source = sourceSection[i0 * sourceDistances[0]];

            target = static_cast<TargetValueType>( applyUnary( op, source ) );
        }
    }
    else if ( nDims == 2 )
    {
        SCAI_REGION( "OpenMP.Section.UnaryOp2" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                TargetValueType& target = targetSection[i0 * targetDistances[0] + i1 * targetDistances[1]];
                const SourceValueType& source = sourceSection[i0 * sourceDistances[0] + i1 * sourceDistances[1]];

                target = static_cast<TargetValueType>( applyUnary( op, source ) );
            }
        }
    }
    else if ( nDims == 3 )
    {
        SCAI_REGION( "OpenMP.Section.UnaryOp3" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizes[2]; ++i2 )
                {
                    TargetValueType& target = targetSection[i0 * targetDistances[0] + i1 * targetDistances[1]
                                                      + i2 * targetDistances[2]];
                    const SourceValueType& source = sourceSection[i0 * sourceDistances[0] + i1 * sourceDistances[1]
                                                      + i2 * sourceDistances[2]];

                    target = static_cast<TargetValueType>( applyUnary( op, source ) );
                }
            }
        }
    }
    else if ( nDims == 4 )
    {
        SCAI_REGION( "OpenMP.Section.UnaryOp4" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizes[2]; ++i2 )
                {
                    for ( IndexType i3 = 0; i3 < sizes[3]; ++i3 )
                    {
                        TargetValueType& target = targetSection[i0 * targetDistances[0] + i1 * targetDistances[1]
                                                        + i2 * targetDistances[2] + i3 * targetDistances[3]];
                        const SourceValueType& source = sourceSection[i0 * sourceDistances[0] + i1 * sourceDistances[1]
                                                        + i2 * sourceDistances[2] + i3 * sourceDistances[3]];

                        target = static_cast<TargetValueType>( applyUnary( op, source ) );
                    }
                }
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "dim = " << nDims << " unsupported number of dims for sections" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPSection::assignScalar( 
    ValueType section[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType distances[],
    ValueType val,
    const BinaryOp op,
    const bool swapOperands )
{
    if ( nDims == 0 )
    {
        SCAI_REGION( "OpenMP.Section.assignScalar0" )

        *section = swapOperands ? applyBinary( val, op, *section ) : applyBinary( *section, op, val );
    }
    else if ( nDims == 1 )
    {
        SCAI_REGION( "OpenMP.Section.assignScalar1" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            ValueType& target = section[i0 * distances[0]];
            target = swapOperands ? applyBinary( val, op, target ) : applyBinary( target, op, val );
        }
    }
    else if ( nDims == 2 )
    {
        SCAI_REGION( "OpenMP.Section.assignScalar2" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                ValueType& target = section[i0 * distances[0] + i1 * distances[1]];
                target = swapOperands ? applyBinary( val, op, target ) : applyBinary( target, op, val );
            }
        }
    }
    else if ( nDims == 3 )
    {
        SCAI_REGION( "OpenMP.Section.assignScalar3" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizes[2]; ++i2 )
                {
                    ValueType& target = section[i0 * distances[0] + i1 * distances[1] + i2 * distances[2]];
                    target = swapOperands ? applyBinary( val, op, target ) : applyBinary( target, op, val );
                }
            }
        }
    }
    else if ( nDims == 4 )
    {
        SCAI_REGION( "OpenMP.Section.assignScalar4" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizes[2]; ++i2 )
                {
                    for ( IndexType i3 = 0; i3 < sizes[3]; ++i3 )
                    {
                        ValueType& target = section[i0 * distances[0] + i1 * distances[1] + i2 * distances[2] + i3 * distances[3]];
                        target = swapOperands ? applyBinary( val, op, target ) : applyBinary( target, op, val );
                    }
                }
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "dim = " << nDims << " unsupported number of dims for sections" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPSection::UnaryOp( 
    ValueType section[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType distances[],
    const common::UnaryOp op )
{
    if ( nDims == 1 )
    {
        SCAI_REGION( "OpenMP.Section.applyUnary1" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            ValueType& target = section[i0 * distances[0]];
            target = applyUnary( op, target );
        }
    }
    else if ( nDims == 2 )
    {
        SCAI_REGION( "OpenMP.Section.applyUnary2" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                ValueType& target = section[i0 * distances[0] + i1 * distances[1]];
                target = applyUnary( op, target );
            }
        }
    }
    else if ( nDims == 3 )
    {
        SCAI_REGION( "OpenMP.Section.applyUnary3" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizes[2]; ++i2 )
                {
                    ValueType& target = section[i0 * distances[0] + i1 * distances[1] + i2 * distances[2]];
                    target = applyUnary( op, target );
                }
            }
        }
    }
    else if ( nDims == 4 )
    {
        SCAI_REGION( "OpenMP.Section.applyUnary4" )

        #pragma omp parallel for
        for ( IndexType i0 = 0; i0 < sizes[0]; ++i0 )
        {
            for ( IndexType i1 = 0; i1 < sizes[1]; ++i1 )
            {
                for ( IndexType i2 = 0; i2 < sizes[2]; ++i2 )
                {
                    for ( IndexType i3 = 0; i3 < sizes[3]; ++i3 )
                    {
                        ValueType& target = section[i0 * distances[0] + i1 * distances[1] + i2 * distances[2]];
                        target = applyUnary( op, target );
                    }
                }
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "dim = " << nDims << " unsupported number of dims for sections" )
    }
}

/* --------------------------------------------------------------------------- */
/*    struct ArrayKernels methods                                              */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPSection::ArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register SectionKernel OpenMP routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<SectionKernelTrait::assign<ValueType> >( assign, ctx, flag );
    KernelRegistry::set<SectionKernelTrait::assignScalar<ValueType> >( assignScalar, ctx, flag );
    KernelRegistry::set<SectionKernelTrait::UnaryOp<ValueType> >( UnaryOp, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPSection::BinOpKernels<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;

    SCAI_LOG_DEBUG( logger, "register SectionKernel OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<SectionKernelTrait::unaryOp<ValueType, OtherValueType> >( unaryOp, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPSection::OpenMPSection()
{
    SCAI_LOG_INFO( logger, "register section OpenMP kernels for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    kregistry::mepr::RegistratorV<ArrayKernels, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<BinOpKernels, SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPSection::~OpenMPSection()
{
    SCAI_LOG_INFO( logger, "unregister section OpenMP kernels for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<ArrayKernels, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<BinOpKernels, SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPSection OpenMPSection::guard;

} /* end namespace utilskernel */

} /* end namespace scai */
