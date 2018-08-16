/**
 * @file LAMAKernel.hpp
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
 * @brief Derived class of KernelContexFunction for more convenient use in LAMA
 * @author Thomas Brandes
 * @date 20.10.2015
 */
#pragma once

#include <scai/kregistry/KernelContextFunction.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/common/macros/assert.hpp>


namespace scai
{
namespace utilskernel
{

/** Define class LAMAKernel that deals always with KernelTrait;
    it also combines it with ContexPtr instead of ContextType
*/

template<typename KernelTrait>
class LAMAKernel : public kregistry::KernelTraitContextFunction<KernelTrait>
{
public:

    typedef typename KernelTrait::FuncType ContextFunctionType;

    LAMAKernel() : kregistry::KernelTraitContextFunction<KernelTrait>()
    {
    }

    /** more convenient now is to access the routine by context pointer */

    ContextFunctionType operator[] ( const hmemo::ContextPtr& context )
    {
        return kregistry::KernelTraitContextFunction<KernelTrait>::operator[]( context->getType() );
    }

    /** Update context if function is not supported on preferred context.
     *
     *  @param[in,out] context where the kernel should be executed
     */

    void getSupportedContext( hmemo::ContextPtr& context ) const
    {
        SCAI_ASSERT_DEBUG( context.get(), "NULL context" );
        common::ContextType defCtx = context->getType();
        common::ContextType runCtx = kregistry::_ContextFunction::validContext( defCtx );

        if ( runCtx != defCtx )
        {
            context = hmemo::Context::getHostPtr();
        }
    }

    /** Update context if function is not supported on preferred context.
     *
     *  @param[in,out] context where the kernel should be executed (in) and can be executed (out)
     *  @param[in] other another kernel function that is also called
     */

    void getSupportedContext( hmemo::ContextPtr& context, const kregistry::_ContextFunction& other ) const
    {
        SCAI_ASSERT_DEBUG( context.get(), "NULL context" );
        common::ContextType defCtx = context->getType();
        common::ContextType runCtx = kregistry::_ContextFunction::validContext( other, defCtx );

        if ( runCtx == defCtx )
        {
            // context is fine, this and other are supported
        }
        else if ( runCtx == common::ContextType::Host )
        {
            context = hmemo::Context::getHostPtr();  // do it on host
        }
        else
        {
            COMMON_THROWEXCEPTION( "illegal registry" )
        }
    }

    void getSupportedContext( hmemo::ContextPtr& context, const kregistry::_ContextFunction& other1,
                              const kregistry::_ContextFunction& other2 ) const
    {
        SCAI_ASSERT_DEBUG( context.get(), "NULL context" );
        common::ContextType defCtx = context->getType();
        common::ContextType runCtx = kregistry::_ContextFunction::validContext( other1, other2, defCtx );

        if ( runCtx == defCtx )
        {
            // context is fine, this and other are supported
        }
        else if ( runCtx == common::ContextType::Host )
        {
            context = hmemo::Context::getHostPtr();  // do it on host
        }
        else
        {
            COMMON_THROWEXCEPTION( "illegal registry" )
        }
    }
};

} /* end namespace utilskernel */

} /* end namespace scai */
