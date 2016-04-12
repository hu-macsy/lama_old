/**
 * @file LAMAKernel.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
    it also combines it with ContexPtr instead of context::ContextType
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

    ContextFunctionType operator[] ( hmemo::ContextPtr context )
    {
        return kregistry::KernelTraitContextFunction<KernelTrait>::operator[]( context->getType() );
    }

    hmemo::ContextPtr getValidContext( hmemo::ContextPtr defaultContext ) 
    {
        SCAI_ASSERT_DEBUG( defaultContext.get(), "NULL context" );

        common::context::ContextType defCtx = defaultContext->getType();
        common::context::ContextType runCtx = kregistry::_ContextFunction::validContext( defCtx );

        if ( runCtx == defCtx )
        {
            return defaultContext;
        }
        else
        {
            return hmemo::Context::getHostPtr();  // do it on host
        }
    }

    /** Change context to host if function is not supported on preferred context. */

    void getSupportedContext( hmemo::ContextPtr& context ) const
    {
        SCAI_ASSERT_DEBUG( context.get(), "NULL context" );

        common::context::ContextType defCtx = context->getType();
        common::context::ContextType runCtx = kregistry::_ContextFunction::validContext( defCtx );

        if ( runCtx != defCtx )
        {
            context = hmemo::Context::getHostPtr();
        }
    }

    void getSupportedContext( hmemo::ContextPtr& context, kregistry::_ContextFunction other ) const
    {
        SCAI_ASSERT_DEBUG( context.get(), "NULL context" );

        common::context::ContextType defCtx = context->getType();
        common::context::ContextType runCtx = kregistry::_ContextFunction::validContext( other, defCtx );

        if ( runCtx == defCtx )
        {
            // context is fine, this and other are supported
        }
        else if ( runCtx == common::context::Host )
        {
            context = hmemo::Context::getHostPtr();  // do it on host
        }
        else
        {
            COMMON_THROWEXCEPTION( "illegal registry" )
        }
    }

    hmemo::ContextPtr getValidContext( kregistry::_ContextFunction other1, 
                                       kregistry::_ContextFunction other2,
                                       hmemo::ContextPtr defaultContext    )
    {
        SCAI_ASSERT_DEBUG( defaultContext.get(), "NULL context" );

        common::context::ContextType defCtx = defaultContext->getType();
        common::context::ContextType runCtx = kregistry::_ContextFunction::validContext( other1, other2, defCtx );

        if ( runCtx == defCtx )
        {
            return defaultContext;
        }
        else if ( runCtx == common::context::Host )
        {
            return hmemo::Context::getHostPtr();  // do it on host
        }
        else
        {
            COMMON_THROWEXCEPTION( "illegal registry" )
        }
    }
};

}

}
