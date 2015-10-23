/**
 * @file kernel_registry.hpp
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
 * @brief Interface class for context dependent operations to be implemented.
 * @author Thomas Brandes
 * @date 27.04.2011
 * @since 1.0.0
 */
#pragma once

#include <scai/kregistry/KernelContextFunction.hpp>

namespace scai
{
namespace lama
{

/** Define class LAMAKernel that deals always with KernelTraits;
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

    ContextFunctionType operator[] ( hmemo::ContextPtr context )
    {
        return kregistry::KernelTraitContextFunction<KernelTrait>::operator[]( context->getType() );
    }

    hmemo::ContextPtr getValidContext( hmemo::ContextPtr defaultContext )
    {
        common::ContextType defCtx = defaultContext->getType();
        common::ContextType runCtx = kregistry::_ContextFunction::validContext( defCtx );

        if ( runCtx == defCtx )
        {
            return defaultContext;
        }
        else
        {
            return hmemo::Context::getHostPtr();  // do it on host
        }
    }

    hmemo::ContextPtr getValidContext( kregistry::_ContextFunction other, hmemo::ContextPtr defaultContext )
    {
        common::ContextType defCtx = defaultContext->getType();
        common::ContextType runCtx = kregistry::_ContextFunction::validContext( other, defCtx );

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

    hmemo::ContextPtr getValidContext( kregistry::_ContextFunction other1, 
                                       kregistry::_ContextFunction other2,
                                       hmemo::ContextPtr defaultContext    )
    {
        common::ContextType defCtx = defaultContext->getType();
        common::ContextType runCtx = kregistry::_ContextFunction::validContext( other1, other2, defCtx );

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

/** Help routine to update context. */

static inline hmemo::ContextPtr getValidContext( hmemo::ContextPtr defaultContext, scai::kregistry::_ContextFunction f )
{
    common::ContextType defCtx = defaultContext->getType();
    common::ContextType runCtx = f.validContext( defCtx );

    if ( runCtx == defCtx )
    {
        return defaultContext;
    }
    else
    {
        return hmemo::Context::getHostPtr();  // do it on host
    }
}

}

}
