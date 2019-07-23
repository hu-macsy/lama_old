/**
 * @file ContextStack.hpp
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
 * @brief Class for stack of context used for scoped execution context
 * @author Thomas Brandes
 * @date 25.04.2019
 */
#pragma once

#include <scai/hmemo/Context.hpp>
#include <stack>

namespace scai
{

namespace hmemo
{

/**
 *  @brief Singleton class for the stack used for scope of context.
 */
class ContextStack 
{
public:

    static void push( ContextPtr comm )
    {
        instance.push( comm );
    }

    static void pop()
    {
        instance.pop();
    }

    static bool empty()
    {
        return instance.empty();
    }

    static ContextPtr top()
    {
        return instance.top();
    }

private:
 
    ContextStack();

    static std::stack<ContextPtr> instance;
};

/* -------------------------------------------------------------------------- */

/** 
 *   @brief Help class to guarantee that a pushed context is removed from the stack
 *          at the end of a scope.
 */
class ScopedContextRecord 
{
public:

    /** Constructor to push a context */

    ScopedContextRecord( ContextPtr ctx )
    {
        ContextStack::push( ctx );
    }

    /** Destructor pops a context */

    ~ScopedContextRecord()
    {
        ContextStack::pop();
    }

    // disable all other default constructors

    ScopedContextRecord() = delete;

    ScopedContextRecord( const ScopedContextRecord& ) = delete;

    // disable all default assginments

    ScopedContextRecord& operator= ( const ScopedContextRecord& ) = delete;
};

/* -------------------------------------------------------------------------- */

}   // namespace hmemo

}   // namespace scai

/* -------------------------------------------------------------------------- */

/** Macro that defines a new current ctx for the actual scope
 *
 *  \code
 *      auto ctx = Context::getContextPtr( common::ContextType::CUDA, device );
 *      {
 *          SCAI_HMEMO_CONTEXT( ctx )
 *          ...                        // default context changed
 *      }
 *  \endcode
 */
#define SCAI_HMEMO_CONTEXT( ctx ) scai::hmemo::ScopedContextRecord SCAI_CtxRec_( ctx );

