/**
 * @file KernelContextFunction.hpp
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
 * @brief Class that stores a context function from the kernel registry
 * @author Thomas Brandes
 * @date 14.10.2015
 */
#pragma once

#include <scai/kregistry/ContextFunction.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

namespace scai
{

namespace kregistry
{

template<typename FunctionType>
class COMMON_DLL_IMPORTEXPORT KernelContextFunction : public ContextFunction<FunctionType>
{
public:

    /** Constructor by name, function type is given by the template name */

    KernelContextFunction( const char* name ) :

        ContextFunction<FunctionType>(),
        mName ( name )

    {
        // get this context function pointers via the kernel registry
        KernelRegistry::get( *this, mName );
    }

    FunctionType operator[] ( common::ContextType ctx )
    {
        FunctionType fn = ContextFunction<FunctionType>::get( ctx );

        if ( fn == NULL )
        {
            // Throw exception
            SCAI_THROWEXCEPTION( KernelRegistryException,
                                 "Context function " << mName << " - " << typeid( FunctionType ).name()
                                 << " not available for context = " << ctx
                                 << ", registered is " <<  this->printIt() )
        }

        return fn;
    }

private:

    const char* mName;   // keep the name for error messages

    using _ContextFunction::mContextFuncArray;
};

/**
 * Template class for ContextFunction by using a Kernel Trait
 *
 * @tparam KernelTrait struct that constains signature (function type defintion) and name of the context function.
 *
 * \code
 *     // Example of KernelTrait
 *     struct isSorted
 *     {
 *         bool ( *FuncType ) ( const double* array, int n, bool ascending );
 *         static inline const char* getId() { return "isSorted" };
 *     };
 * \endcode
 */

template<typename KernelTrait>
class KernelTraitContextFunction : public KernelContextFunction<typename KernelTrait::FuncType>
{
public:

    /** Take over the function type of the trait. */

    typedef typename KernelTrait::FuncType ContextFunctionType;

    KernelTraitContextFunction() :

        KernelContextFunction<ContextFunctionType>( KernelTrait::getId() )
    {
    }

    using KernelContextFunction<ContextFunctionType>::operator[];
};

} /* end namespace kregistry */

} /* end namespace scai */
