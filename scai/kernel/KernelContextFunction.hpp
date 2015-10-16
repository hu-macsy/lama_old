/**
 * @file KernelContextFunction.hpp
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
 * @brief Class that stores a context function from the kernel interface
 * @author Thomas Brandes
 * @date 14.10.2015
 */
#pragma once

#include <scai/kernel/ContextFunction.hpp>
#include <scai/kernel/KernelInterface.hpp>

namespace scai
{

namespace interface
{

template<typename FunctionType> 
class KernelContextFunction : public ContextFunction<FunctionType>
{
public:

    /** Constructor by name, function type is given by the template name */

    KernelContextFunction( const char* name ) :
 
        ContextFunction<FunctionType>(),
        mName ( name ) 

    {
        // get this context function pointers via the kernel interface

        KernelInterface::get( *this, mName );
    }

    FunctionType operator[] ( ContextType ctx )
    {
        FunctionType fn = ContextFunction<FunctionType>::get( ctx );

        if ( fn == NULL )
        {
            // Throw exception

            COMMON_THROWEXCEPTION( "Context function " << mName << " not available for context = " << ctx )
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
 * @tname <KernelTrait>  struct that constains signature and name of the context function.
 *
 * \begincode
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

    typedef typename KernelTrait::FuncType ContextFunctionType;

    KernelTraitContextFunction() : 
 
        KernelContextFunction<ContextFunctionType>( KernelTrait::getId() )
    {
    }

    using KernelContextFunction<ContextFunctionType>::operator[];
};

} /* end namespace interface */

} /* end namespace scai */
