/**
 * @file interface.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Macros used for LAMA interface
 * @author Thomas Brandes
 * @date 02.04.2013
 * @since 1.0.0
 */
#ifndef LAMA_INTERFACE_MACROS_HPP_
#define LAMA_INTERFACE_MACROS_HPP_

#include <lama/Scalar.hpp>

/** Macro that creates for a function pointer type a table that contains pointer for each type.
 */

#define LAMA_INTERFACE_DEFINE( structname, functionname )                                     \
                                                                                              \
    /** Method that returns the function pointer for the correspondingly typed version. */    \
                                                                                              \
    structname::functionname                                                                  \
    functionname () const                                                                     \
    {                                                                                         \
        return functionname##_Table;                                                          \
    }                                                                                         \
                                                                                              \
    /** Register a function pointer in the type table. */                                     \
                                                                                              \
    void functionname##_add( structname::functionname functionPtr )                           \
    {                                                                                         \
        functionname##_Table = functionPtr;                                                   \
    }                                                                                         \
                                                                                              \
    structname::functionname functionname##_Table;

#define LAMA_INTERFACE_DEFINE_T( structname, functionname )                                   \
                                                                                              \
    /** Method that returns the function pointer for the correspondingly typed version. */    \
                                                                                              \
    template<typename T>                                                                      \
    typename structname<T>::functionname                                                      \
    functionname () const                                                                     \
    {                                                                                         \
        return ( typename structname<T>::functionname )                                       \
               functionname##_Table[ Scalar::getType<T>() ];                                  \
    }                                                                                         \
                                                                                              \
    /** Register a function pointer in the type table. */                                     \
                                                                                              \
    template<typename T>                                                                      \
    void functionname##_add( typename structname<T>::functionname functionPtr )               \
    {                                                                                         \
        functionname##_Table[ Scalar::getType<T>() ]                                          \
         = ( void (*) () ) functionPtr;                                                       \
    }                                                                                         \
                                                                                              \
    void ( *functionname##_Table[ Scalar::UNKNOWN ] ) ();

#define LAMA_INTERFACE_DEFINE_TT( structname, functionname )                                  \
                                                                                              \
    /** Method that returns the function pointer for the correspondingly typed version. */    \
                                                                                              \
    template<typename T1, typename T2>                                                        \
    typename structname<T1, T2>::functionname                                                 \
    functionname () const                                                                     \
    {                                                                                         \
        return ( typename structname<T1,T2>::functionname )                                   \
               functionname##_Table[ Scalar::getType<T1>()][ Scalar::getType<T2>() ];         \
    }                                                                                         \
                                                                                              \
    /** Register a function pointer in the type table. */                                     \
                                                                                              \
    template<typename T1, typename T2>                                                        \
    void functionname##_add( typename structname<T1, T2>::functionname functionPtr )          \
    {                                                                                         \
        functionname##_Table[ Scalar::getType<T1>() ][ Scalar::getType<T2>() ]                \
           = ( void (*) () ) functionPtr;                                                     \
    }                                                                                         \
                                                                                              \
    void ( *functionname##_Table[ Scalar::UNKNOWN ][ Scalar::UNKNOWN] ) ();

/** This macro registers in an interface a function pointer variable 
 *  with a corresponding function.
 */

#define LAMA_INTERFACE_REGISTER( interface, function )                                      \
    interface.function##_add( function );

#define LAMA_INTERFACE_REGISTER_T( interface, function, T )                                 \
    interface.function##_add<T>( function<T> );

#define LAMA_INTERFACE_REGISTER_TT( interface, function, T1, T2 )                           \
    interface.function##_add<T1,T2>( function<T1,T2> );

#endif // LAMA_INTERFACE_MACROS_HPP_
