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

#define LAMA_INTERFACE_DEFINE( structname, functionname )                                             \
                                                                                                      \
    /** Method that returns the function pointer for the correspondingly typed version. */            \
                                                                                                      \
    structname::functionname                                                                          \
    functionname () const                                                                             \
    {                                                                                                 \
        return functionname##_Table;                                                                  \
    }                                                                                                 \
                                                                                                      \
    /** Register a function pointer in the type table. */                                             \
                                                                                                      \
    void functionname##_add( structname::functionname functionPtr,                                    \
                             bool replace )                                                           \
    {                                                                                                 \
        if ( replace || !functionname##_Table )                                                       \
        {                                                                                             \
            functionname##_Table = functionPtr;                                                       \
        }                                                                                             \
    }                                                                                                 \
                                                                                                      \
    structname::functionname functionname##_Table;

#define LAMA_INTERFACE_DEFINE_T( structname, functionname )                                           \
                                                                                                      \
    /** Method that returns the function pointer for the correspondingly typed version. */            \
                                                                                                      \
    template<typename ValueType>                                                                      \
    typename structname<ValueType>::functionname                                                      \
    functionname () const                                                                             \
    {                                                                                                 \
        return ( typename structname<ValueType>::functionname )                                       \
               functionname##_Table[ Scalar::getType<ValueType>() ];                                  \
    }                                                                                                 \
                                                                                                      \
    /** Register a function pointer in the type table. */                                             \
                                                                                                      \
    template<typename ValueType>                                                                      \
    void functionname##_add( typename structname<ValueType>::functionname functionPtr,                \
                             bool replace )                                                           \
    {                                                                                                 \
        if ( replace || !functionname##_Table[ Scalar::getType<ValueType>() ] )                       \
        {                                                                                             \
            functionname##_Table[ Scalar::getType<ValueType>() ]                                      \
             = ( void (*) () ) functionPtr;                                                           \
        }                                                                                             \
    }                                                                                                 \
                                                                                                      \
    void ( *functionname##_Table[ Scalar::UNKNOWN ] ) ();

#define LAMA_INTERFACE_DEFINE_TT( structname, functionname )                                          \
                                                                                                      \
    /** Method that returns the function pointer for the correspondingly typed version. */            \
                                                                                                      \
    template<typename ValueType1, typename ValueType2>                                                \
    typename structname<ValueType1, ValueType2>::functionname                                         \
    functionname () const                                                                             \
    {                                                                                                 \
        return ( typename structname<ValueType1,ValueType2>::functionname )                           \
               functionname##_Table[ Scalar::getType<ValueType1>()][ Scalar::getType<ValueType2>() ]; \
    }                                                                                                 \
                                                                                                      \
    /** Register a function pointer in the type table. */                                             \
                                                                                                      \
    template<typename ValueType1, typename ValueType2>                                                \
    void functionname##_add( typename structname<ValueType1, ValueType2>::functionname functionPtr,   \
                             bool replace )                                                           \
    {                                                                                                 \
        if ( replace || !functionname##_Table[ Scalar::getType<ValueType1>() ]                        \
                                             [ Scalar::getType<ValueType2>() ] )                      \
        {                                                                                             \
            functionname##_Table[ Scalar::getType<ValueType1>() ][ Scalar::getType<ValueType2>() ]    \
               = ( void (*) () ) functionPtr;                                                         \
        }                                                                                             \
    }                                                                                                 \
                                                                                                      \
    void ( *functionname##_Table[ Scalar::UNKNOWN ][ Scalar::UNKNOWN] ) ();

/** This macro registers in an interface a function pointer variable 
 *  with a corresponding function.
 *
 *  REGISTER : will not overwrite existing entries
 *  REGISTER1 : overwrites existing entries
 */

#define LAMA_INTERFACE_REGISTER( interface, function )                                                \
    interface.function##_add( function,false );

#define LAMA_INTERFACE_REGISTER1( interface, function )                                               \
    interface.function##_add( function, true );

#define LAMA_INTERFACE_REGISTER_T( interface, function, ValueType )                                   \
    interface.function##_add<ValueType>( function<ValueType>, false );

#define LAMA_INTERFACE_REGISTER1_T( interface, function, ValueType )                                  \
    interface.function##_add<ValueType>( function<ValueType>, true );

#define LAMA_INTERFACE_REGISTER_TT( interface, function, ValueType1, ValueType2 )                     \
    interface.function##_add<ValueType1,ValueType2>( function<ValueType1,ValueType2>, false );

#define LAMA_INTERFACE_REGISTER1_TT( interface, function, ValueType1, ValueType2 )                    \
    interface.function##_add<ValueType1,ValueType2>( function<ValueType1,ValueType2>, true );

#endif // LAMA_INTERFACE_MACROS_HPP_
