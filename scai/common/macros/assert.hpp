/**
 * @file assert.hpp
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
 * @brief Definition of macros for assertions.
 * @author eschricker
 * @date 31.08.2015
 */
#pragma once

// local library
// #include <scai/common/macros/throw.hpp>
#include <scai/common/Utils.hpp>
#include <scai/common/exception/AssertException.hpp>

#define UNUSED_EXP( expr )    \
    if ( false )              \
    {                         \
        ( void )( expr );     \
    } 

#define UNUSED_EXP2( expr1, expr2 )    \
    if ( false )                       \
    {                                  \
        ( void )( expr1 );             \
        ( void )( expr2 );             \
    }

#define UNUSED_STRING( ... )                                                    \
    if ( false )                                                                \
    {                                                                           \
        std::cout<< __VA_ARGS__;                                                \
    }

/**
 * @brief The macro SCAI_ASSERT checks a condition and throws an exception
 *        when the condition fails.
 *
 * @param[in] cond  boolean expression that is checked
 * @param[in] msg   message to indicate reason for the exception
 * @throws    common::Exception (derived from std::exception)
 *
 * The message text will also contain the condition as string, file location
 * and the current call stack.
 */

#define SCAI_ASSERT( cond, msg )                                               \
    {                                                                              \
        if (!(cond))                                                               \
        {                                                                          \
            std::ostringstream errorStr;                                           \
            errorStr << "Assertion failed in line " << __LINE__;                   \
            errorStr << " of file " << __FILE__ << "\n";                           \
            errorStr << "    Condition: " << #cond << "\n";                        \
            errorStr << "    Message: " << msg << "\n";                            \
            scai::common::AssertException::addCallStack( errorStr );               \
            throw scai::common::AssertException( errorStr.str() );                 \
        }                                                                          \
    }

/**
 * @brief The macro SCAI_ASSERT_OPERATOR is used to create the other assert macros
 *
 * @param[in] operator  the operator used to evaluate the expressions
 * @param[in] exp1      first expression for comparison
 * @param[in] exp2      second expression for comparison
 * @param[in] msg       message to indicate reason for the exception
 * @throws    common::Exception (derived from std::exception)
 *
 * Attention: the operator must be defined for the two expressions.
 *
 * The message text will also contain the expressions and their values as string,
 * file location and the current call stack.
 */

#define SCAI_ASSERT_OPERATOR( operator, exp1, exp2, msg )                          \
    {                                                                                  \
        if ( ! ( ( exp1 ) operator ( exp2 ) ) )                                        \
        {                                                                              \
            std::ostringstream errorStr;                                               \
            errorStr << "Assert exp_1 "#operator" exp_2 failed in line " << __LINE__;  \
            errorStr << " of file " << __FILE__ << "\n";                               \
            errorStr << "    Message: " << msg << "\n";                                \
            errorStr << "    exp_1: " << #exp1 " = " << exp1 << "\n";                  \
            errorStr << "    exp_2: " << #exp2 " = " << exp2 << "\n";                  \
            scai::common::AssertException::addCallStack( errorStr );                   \
            throw scai::common::AssertException( errorStr.str() );                     \
        }                                                                              \
    }

#define SCAI_ASSERT_VALID_INDEX( index, size, msg )                                    \
    {                                                                                  \
        if ( ! ( scai::common::Utils::validIndex( index, size ) ) )                    \
        {                                                                              \
            std::ostringstream errorStr;                                               \
            errorStr << "Assert valid index failed in line " << __LINE__;              \
            errorStr << " of file " << __FILE__ << "\n";                               \
            errorStr << "    Message: " << msg << "\n";                                \
            errorStr << "    index: " << #index " = " << index << "\n";                \
            errorStr << "    size : " << #size " = " << size << "\n";                  \
            scai::common::AssertException::addCallStack( errorStr );                   \
            throw scai::common::AssertException( errorStr.str() );                     \
        }                                                                              \
    }

#define SCAI_ASSERT_EQUAL( exp1, exp2, msg )       \
    SCAI_ASSERT_OPERATOR( ==, exp1, exp2, msg )

#define SCAI_ASSERT_UNEQUAL( exp1, exp2, msg )     \
    SCAI_ASSERT_OPERATOR( !=, exp1, exp2, msg )

#define SCAI_ASSERT_LT( exp1, exp2, msg )          \
    SCAI_ASSERT_OPERATOR( <, exp1, exp2, msg )

#define SCAI_ASSERT_LE( exp1, exp2, msg )          \
    SCAI_ASSERT_OPERATOR( <=, exp1, exp2, msg )

#define SCAI_ASSERT_GT( exp1, exp2, msg )          \
    SCAI_ASSERT_OPERATOR( >, exp1, exp2, msg )

#define SCAI_ASSERT_GE( exp1, exp2, msg )          \
    SCAI_ASSERT_OPERATOR( >=, exp1, exp2, msg )

/*
 * Check if ASSERT_LEVEL is set, if not print warning and use default ASSERT_LEVEL (ERROR)
 */
#if !defined( SCAI_ASSERT_LEVEL_OFF ) && !defined( SCAI_ASSERT_LEVEL_DEBUG ) && !defined( SCAI_ASSERT_LEVEL_ERROR )
// turned off for master branch
// #pragma message( "Please specify SCAI_ASSERT_LEVEL_xxx with xxx = DEBUG, ERROR, or OFF" )
// #pragma message( "Will use SCAI_ASSERT_LEVEL_ERROR by default." )

#define SCAI_ASSERT_LEVEL_ERROR
#endif

#ifndef SCAI_ASSERT_LEVEL_OFF
#ifndef SCAI_CHECK_ASSERTS
/**
 * @brief The macro SCAI_CHECK_ASSERTS is used to control the checking of
 *        asserts.
 *
 * The macro LAMACHECKASSERTS is used to control the checking of asserts.
 * LAMACHECKASSERTS will be automatically defined if NDEBUG (for debug
 * builds) is not defined.
 * If LAMACHECKASSERTS is not defined the assertions will not be checked.
 */
#define SCAI_CHECK_ASSERTS
#endif
#endif // NDEBUG

/*
 * Definition of Assert levels
 */

/*
 * SCAI_ASSERT_LEVEL = OFF
 */
#if defined( SCAI_ASSERT_LEVEL_OFF )
#define SCAI_ASSERT_DEBUG( exp, msg ) UNUSED_EXP( exp ); UNUSED_STRING( msg );
#define SCAI_ASSERT_ERROR( exp, msg ) UNUSED_EXP( exp ); UNUSED_STRING( msg );

#define SCAI_ASSERT_EQUAL_DEBUG( exp1, exp2 ) UNUSED_EXP2( exp1, exp2 );
#define SCAI_ASSERT_EQUAL_ERROR( exp1, exp2 ) UNUSED_EXP2( exp1, exp2 );

#define SCAI_ASSERT_EQ_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );
#define SCAI_ASSERT_EQ_ERROR( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );

#define SCAI_ASSERT_NE_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );
#define SCAI_ASSERT_NE_ERROR( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );

#define SCAI_ASSERT_LT_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );
#define SCAI_ASSERT_LT_ERROR( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );

#define SCAI_ASSERT_LE_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );
#define SCAI_ASSERT_LE_ERROR( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );

#define SCAI_ASSERT_GT_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );
#define SCAI_ASSERT_GT_ERROR( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );

#define SCAI_ASSERT_GE_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );
#define SCAI_ASSERT_GE_ERROR( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 ); UNUSED_STRING( msg );

#define SCAI_ASSERT_VALID_INDEX_DEBUG( index, size, msg ) UNUSED_EXP2( index, size ); UNUSED_STRING( msg );
#define SCAI_ASSERT_VALID_INDEX_ERROR( index, size, msg ) UNUSED_EXP2( index, size ); UNUSED_STRING( msg );

/*
 * SCAI_ASSERT_LEVEL = ERROR
 */
#elif defined( SCAI_ASSERT_LEVEL_ERROR )

#define SCAI_ASSERT_DEBUG( exp, msg ) UNUSED_EXP( exp );
#define SCAI_ASSERT_ERROR( exp, msg ) SCAI_ASSERT( exp, msg );

#define SCAI_ASSERT_EQUAL_DEBUG( exp1, exp2 ) UNUSED_EXP2( exp1, exp2 );
#define SCAI_ASSERT_EQUAL_ERROR( exp1, exp2 ) SCAI_ASSERT_EQUAL( exp1, exp2, "" );

#define SCAI_ASSERT_EQ_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 );
#define SCAI_ASSERT_EQ_ERROR( exp1, exp2, msg ) SCAI_ASSERT_EQUAL( exp1, exp2, msg );

#define SCAI_ASSERT_NE_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 );
#define SCAI_ASSERT_NE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_UNEQUAL( exp1, exp2, msg );

#define SCAI_ASSERT_LT_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 );
#define SCAI_ASSERT_LT_ERROR( exp1, exp2, msg ) SCAI_ASSERT_LT( exp1, exp2, msg );

#define SCAI_ASSERT_LE_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 );
#define SCAI_ASSERT_LE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_LE( exp1, exp2, msg );

#define SCAI_ASSERT_GT_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 );
#define SCAI_ASSERT_GT_ERROR( exp1, exp2, msg ) SCAI_ASSERT_GT( exp1, exp2, msg );

#define SCAI_ASSERT_GE_DEBUG( exp1, exp2, msg ) UNUSED_EXP2( exp1, exp2 );
#define SCAI_ASSERT_GE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_GE( exp1, exp2, msg );

#define SCAI_ASSERT_VALID_INDEX_DEBUG( index, size, msg ) UNUSED_EXP2( index, size );
#define SCAI_ASSERT_VALID_INDEX_ERROR( index, size, msg ) SCAI_ASSERT_VALID_INDEX( index, size, msg );

/*
 * SCAI_ASSERT_LEVEL = DEBUG
 */
#elif defined( SCAI_ASSERT_LEVEL_DEBUG )

#define SCAI_ASSERT_DEBUG( exp, msg ) SCAI_ASSERT( exp, msg );
#define SCAI_ASSERT_ERROR( exp, msg ) SCAI_ASSERT( exp, msg );

#define SCAI_ASSERT_EQUAL_DEBUG( exp1, exp2 ) SCAI_ASSERT_EQUAL( exp1, exp2, "" );
#define SCAI_ASSERT_EQUAL_ERROR( exp1, exp2 ) SCAI_ASSERT_EQUAL( exp1, exp2, "" );

#define SCAI_ASSERT_EQ_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_EQUAL( exp1, exp2, msg );
#define SCAI_ASSERT_EQ_ERROR( exp1, exp2, msg ) SCAI_ASSERT_EQUAL( exp1, exp2, msg );

#define SCAI_ASSERT_NE_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_UNEQUAL( exp1, exp2, msg );
#define SCAI_ASSERT_NE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_UNEQUAL( exp1, exp2, msg );

#define SCAI_ASSERT_LT_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_LT( exp1, exp2, msg );
#define SCAI_ASSERT_LT_ERROR( exp1, exp2, msg ) SCAI_ASSERT_LT( exp1, exp2, msg );

#define SCAI_ASSERT_LE_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_LE( exp1, exp2, msg );
#define SCAI_ASSERT_LE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_LE( exp1, exp2, msg );

#define SCAI_ASSERT_GT_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_GT( exp1, exp2, msg );
#define SCAI_ASSERT_GT_ERROR( exp1, exp2, msg ) SCAI_ASSERT_GT( exp1, exp2, msg );

#define SCAI_ASSERT_GE_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_GE( exp1, exp2, msg );
#define SCAI_ASSERT_GE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_GE( exp1, exp2, msg );

#define SCAI_ASSERT_VALID_INDEX_DEBUG( index, size, msg ) SCAI_ASSERT_VALID_INDEX( index, size, msg );
#define SCAI_ASSERT_VALID_INDEX_ERROR( index, size, msg ) SCAI_ASSERT_VALID_INDEX( index, size, msg );

#endif
