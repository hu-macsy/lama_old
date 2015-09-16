/*
 * Assert.hpp
 *
 *  Created on: Aug 31, 2015
 *      Author: eschricker
 */

#pragma once

// local library
#include <scai/common/exception/AssertException.hpp>

template<typename ValueType>
static inline void unused( ValueType const )
{
    // do nothing
}

template<typename ValueType1,typename ValueType2>
static inline void unused( const ValueType1&, const ValueType2& )
{
    // do nothing
}

#define UNUSED_STRING( ... )                                                \
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
	#pragma message( "Please specify SCAI_ASSERT_LEVEL_xxx with xxx = DEBUG, ERROR, or OFF" )
	#pragma message( "Will use SCAI_ASSERT_LEVEL_ERROR by default." )

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
	#define SCAI_ASSERT_DEBUG( exp, msg ) unused( exp ); UNUSED_STRING( msg );
	#define SCAI_ASSERT_ERROR( exp, msg ) unused( exp ); UNUSED_STRING( msg );

	#define SCAI_ASSERT_EQUAL_DEBUG( exp1, exp2 ) unused( exp1, exp2 );
	#define SCAI_ASSERT_EQUAL_ERROR( exp1, exp2 ) unused( exp1, exp2 );

	#define SCAI_ASSERT_UNEQUAL_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );
	#define SCAI_ASSERT_UNEQUAL_ERROR( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );

	#define SCAI_ASSERT_LT_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );
	#define SCAI_ASSERT_LT_ERROR( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );

	#define SCAI_ASSERT_LE_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );
	#define SCAI_ASSERT_LE_ERROR( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );

	#define SCAI_ASSERT_GT_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );
	#define SCAI_ASSERT_GT_ERROR( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );

	#define SCAI_ASSERT_GE_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );
	#define SCAI_ASSERT_GE_ERROR( exp1, exp2, msg ) unused( exp1, exp2 ); UNUSED_STRING( msg );

/*
 * SCAI_ASSERT_LEVEL = ERROR
 */
#elif defined( SCAI_ASSERT_LEVEL_ERROR )

	#define SCAI_ASSERT_DEBUG( exp, msg ) unused( exp );
	#define SCAI_ASSERT_ERROR( exp, msg ) SCAI_ASSERT( exp, msg );

	#define SCAI_ASSERT_EQUAL_DEBUG( exp1, exp2 ) unused( exp1, exp2 );
	#define SCAI_ASSERT_EQUAL_ERROR( exp1, exp2 ) SCAI_ASSERT_EQUAL( exp1, exp2, "" );

	#define SCAI_ASSERT_UNEQUAL_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 );
	#define SCAI_ASSERT_UNEQUAL_ERROR( exp1, exp2, msg ) SCAI_ASSERT_UNEQUAL( exp1, exp2, msg );

	#define SCAI_ASSERT_LT_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 );
	#define SCAI_ASSERT_LT_ERROR( exp1, exp2, msg ) SCAI_ASSERT_LT( exp1, exp2, msg );

	#define SCAI_ASSERT_LE_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 );
	#define SCAI_ASSERT_LE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_LE( exp1, exp2, msg );

	#define SCAI_ASSERT_GT_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 );
	#define SCAI_ASSERT_GT_ERROR( exp1, exp2, msg ) SCAI_ASSERT_GT( exp1, exp2, msg );

	#define SCAI_ASSERT_GE_DEBUG( exp1, exp2, msg ) unused( exp1, exp2 );
	#define SCAI_ASSERT_GE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_GE( exp1, exp2, msg );

/*
 * SCAI_ASSERT_LEVEL = DEBUG
 */
#elif defined( SCAI_ASSERT_LEVEL_DEBUG )

	#define SCAI_ASSERT_DEBUG( exp, msg ) SCAI_ASSERT( exp, msg );
	#define SCAI_ASSERT_ERROR( exp, msg ) SCAI_ASSERT( exp, msg );

	#define SCAI_ASSERT_EQUAL_DEBUG( exp1, exp2 ) SCAI_ASSERT_EQUAL( exp1, exp2, "" );
	#define SCAI_ASSERT_EQUAL_ERROR( exp1, exp2 ) SCAI_ASSERT_EQUAL( exp1, exp2, "" );

	#define SCAI_ASSERT_UNEQUAL_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_UNEQUAL( exp1, exp2, msg );
	#define SCAI_ASSERT_UNEQUAL_ERROR( exp1, exp2, msg ) SCAI_ASSERT_UNEQUAL( exp1, exp2, msg );

	#define SCAI_ASSERT_LT_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_LT( exp1, exp2, msg );
	#define SCAI_ASSERT_LT_ERROR( exp1, exp2, msg ) SCAI_ASSERT_LT( exp1, exp2, msg );

	#define SCAI_ASSERT_LE_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_LE( exp1, exp2, msg );
	#define SCAI_ASSERT_LE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_LE( exp1, exp2, msg );

	#define SCAI_ASSERT_GT_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_GT( exp1, exp2, msg );
	#define SCAI_ASSERT_GT_ERROR( exp1, exp2, msg ) SCAI_ASSERT_GT( exp1, exp2, msg );

	#define SCAI_ASSERT_GE_DEBUG( exp1, exp2, msg ) SCAI_ASSERT_GE( exp1, exp2, msg );
	#define SCAI_ASSERT_GE_ERROR( exp1, exp2, msg ) SCAI_ASSERT_GE( exp1, exp2, msg );
#endif
