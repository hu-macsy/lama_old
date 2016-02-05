/**
 * @file scai/common/test/TestMacros.hpp
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
 * @brief Additional Macros used for testing of LAMA with Boost Test.
 * @author Jiri Kraus
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once

// local library
#include <scai/common/Printable.hpp>
#include <scai/common/exception/UnsupportedException.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

// boost
#include <boost/assign/list_of.hpp>
#include <boost/preprocessor.hpp>

#include <boost/test/detail/unit_test_parameters.hpp>

// std
#include <string>
#include <map>

/*
 * log levels are defined in boost/test/detail/log_level.hpp
 */

#define IF_LOG_LEVEL_IS_TEST_SUITE \
	if ( boost::unit_test::runtime_config::log_level() == boost::unit_test::log_test_units )

//template<typename ValueType>
//inline scai::lama::Scalar scalarEps();
//
//template<>
//inline scai::lama::Scalar scalarEps<float>()
//{
//    return scai::lama::Scalar( 1E-8f );
//}
//
//template<>
//inline scai::lama::Scalar scalarEps<double>()
//{
//    return scai::lama::Scalar( 1E-16 );
//}

/*
 * @brief HelperMacro SCAI_CHECK_SCALAR_CLOSE( x, y, type, percent_eps )
 *
 * Extended macro BOOST_CHECK_CLOSE(left,right,procent_tolerance) from Boost.Test.
 * Checks if the difference between x and y is smaller then eps.
 * Transform Scalar into ValueType, calls BOOST_CHECK_CLOSE with ValueTypes.
 *
 * @param x             Scalar
 * @param y             Scalar
 * @param type          value type to be used for check
 * @param percent_eps   Epsilon[%]
 *
 */

#define SCAI_CHECK_CLOSE( x, y, tolerance ) 																	    \
	{																											    \
		ValueType diff = (x) - (y);																				    \
		BOOST_CHECK( scai::common::Math::abs( diff ) < static_cast<ValueType>( tolerance ) ) ;		\
	}

/*
 * @brief HelperMacro LAMA_WRITEAT_TEST( printable )
 *
 * This macro checks if a output will be created by writing an object
 * into a stream. The length of this output must be greater than 0.
 * This object must be inherited from class Printable.
 *
 * @param printable     object of type printable
 */

#define LAMA_WRITEAT_TEST( printable )                                                                                 \
    {                                                                                                                  \
        scai::common::Printable* p = dynamic_cast<scai::common::Printable*>( &printable );                             \
        std::stringstream mStream;                                                                                     \
        p->writeAt( mStream );                                                                                         \
        std::string mString = mStream.str();                                                                           \
        BOOST_CHECK( mString.length() > 0 );                                                                           \
    }

/*
 * @brief HelperMacro LAMA_WRITEAT_PTR_TEST( printable )
 *
 * This macro checks if a output will be created by writing an object
 * into a stream. The length of this output must be greater than 0.
 * This object must be a pointer to an object inherited from class
 * Printable.
 *
 * @param printable     pointer to an object of type printable
 */

#define LAMA_WRITEAT_PTR_TEST( printable )         \
    {                                              \
        std::stringstream mStream;                 \
        printable->writeAt( mStream );             \
        std::string mString = mStream.str();       \
        BOOST_CHECK( mString.length() > 0 );       \
    }

/*
 * @brief HelperMacro LAMA_AUTO_TEST_CASE_T( name, classname )
 *
 * This macro creates a small boost test auto case for all value types (without context).
 * The test case name is based on the name of the given test method.
 *
 * @param name          name of test method, which will invoke.
 * @param classname     name of the given test class.
 */
#define LAMA_AUTO_TEST_CASE_T( name, classname)                                                                        \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
        const std::string lama_name = #name;                                                                           \
        const std::string lama_classname = #classname;                                                                 \
        scai::lama::classname::name<float>( );                                                                         \
        scai::lama::classname::name<double>( );                                                                        \
    }

/*
 * @brief HelperMacro LAMA_AUTO_TEST_CASE_T( name, classname )
 *
 * This macro creates a small boost test auto case for all value types (without context).
 * The test case name is based on the name of the given test method.
 *
 * @param name          name of test method, which will invoke.
 * @param classname     name of the given test class.
 */
#define LAMA_AUTO_TEST_CASE_TL( name, classname)                                                                       \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
        const std::string lama_name = #name;                                                                           \
        const std::string lama_classname = #classname;                                                                 \
        scai::lama::classname::name<float>( logger );                                                                  \
        scai::lama::classname::name<double>( logger );                                                                 \
    }

    /*
     * @brief HelperMacro COMMONTESTCASEINVOKER( object_name, method_name )
     *
     * This macro represents the invoke of a common used test method.
     *
     * @param object_name     name of the object, which consists the test method.
     * @param method_name     name of method that will be called.
     */

#define COMMONTESTCASEINVOKER( object_name, method_name )            \
    {                                                                \
        if ( testcase == #method_name )                              \
        {                                                            \
        	try                                                      \
            {                                                        \
        		object_name.method_name();                           \
            }                                                        \
            catch(scai::common::UnsupportedException const& ex)      \
            {                                                        \
            	std::cout << "failed on " #method_name << std::endl; \
            }                                                        \
        }                                                            \
	}

    /*
     * @brief HelperMacro COMMONTESTCASEINVOKER_TEMPLATE( object_name, method_name, ValueType )
     *
     * This macro represents the invoke of a common used test method.
     *
     * @param object_name     name of the object, which consists the test method.
     * @param method_name     name of method that will be called.
     */

#define COMMONTESTCASEINVOKER_TEMPLATE( object_name, method_name, ValueType )   \
    {                                                                           \
		if ( testcase == #method_name )                                         \
			object_name.method_name<ValueType>();                               \
	}

    /*
     * @brief HelperMacro LAMA_COMMON_TEST_CASE_TM( classname, templatename, methodname )
     *
     * This macro creates a templated test method, which can be used by test cases.
     *
     * @param classname       name of the test class.
     * @param templatename    name of a template.
     * @param methodname      name of the test method, that will be created.
     */
#define LAMA_COMMON_TEST_CASE_TM( classname, templatename, methodname )                                                \
    template<typename templatename>                                                                                    \
    void classname::methodname()                                                                                       \
    {                                                                                                                  \
        const std::string lama_common_testcase_method = #methodname;                                                   \
        IF_LOG_LEVEL_IS_TEST_SUITE                                                                                     \
            BOOST_TEST_MESSAGE( "    Entering common test case \"" + lama_common_testcase_method + "\" " );            \
        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE_TM_END()
         *
         * This macro closes this testmethod.
         *
         */
#define LAMA_COMMON_TEST_CASE_TM_END();                                                                                \
	IF_LOG_LEVEL_IS_TEST_SUITE                                                                                         \
        BOOST_TEST_MESSAGE( "    Leaving common test case \"" + lama_common_testcase_method + "\" " );                 \
    }

        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE( classname, methodname )
         *
         * This macro creates a test method, which is used by many test cases.
         *
         * @param classname       name of the test class.
         * @param methodname      name of the test method, that will be created.
         */
#define LAMA_COMMON_TEST_CASE( classname, methodname )                                                                 \
    void classname::methodname()                                                                                       \
    {                                                                                                                  \
        const std::string lama_common_testcase_method = #methodname;                                                   \
        IF_LOG_LEVEL_IS_TEST_SUITE                                                                                     \
            BOOST_TEST_MESSAGE( "    Entering common test case \"" + lama_common_testcase_method + "\" " );            \

        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE_END()
         *
         * This macro closes this test method.
         *
         */
#define LAMA_COMMON_TEST_CASE_END();                                                                                   \
	IF_LOG_LEVEL_IS_TEST_SUITE                                                                                   	   \
        BOOST_TEST_MESSAGE( "    Leaving common test case \"" + lama_common_testcase_method + "\" " );                 \
    }

        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE_RUNNER( classname )
         *
         * This macro creates the runTests()-Method to invoke all common test methods.
         *
         */

#define LAMA_COMMON_TEST_CASE_RUNNER( classname );                                                                     \
    void classname::runTests()

        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE( classname, templatename, methodname )
         *
         * This macro creates a test method, which is used by many test cases.
         *
         * @param classname       name of the test class.
         * @param templatename    name of the template.
         * @param methodname      name of the test method, that will be created.
         */

#define LAMA_COMMON_TEST_CASE_TEMPLATE( classname, templatename, methodname )                                          \
    template<typename templatename>                                                                                    \
    void classname<templatename>::methodname()                                                                         \
    {                                                                                                                  \
        std::string lama_common_testcase_method = #methodname;                                                         \
        std::ostringstream omsg;                                                                                       \
        omsg << scai::common::getScalarType<templatename>();                                                           \
        std::string lama_common_testcase_template = omsg.str();                                                        \
        IF_LOG_LEVEL_IS_TEST_SUITE                                                                                     \
        {                                                                                                              \
            SCAI_LOG_INFO( logger, "    Entering common test case \"" + lama_common_testcase_method + "<" +            \
                           lama_common_testcase_template + ">\" " );                                                   \
        }



        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE_END()
         *
         * This macro closes this test method.
         *
         */
#define LAMA_COMMON_TEST_CASE_TEMPLATE_END();                                                                    \
		IF_LOG_LEVEL_IS_TEST_SUITE                                                                               \
        	SCAI_LOG_INFO( logger, "    Leaving common test case \"" + lama_common_testcase_method + "\" " );    \
    }

        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE_RUNNER( classname )
         *
         * This macro creates the runTests()-Method to invoke all common test methods.
         *
         */

#define LAMA_COMMON_TEST_CASE_RUNNER_TEMPLATE( classname );                                                            \
    template<typename StorageType>                                                                                     \
    void classname<StorageType>::runTests()

#if defined(SCAI_ASSERT_LEVEL_OFF)

// Do not check for exception if ASSERTions are switched off

#define SCAI_CHECK_THROW( stmt, exception )

#else

// use Boost macro for throwing exception

#define SCAI_CHECK_THROW( stmt, exception )    \
    BOOST_CHECK_THROW( stmt, exception )

#endif
