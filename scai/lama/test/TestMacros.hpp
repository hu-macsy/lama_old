/**
 * @file scai/lama/test/TestMacros.hpp
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

#include <scai/common/test/TestMacros.hpp>
#include <scai/hmemo/test/TestMacros.hpp>
#include <scai/kregistry/test/TestMacros.hpp>

#include <scai/lama/Scalar.hpp>

#include <boost/test/detail/unit_test_parameters.hpp>

#ifdef SCAI_CHECK_CLOSE
    #undef  SCAI_CHECK_CLOSE
#endif

/*
 * @brief HelperMacro SCAI_CHECK_CLOSE( x, y, tolerance )
 *
 * Extends macro BOOST_CHECK_CLOSE( x, y, tolerance) from Boost.Test.
 * Extends macro SCAI_CHECK_CLOSE( x, y, tolerance) from common
 * now works also if x or y is a Scalar
 *
 * @param ValueType     value type to be used for check
 * @param x             Scalar
 * @param y             Scalar
 * @param percent_eps   Epsilon[%]
 *
 */

#define SCAI_CHECK_CLOSE( x, y, tolerance )                                                                   \
    {                                                                                                         \
        scai::lama::Scalar xScalar = scai::lama::Scalar( x );                                                 \
        scai::lama::Scalar yScalar = scai::lama::Scalar( y );                                                 \
        ScalarRepType xVal = xScalar.getValue<ScalarRepType>();                                               \
        ScalarRepType yVal = yScalar.getValue<ScalarRepType>();                                               \
        BOOST_CHECK_CLOSE( scai::common::Math::real( xVal ), scai::common::Math::real( yVal ), tolerance );   \
        BOOST_CHECK_CLOSE( scai::common::Math::imag( xVal ), scai::common::Math::imag( yVal ), tolerance ) ;  \
    }

#define SCAI_CHECK_SMALL( x, ValueType, eps )                   \
        BOOST_CHECK_SMALL( x, static_cast<ValueType>(eps) );    \

#define SCAI_CHECK_SMALL_EPS( x, ValueType )                                        \
    SCAI_CHECK_SMALL( x, ValueType, scai::common::TypeTraits<ValueType>::small() )

/*
 * @brief HelperMacro SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps )
 *
 * Extended Macro BOOST_CHECK_SMALL( x, eps ) from Boost.Test for
 * Scalar class of LAMA. Transforms Scalar x into ValueType,
 * and calls BOOST_CHECK_SMALL with arguments of type ValueType.
 *
 * @param x             Scalar
 * @param ValueType     type of Scalar x used for test
 * @param eps           Epsilon
 *
 * Static cast is used to convert eps to the right ValueType.
 */

/*#define SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps )                     \
    {                                                                    \
        ValueType xHelper = (x).getValue<ValueType >();                  \
        BOOST_CHECK_SMALL( xHelper, static_cast<ValueType >( eps ) );    \
    }*/

#define SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps )                     \
        SCAI_CHECK_SMALL( (x).getValue<ValueType>(), ValueType, eps )    \

/*
 * @brief HelperMacro SCAI_CHECK_SCALAR_SMALL_EPS( x, ValueType )
 *
 * Same as SCAI_CHECK_SCALAR_SMALL but with default eps value.
 *
 * @param x             Scalar
 * @param ValueType     type of Scalar to be used for test
 */

#define SCAI_CHECK_SCALAR_SMALL_EPS( x, ValueType )                  \
    SCAI_CHECK_SCALAR_SMALL( x, ValueType, scai::common::TypeTraits<ValueType>::small() )

/*
 * log levels are defined in boost/test/detail/log_level.hpp
 */

#define IF_LOG_LEVEL_IS_TEST_SUITE \
    if ( boost::unit_test::runtime_config::log_level() == boost::unit_test::log_test_units )

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
    IF_LOG_LEVEL_IS_TEST_SUITE                                                                                         \
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


#define LAMA_COMMON_TEST_CASE_RUNNER_TEMPLATE( classname );                                                            \
    template<typename StorageType>                                                                                     \
    void classname<StorageType>::runTests()
