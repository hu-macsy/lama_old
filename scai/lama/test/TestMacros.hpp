/**
 * @file scai/lama/test/TestMacros.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Additional Macros used for testing of LAMA with Boost Test.
 * @author Jiri Kraus
 * @date 06.04.2011
 */

#pragma once

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/typeloop.hpp>
#include <scai/logging.hpp>
#include <scai/hmemo/Context.hpp>
#include <scai/kregistry/test/TestMacros.hpp>
#include <scai/lama/Scalar.hpp>

#if BOOST_VERSION > 105900
	#include <boost/test/unit_test_parameters.hpp>
#else
	#include <boost/test/detail/unit_test_parameters.hpp>
#endif

#include <boost/test/test_tools.hpp>
#include <boost/test/detail/global_typedef.hpp>
#include <boost/test/detail/log_level.hpp>

#include <boost/test/detail/suppress_warnings.hpp>


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

/*
#define SCAI_CHECK_SMALL( x, ValueType, eps )                   \
        BOOST_CHECK_SMALL( x, static_cast<ValueType>(eps) );    \
*/

#define SCAI_CHECK_SMALL( x, ValueType, eps )                   \
        BOOST_CHECK( x < static_cast<ValueType>(eps) );    \

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

#if BOOST_VERSION > 105900
    #define IF_LOG_LEVEL_IS_TEST_SUITE \
    	if ( boost::unit_test::runtime_config::LOG_LEVEL == "test_suite" )
#else
    #define IF_LOG_LEVEL_IS_TEST_SUITE \
	if ( boost::unit_test::runtime_config::log_level() == boost::unit_test::log_test_units )
#endif 
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


/** This macro runs method<ValueType>( context ) where ValueType is given by the I-th arithmetic Host Type
 *  (skips the run for long double types on CUDA as not supported there).
 *
 *  KernelRegistryException is caught with a correpsonding warn message on logger
 */

#define LAMA_RUN_TEST( method, ValueType )                                                                      \
    try                                                                                                         \
    {                                                                                                           \
        if ( context->getType() == scai::common::context::CUDA                                                  \
                               || context->getType() == scai::common::context::MIC )                            \
        {                                                                                                       \
            switch( scai::common::getScalarType<ValueType>() )                                                  \
            {                                                                                                   \
                case scai::common::scalar::LONG_DOUBLE:                                                         \
                case scai::common::scalar::LONG_DOUBLE_COMPLEX:                                                 \
                    return;                                                                                     \
                default:                                                                                        \
                     ;                                                                                          \
             }                                                                                                  \
        }                                                                                                       \
        method<ValueType>( context );                                                                           \
    }                                                                                                           \
    catch ( scai::kregistry::KernelRegistryException& )                                                         \
    {                                                                                                           \
        SCAI_LOG_WARN( logger, #method << "<" << PRINT_STRING( ValueType ) << "> cannot run on "                \
                       << context->getType() << ", corresponding function not implemented yet."        )        \
        return;                                                                                                 \
    }

/*
 * @brief HelperMacro LAMA_AUTO_TEST_CASE_CT( name, classname )
 *
 * This macro creates a boost test auto case, which uses all possible contexts.
 * The test case name is based on the name of the given test method.
 *
 * @param name          name of test method, which will invoke.
 * @param classname     name of the given test class.
*/
#define LAMA_AUTO_TEST_CASE_CT( name, classname, namespacename )                                                                        \
                                                                                                                                        \
    BOOST_AUTO_TEST_CASE( name )                                                                                                        \
    {                                                                                                                                   \
            ContextPtr context = Context::getContextPtr();                                                                              \
            const std::string lama_name = #name;                                                                                        \
            const std::string lama_classname = #classname;                                                                              \
            SCAI_COMMON_TYPELOOP_LVL2( SCAI_ARITHMETIC_HOST_CNT, namespacename::classname::name, LAMA_RUN_TEST, SCAI_ARITHMETIC_HOST )  \
    }

