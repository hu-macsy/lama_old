/**
 * @file scai/hmemo/test/TestMacros.hpp
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
 * @brief Additional Macros used for testing of hmemo with Boost Test.
 * @author Jiri Kraus
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/HArray.hpp>

#include <scai/common/ContextType.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <boost/test/unit_test.hpp>

/* --------------------------------------------------------------------- */

/** Fixture to be used for BOOST_GLOBAL_FIXTURE     
 *
 *  provides access to testContext used as context at which tests should run
 *
 *  Purpose: use global Fixture avoids init/free of context device for each single test
 *
 *  Note: static variable ContextFix::testContext must be defined in cpp file.
 */
struct ContextFix
{
    ContextFix()
    {   
        testContext = scai::hmemo::Context::getContextPtr();
        // BOOST_TEST_MESSAGE( "Setup ContextFix: test context = " << *testContext ); 
    }

    ~ContextFix()
    {
        // BOOST_TEST_MESSAGE( "Teardown ContextFix" ); 
        testContext.reset();
    }
    
    static scai::hmemo::ContextPtr testContext;
};

// ***********************************************************
// ***********************************************************
// ***********************************************************
// Rest of this file is old stuff and considered to be deleted
// ***********************************************************
// ***********************************************************
// ***********************************************************

template<typename ValueType>
static inline void initArrayVals( scai::hmemo::HArray<ValueType>& dst, const ValueType src[], const IndexType size)
{
    dst.init( src, size );
}

template<typename ValueType>
static inline void initArrayVal( scai::hmemo::HArray<ValueType>& dst, const ValueType value, const IndexType size)
{
    dst.init( value, size );
}

/** This macro runs method<ValueType>( context ) where ValueType is given by the I-th arithmetic Host Type
 *  (skips the run for long double types on CUDA as not supported there).
 *
 *  KernelRegistryException is caught with a correpsonding warn message on logger
 */

#define LAMA_RUN_TEST(z, I, method )                                                                            \
    try                                                                                                         \
    {                                                                                                           \
        if ( context->getType() == scai::common::context::CUDA                                                  \
                               || context->getType() == scai::common::context::MIC )                                           \
        {                                                                                                       \
            switch( scai::common::getScalarType<ARITHMETIC_HOST_TYPE_##I>() )                                   \
            {                                                                                                   \
                case scai::common::scalar::LONG_DOUBLE:                                                         \
                case scai::common::scalar::LONG_DOUBLE_COMPLEX:                                                 \
                    return;                                                                                     \
                default:                                                                                        \
                     ;                                                                                          \
             }                                                                                                  \
        }                                                                                                       \
        method<ARITHMETIC_HOST_TYPE_##I>( context );                                                            \
    }                                                                                                           \
    catch ( scai::kregistry::KernelRegistryException& )                                                         \
    {                                                                                                           \
        SCAI_LOG_WARN( logger, #method << "<" << PRINT_STRING( ARITHMETIC_HOST_TYPE_##I ) << "> cannot run on " \
                       << context->getType() << ", corresponding function not implemented yet."        )        \
        return;                                                                                                 \
    }

/*
 * @brief HelperMacro LAMA_AUTO_TEST_CASE( name, classname )
 *
 * This macro creates a boost test auto case, which uses all possible contexts.
 * The test case name is based on the name of the given test method.
 *
 * @param name          name of test method, which will invoke.
 * @param classname     name of the given test class.
 */
#define LAMA_AUTO_TEST_CASE( name, classname )                                                                         \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
            ContextPtr context = Context::getContextPtr();                                                             \
            IF_LOG_LEVEL_IS_TEST_SUITE                                                                                 \
                SCAI_LOG_INFO( logger, "    Entering context: " << context->getType() );                               \
			scai::lama::classname::name( context );                                                                    \
        }                                                                                                              \
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
#define LAMA_AUTO_TEST_CASE_CT( name, classname, namespacename )                                        \
                                                                                                        \
    BOOST_AUTO_TEST_CASE( name )                                                                        \
    {                                                                                                   \
            ContextPtr context = Context::getContextPtr();                                              \
            const std::string lama_name = #name;                                                        \
            const std::string lama_classname = #classname;                                              \
            BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_RUN_TEST, namespacename::classname::name )  \
    }

    /*
     * @brief HelperMacro LAMA_AUTO_TEST_CASE_CTL( name, classname )
     *
     * This macro creates a boost test auto case, which uses all possible contexts.
     * The test case name is based on the name of the given test method.
     *
     * @param name          name of test method, which will invoke.
     * @param classname     name of the given test class.
     */
#define LAMA_AUTO_TEST_CASE_CTL( name, classname )                                                                     \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
            ContextPtr context = Context::getContextPtr();                                                             \
            BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE__CNT, LAMA_RUN_TESTL, scai::lama::classname::name )                  \
    }


    /*
     * @brief HelperMacro LAMA_AUTO_TEST_CASE_CTT( name, classname )
     *
     * This macro creates a boost test auto case, which uses all possible contexts.
     * The test case name is based on the name of the given test method. All combinations
     * of ValueTypes (at the moment float and double) will be executed.
     *
     * @param name       name of test method, which will invoke.
     * @param classname  name of the given test class.
     */

#define LAMA_AUTO_TEST_CASE_CTT( name, classname )                                                                     \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
            ContextPtr context = Context::getContextPtr();                                                             \
            IF_LOG_LEVEL_IS_TEST_SUITE                                                                                 \
            {                                                                                                          \
                SCAI_LOG_INFO( logger, "    Entering context: " << context->getType() );                               \
            }                                                                                                          \
            const std::string lama_name = #name;                                                                       \
            const std::string lama_classname = #classname;                                                             \
            try                                                                                                        \
            {                                                                                                          \
                       scai::lama::classname::name<float, float>( context );                                           \
            }                                                                                                          \
            catch ( scai::kregistry::KernelRegistryException& )                                                        \
            {                                                                                                          \
                SCAI_LOG_WARN( logger, lama_classname << "::" << lama_name << "<float, float> cannot run on  "         \
                               << context->getType() << ", corresponding function not implemented yet." );             \
                return;                                                                                                \
            }                                                                                                          \
            try                                                                                                        \
            {                                                                                                          \
                       scai::lama::classname::name<double, double>( context );                                         \
            }                                                                                                          \
            catch ( scai::kregistry::KernelRegistryException& )                                                        \
            {                                                                                                          \
                SCAI_LOG_WARN( logger, lama_classname << "::" << lama_name << "<double, double> cannot run on  "       \
                               << context->getType() << ", corresponding function not implemented yet." );             \
                return;                                                                                                \
            }                                                                                                          \
            try                                                                                                        \
            {                                                                                                          \
                       scai::lama::classname::name<float, double>( context );                                          \
            }                                                                                                          \
            catch ( scai::kregistry::KernelRegistryException& )                                                        \
            {                                                                                                          \
                SCAI_LOG_WARN( logger, lama_classname << "::" << lama_name << "<float, double> cannot run on  "        \
                               << context->getType() << ", corresponding function not implemented yet." );             \
                return;                                                                                                \
            }                                                                                                          \
            try                                                                                                        \
            {                                                                                                          \
                       scai::lama::classname::name<double, float>( context );                                          \
            }                                                                                                          \
            catch ( scai::kregistry::KernelRegistryException& )                                                        \
            {                                                                                                          \
                SCAI_LOG_WARN( logger, lama_classname << "::" << lama_name << "<double, float> cannot run on  "        \
                               << context->getType() << ", corresponding function not implemented yet." );             \
                return;                                                                                                \
            }                                                                                                          \
    }

    /*
     * @brief HelperMacro LAMA_AUTO_TEST_CASE_CTDUMMY( name, classname, logger )
     *
     * This macro creates a boost test auto case, which uses all possible contexts. The test case name is based
     * on the name of the given test method. Uses a dummy type to execute the test, that is needed for methods
     * that needs templating for other reasons (usage of the lama interface).
     *
     * @param name          name of test method, which will invoke.
     * @param classname     name of the given test class.
     */
#define LAMA_AUTO_TEST_CASE_CTDUMMY( name, classname )                                                                 \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
            ContextPtr context = Context::getContextPtr();                                                             \
            IF_LOG_LEVEL_IS_TEST_SUITE                                                                                 \
            {                                                                                                          \
                SCAI_LOG_INFO( logger, "    Entering context: " << context->getType() );                               \
            }                                                                                                          \
            const std::string lama_name = #name;                                                                       \
            const std::string lama_classname = #classname;                                                             \
            try                                                                                                        \
            {                                                                                                          \
                       scai::lama::classname::name<name>( context );                                                   \
            }                                                                                                          \
            catch ( scai::kregistry::KernelRegistryException& )                                                        \
            {                                                                                                          \
                SCAI_LOG_WARN( logger, lama_classname << "::" << lama_name << " cannot run on  "                       \
                               << context->getType() << ", corresponding function not implemented yet." );             \
                return;                                                                                                \
            }                                                                                                          \
    }

#define LAMA_AUTO_TEST_CASE_CDUMMY( name, classname )                                                                  \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
            ContextPtr context = Context::getContextPtr();                                                             \
            IF_LOG_LEVEL_IS_TEST_SUITE                                                                                 \
            {                                                                                                          \
                SCAI_LOG_INFO( logger, "    Entering context: " << context->getType() );                               \
            }                                                                                                          \
            const std::string lama_name = #name;                                                                       \
            const std::string lama_classname = #classname;                                                             \
            try                                                                                                        \
            {                                                                                                          \
                scai::lama::classname::name( context );                                                                \
            }                                                                                                          \
            catch ( scai::kregistry::KernelRegistryException& )                                                        \
            {                                                                                                          \
                SCAI_LOG_WARN( logger, lama_classname << "::" << lama_name << " cannot run on  "                       \
                               << context->getType() << ", corresponding function not implemented yet." );             \
                return;                                                                                                \
            }                                                                                                          \
    }
