/**
 * @file TestMacros.hpp
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
 * @brief Additional Macros used for testing of LAMA with Boost Test.
 * @author Jiri Kraus
 * @date 06.04.2011
 * @since 1.0.0
 */

#include <string>
#include <map>

#include <lama/Scalar.hpp>
#include <lama/Printable.hpp>
#include <lama/Context.hpp>

#include <boost/assign/list_of.hpp>

#ifndef LAMA_TESTMACROS_HPP_
#define LAMA_TESTMACROS_HPP_

extern std::string loglevel_argument;

/**
 * @brief eps() returns the desired precision of calculations with T in tests.
 *
 * @return the desired precision for calculations with T in tests
 */
template<typename T>
inline T eps();

/**
 * @brief eps<double>() returns the desired precision of calculations in single
 *        precision.
 *
 * @return the desired precision for calculations with single precision in tests
 */
template<>
inline float eps<float>()
{
    return 1E-5f;
}

/**
 * @brief eps<double>() returns the desired precision of calculations in double
 *        precision.
 *
 * @return the desired precision for calculations with double precision in tests
 */
template<>
inline double eps<double>()
{
    return 1E-5f;
}

template<typename T>
inline lama::Scalar scalarEps();

template<>
inline lama::Scalar scalarEps<float>()
{
    return lama::Scalar( 1E-8f );
}

template<>
inline lama::Scalar scalarEps<double>()
{
    return lama::Scalar( 1E-16 );
}

/**
 * @brief getEnvContext() returns the current context from enviroment variable LAMA_TEST_CONTEXT
 * If this variable is not set or contains a context, which is not available, all contexts will be used.
 *
 * @return the current context from enviroment variable LAMA_TEST_CONTEXT as a string
 */

inline std::string getEnvContext()
{
    std::string context;  // default is using each context

    const char* envContext = getenv( "LAMA_TEST_CONTEXT" );

    if ( !envContext )
    {
        context = "*";
        return context;
    }
    
    context = envContext;

    // to upper for convenience

    for ( std::string::iterator p = context.begin(); context.end() != p; ++p )
    {
        *p = toupper( *p );
    }

    if ( context == "HOST" )
    { 
        context = "Host";    // set back 
    }
    else if ( context == "CUDA" || context == "OPENCL" )
    {
        // good values
    }
    else
    {
        context = "*";
    }
    
    return context;
}

/**
 * @brief mapEnvContexttoContextType() returns the current context as a ContextType from a std::string
 *
 * @return the current context as a ContextType from a std::string
 */

inline lama::ContextType mapEnvContexttoContextType( std::string contextname )
{
    lama::ContextType myContext;
    std::map<std::string,lama::ContextType> contextmap = 
       boost::assign::map_list_of ( "Host", lama::Context::Host )
                                  ( "CUDA", lama::Context::CUDA )
                                  ( "OPENCL", lama::Context::OpenCL )
                                  ( "MaxContext", lama::Context::MaxContext );
    myContext = contextmap[contextname];
    return myContext;
}

/*
 * @brief HelperMakro LAMA_CHECK_SCALAR_CLOSE( x, y, type, percent_eps )
 *
 * Extended Makro BOOST_CHECK_CLOSE(left,right,procent_tolerance) from Boost.Test.
 * Checks if the difference between x and y is smaller then eps.
 * Transform Scalar into ValueType, calls BOOST_CHECK_CLOSE with ValueTypes.
 *
 * @param x             Scalar
 * @param y             Scalar
 * @param type          value type to be used for check
 * @param percent_eps   Epsilon[%]
 *
 */

#define LAMA_CHECK_SCALAR_CLOSE( x, y, ValueType, tolerance )                                                          \
    ValueType xHelper = x.getValue<ValueType>();                                                                   \
    ValueType yHelper = y.getValue<ValueType>();                                                                   \
    BOOST_CHECK_CLOSE( xHelper, yHelper, tolerance );

/*
 * @brief HelperMakro LAMA_CHECK_SCALAR_SMALL( x, ValueType, eps )
 *
 * Extended Makro BOOST_CHECK_SMALL( x, eps ) from Boost.Test for
 * Scalar class of LAMA. Transforms Scalar x into ValueType,
 * and calls BOOST_CHECK_SMALL with arguments of type ValueType.
 *
 * @param x             Scalar
 * @param ValueType     type of Scalar x used for test
 * @param eps           Epsilon
 *
 * Static cast is used to convert eps to the right ValueType.
 */

#define LAMA_CHECK_SCALAR_SMALL( x, ValueType, eps )                                                                   \
    ValueType xHelper = (x).getValue<ValueType>();                                                                 \
    BOOST_CHECK_SMALL( xHelper, static_cast<ValueType>( eps ) );

/*
 * @brief HelperMakro LAMA_CHECK_SCALAR_SMALL_EPS( x, ValueType )
 *
 * Same as LAMA_CHECK_SCALAR_SMALL but with default eps value.
 *
 * @param x             Scalar
 * @param ValueType     type of Scalar to be used for test
 */

#define LAMA_CHECK_SCALAR_SMALL_EPS( x, ValueType )                                                                    \
    LAMA_CHECK_SCALAR_SMALL( x, ValueType, eps<ValueType> () )

/*
 * @brief HelperMakro LAMA_WRITEAT_TEST( printable )
 *
 * This makro checks if a output will be created by writing an object
 * into a stream. The length of this output must be greater than 0.
 * This object must be inherited from class Printable.
 *
 * @param printable     object of type printable
 */

#define LAMA_WRITEAT_TEST( printable )                                                                                 \
    { Printable* p = dynamic_cast<Printable*>( &printable );                                                           \
        std::stringstream mStream;                                                                                         \
        p->writeAt( mStream );                                                                                             \
        std::string mString = mStream.str();                                                                               \
        BOOST_CHECK( mString.length() > 0 ); }

/*
 * @brief HelperMakro LAMA_WRITEAT_PTR_TEST( printable )
 *
 * This makro checks if a output will be created by writing an object
 * into a stream. The length of this output must be greater than 0.
 * This object must be a pointer to an object inherited from class
 * Printable.
 *
 * @param printable     pointer to an object of type printable
 */

#define LAMA_WRITEAT_PTR_TEST( printable )                                                                             \
    { std::stringstream mStream;                                                                                       \
        printable->writeAt( mStream );                                                                                     \
        std::string mString = mStream.str();                                                                               \
        BOOST_CHECK( mString.length() > 0 ); }

/*
 * @brief HelperMakro CONTEXTLOOP()
 *
 * This makro creates a loop iterating throw all existing contexts.
 * All contexts, which are not adressable, will be skipped.
 *
 * If the environment variable LAMA_TEST_CONTEXT or the runtime parameter --context
 * is set to a specific context, just this context will be used.
 *
 */

#define CONTEXTLOOP()                                                                                                  \
    std::list<ContextType> listofcontexts;                                                                             \
    std::list<ContextType>::iterator Iter;                                                                             \
    std::string contexttype;                                                                                           \
    contexttype = getEnvContext();                                                                                     \
    if ( contexttype == "*" )                                                                                          \
    {                                                                                                                  \
        LAMA_LOG_INFO( logger, "LAMA_TEST_CONTEXT is not set or has value '*', so all available contexts will be used." );  \
        for ( ContextType i = Context::Host; i < Context::MaxContext; i = static_cast<ContextType>( i + 1 ) )          \
        {                                                                                                              \
            if ( ContextFactory::hasContext( i ) )                                                                     \
            {                                                                                                          \
                listofcontexts.push_back( i );                                                                         \
                LAMA_LOG_DEBUG( logger, "Context " << i << " is available");                                           \
            }                                                                                                          \
            else                                                                                                       \
                LAMA_LOG_INFO( logger, "The following context will be skipped, because it is not available: " << i );  \
        }                                                                                                              \
    } else {                                                                                                           \
        listofcontexts.push_back( mapEnvContexttoContextType( contexttype ) );                                         \
        LAMA_LOG_INFO( logger, "Environment variable LAMA_TEST_CONTEXT contains context = " << getEnvContext() );       \
    }                                                                                                                  \
    for ( Iter = listofcontexts.begin(); Iter != listofcontexts.end(); Iter++ )

/*
 * @brief HelperMakro GETCONTEXT( loc )
 *
 * This makro creates the adressable context. It is only usable in a CONTEXTLOOP().
 *
 * @param loc          loc is the variable, which contains the context.
 */

#define GETCONTEXT( loc )                                                                                              \
    ContextPtr loc;                                                                                                    \
    loc = ContextFactory::getContext( *Iter );

/*
 * @brief HelperMakro LAMA_AUTO_TEST_CASE( name, classname )
 *
 * This makro creates a boost test auto case, which uses all possible contexts.
 * The test case name is based on the name of the given testmethod.
 *
 * @param name          Name of testmethod, which will invoke
 * @param classname     Name of the given Testclass
 */

#define LAMA_AUTO_TEST_CASE( name, classname )                                                                         \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context );                                                                                     \
            if ( loglevel_argument == "test_suite" )                                                                   \
                BOOST_TEST_MESSAGE( "    Entering context: " << context );                                             \
            lama::classname::name( context );                                                                          \
        }                                                                                                              \
    }

/*
 * @brief HelperMakro LAMA_AUTO_TEST_CASE_TT( name, classname )
 *
 * This makro creates a boost test auto case, which uses all possible contexts.
 * The test case name is based on the name of the given testmethod.
 *
 * @param name          Name of testmethod, which will invoke
 * @param classname     Name of the given Testclass
 */

#define LAMA_AUTO_TEST_CASE_T( name, classname )                                                                       \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context );                                                                                     \
            if ( loglevel_argument == "test_suite" )                                                                   \
                BOOST_TEST_MESSAGE( "    Entering context: " << context->getType() );                                  \
            lama::classname::name<float>( context );                                                                   \
            lama::classname::name<double>( context );                                                                  \
        }                                                                                                              \
    }

/*
 * @brief HelperMakro LAMA_AUTO_TEST_CASE_TT( name )
 *
 * This makro creates a boost test auto case, which uses all possible contexts.
 * The test case name is based on the name of the given testmethod. All combinations
 * of ValueTypes (at the moment float and double) will be executed.
 *
 * @param name     Name of testmethod, which will invoke
 */

#define LAMA_AUTO_TEST_CASE_TT( name, classname )                                                                      \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context );                                                                                     \
            if ( loglevel_argument == "test_suite" )                                                                   \
                BOOST_TEST_MESSAGE( "    Entering context: " << context->getType() );                                  \
            lama::classname::name<float, float>( context );                                                            \
            lama::classname::name<double, double>( context );                                                          \
            lama::classname::name<float, double>( context );                                                           \
            lama::classname::name<double, float>( context );                                                           \
        }                                                                                                              \
    }

/*
 * @brief HelperMakro LAMA_AUTO_TEST_CASE_TDUMMY( name, classname )
 *
 * This makro creates a boost test auto case, which uses all possible contexts. The test case name is based
 * on thename of the given testmethod. Uses a dummy type to execute the test, that is needed for methods
 * that needs templation for other reasons (usage of the lama interface)
 *
 * @param name          Name of testmethod, which will invoke
 * @param classname     Name of the given Testclass
 */

#define LAMA_AUTO_TEST_CASE_TDUMMY( name, classname )                                                                  \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context );                                                                                     \
            if ( loglevel_argument == "test_suite" )                                                                   \
                BOOST_TEST_MESSAGE( "    Entering context: " << context->getType() );                                  \
            lama::classname::name<name>( context );                                                                    \
        }                                                                                                              \
    }

/*
 * @brief HelperMakro COMMONTESTCASEINVOKER( object_name, method_name )
 *
 * This makro represents the invoke of a common used testmethod.
 *
 * @param object_name     Name of the Object, which consists the test_method
 * @param method_name     Name of Method that will be call
 */

#define COMMONTESTCASEINVOKER( object_name, method_name )                                                              \
    { if ( testcase == #method_name ) object_name.method_name(); }

/*
 * @brief HelperMakro COMMONTESTCASEINVOKER_TEMPLATE( object_name, method_name, ValueType )
 *
 * This makro represents the invoke of a common used testmethod.
 *
 * @param object_name     Name of the Object, which consists the test_method
 * @param method_name     Name of Method that will be call
 */

#define COMMONTESTCASEINVOKER_TEMPLATE( object_name, method_name, ValueType )                                          \
    { if ( testcase == #method_name ) object_name.method_name<ValueType>(); }

/*
 * @brief HelperMakro LAMA_COMMON_TEST_CASE_TM( classname, templatename, methodname )
 *
 * This makro creates a Templated testMethod, which can be used by testcases.
 *
 * @param classname       Name of the Testclass
 * @param templatename    Name of a template
 * @param methodname      Name of the Testmethod, that will be created
 */

#define LAMA_COMMON_TEST_CASE_TM( classname, templatename, methodname )                                                \
    template<typename templatename>                                                                                        \
    void classname::methodname()                                                                                           \
    {                                                                                                                      \
        const std::string lama_common_testcase_method = #methodname;                                                       \
        if ( loglevel_argument == "test_suite" )                                                                           \
            BOOST_TEST_MESSAGE( "    Entering common test case \"" + lama_common_testcase_method + "\" " );                \
         
/*
 * @brief HelperMakro LAMA_COMMON_TEST_CASE_TM_END()
 *
 * This makro closes this testmethod.
 *
 */

#define LAMA_COMMON_TEST_CASE_TM_END();                                                                                \
    if ( loglevel_argument == "test_suite" )                                                                           \
        BOOST_TEST_MESSAGE( "    Leaving common test case \"" + lama_common_testcase_method + "\" " );                 \
    }

/*
 * @brief HelperMakro LAMA_COMMON_TEST_CASE( classname, methodname )
 *
 * This makro creates a testmethod, which is used by many test cases.
 *
 * @param classname       Name of the Testclass
 * @param methodname      Name of the Testmethod, that will be created
 */

#define LAMA_COMMON_TEST_CASE( classname, methodname )                                                                 \
    void classname::methodname()                                                                                           \
    {                                                                                                                      \
        const std::string lama_common_testcase_method = #methodname;                                                       \
        if ( loglevel_argument == "test_suite" )                                                                           \
            BOOST_TEST_MESSAGE( "    Entering common test case \"" + lama_common_testcase_method + "\" " );                \
         
/*
 * @brief HelperMakro LAMA_COMMON_TEST_CASE_END()
 *
 * This makro closes this testmethod.
 *
 */

#define LAMA_COMMON_TEST_CASE_END();                                                                                   \
    if ( loglevel_argument == "test_suite" )                                                                           \
        BOOST_TEST_MESSAGE( "    Leaving common test case \"" + lama_common_testcase_method + "\" " );                 \
    }

/*
 * @brief HelperMakro LAMA_COMMON_TEST_CASE_RUNNER( classname )
 *
 * This makro creates the runTests()-Method to invoke all common testmethods.
 *
 */

#define LAMA_COMMON_TEST_CASE_RUNNER( classname );                                                                     \
    void classname::runTests()

/*
 * @brief HelperMakro LAMA_COMMON_TEST_CASE( classname, methodname )
 *
 * This makro creates a testmethod, which is used by many test cases.
 *
 * @param classname       Name of the Testclass
 * @param methodname      Name of the Testmethod, that will be created
 */

#define LAMA_COMMON_TEST_CASE_TEMPLATE( classname, templatename, methodname )                                          \
    extern std::string loglevel_argument;                                                                                  \
    template<typename templatename>                                                                                        \
    void classname<templatename>::methodname()                                                                             \
    {                                                                                                                      \
        std::string lama_common_testcase_method = #methodname;                                                             \
        if ( loglevel_argument == "test_suite" )                                                                           \
            BOOST_TEST_MESSAGE( "    Entering common test case \"" + lama_common_testcase_method + "\" " );

/*
 * @brief HelperMakro LAMA_COMMON_TEST_CASE_END()
 *
 * This makro closes this testmethod.
 *
 */

#define LAMA_COMMON_TEST_CASE_TEMPLATE_END();                                                                          \
    if ( loglevel_argument == "test_suite" )                                                                           \
        BOOST_TEST_MESSAGE( "    Leaving common test case \"" + lama_common_testcase_method  );                        \
    }

/*
 * @brief HelperMakro LAMA_COMMON_TEST_CASE_RUNNER( classname )
 *
 * This makro creates the runTests()-Method to invoke all common testmethods.
 *
 */

#define LAMA_COMMON_TEST_CASE_RUNNER_TEMPLATE( classname );                                                            \
    template<typename StorageType>                                                                                         \
    void classname<StorageType>::runTests()

#if defined(LAMA_ASSERT_LEVEL_OFF)

// Do not check for exception if ASSERTions are switched off

#define LAMA_CHECK_THROW( stmt, exception )

#else

// use Boost macro for throwing exception

#define LAMA_CHECK_THROW( stmt, exception )    \
    BOOST_CHECK_THROW( stmt, exception )

#endif

#endif // LAMA_TESTMACROS_HPP_
