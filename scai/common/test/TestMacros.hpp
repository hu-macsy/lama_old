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

// TODO: no dependencies to project ahead
#include <scai/hmemo/Context.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/kregistry/exception/KernelRegistryException.hpp>

// boost
#include <boost/assign/list_of.hpp>
#include <boost/preprocessor.hpp>

// std
#include <string>
#include <map>

extern std::string loglevel_argument;

/**
 * @brief eps() returns the desired precision of calculations with ValueType in tests.
 *
 * @return the desired precision for calculations with ValueType in tests
 */
template<typename ValueType>
inline ValueType eps();

/**
 * @brief eps<double>() returns the desired precision of calculations in single
 *        precision.
 *
 * @return the desired precision for calculations with single precision in tests
 */
template<>
inline float eps<float>()
{
    return 1E-3f;
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
    return 1E-5;
}

template<>
inline long double eps<long double>()
{
    return 1E-8L;
}

template<>
inline ComplexFloat eps<ComplexFloat>()
{
    return ComplexFloat(1E-3f);
}

template<>
inline ComplexDouble eps<ComplexDouble>()
{
    return ComplexDouble(1E-5);
}

template<>
inline ComplexLongDouble eps<ComplexLongDouble>()
{
    return ComplexLongDouble(1E-8L);
}


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

/**
 * @brief getEnvContext() returns the current context from enviroment variable LAMA_TEST_CONTEXT
 * If this variable is not set or contains a context, which is not available, all contexts will be used.
 *
 * @return the current context from environment variable LAMA_TEST_CONTEXT as a string
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
        *p = static_cast<char>( toupper( *p ) );
    }

    if ( context == "HOST" )
    {
        context = "Host";    // set back
    }
    else if ( context == "CUDA" || context == "OPENCL" || context == "MIC" )
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

inline scai::common::context::ContextType mapEnvContexttoContextType( std::string contextname )
{
    using namespace scai::common::context;

           ContextType myContext;
    std::map<std::string, ContextType> contextmap =
        boost::assign::map_list_of ( "Host", Host )
        ( "CUDA", CUDA )
        ( "OPENCL", OpenCL )
        ( "MIC", MIC )
        ( "MaxContext", MaxContext );
    myContext = contextmap[contextname];
    return myContext;
}

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

#define SCAI_CHECK_CLOSE( x, y, tolerance )                         \
    {                                                               \
        Scalar xScalar = Scalar( x );                               \
        Scalar yScalar = Scalar( y );                               \
        ComplexDouble xVal = xScalar.getValue<ComplexDouble>();     \
        ComplexDouble yVal = yScalar.getValue<ComplexDouble>();     \
        BOOST_CHECK_CLOSE( xVal.real(), yVal.real(), tolerance );   \
        BOOST_CHECK_CLOSE( xVal.imag(), yVal.imag(), tolerance );   \
    }

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

#define SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps )                                                                   \
    ValueType xHelper = (x).getValue<ValueType>();                                                                     \
    BOOST_CHECK_SMALL( xHelper, static_cast<ValueType>( eps ) );

/*
 * @brief HelperMacro SCAI_CHECK_SCALAR_SMALL_EPS( x, ValueType )
 *
 * Same as SCAI_CHECK_SCALAR_SMALL but with default eps value.
 *
 * @param x             Scalar
 * @param ValueType     type of Scalar to be used for test
 */

#define SCAI_CHECK_SCALAR_SMALL_EPS( x, ValueType )                                                                    \
    SCAI_CHECK_SCALAR_SMALL( x, ValueType, eps<ValueType> () )

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
    { scai::common::Printable* p = dynamic_cast<scai::common::Printable*>( &printable );                                                           \
        std::stringstream mStream;                                                                                     \
        p->writeAt( mStream );                                                                                         \
        std::string mString = mStream.str();                                                                           \
        BOOST_CHECK( mString.length() > 0 ); }

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

#define LAMA_WRITEAT_PTR_TEST( printable )                                                                             \
    { std::stringstream mStream;                                                                                       \
        printable->writeAt( mStream );                                                                                 \
        std::string mString = mStream.str();                                                                           \
        BOOST_CHECK( mString.length() > 0 ); }

/*
 * @brief HelperMacro CONTEXTLOOP()
 *
 * This macro creates a loop iterating throw all existing contexts.
 * All contexts, which are not adressable, will be skipped.
 *
 * If the environment variable LAMA_TEST_CONTEXT or the runtime parameter --context
 * is set to a specific context, just this context will be used.
 */

#define CONTEXTLOOP()                                                                                                  \
    std::list<scai::common::context::ContextType> listofcontexts;                                                      \
    std::list<scai::common::context::ContextType>::iterator Iter;                                                      \
    std::string contexttype;                                                                                           \
    contexttype = getEnvContext();                                                                                     \
    if ( contexttype == "*" )                                                                                          \
    {                                                                                                                  \
        SCAI_LOG_INFO( logger, "LAMA_TEST_CONTEXT is not set or has value '*', so all available contexts will be used." );  \
        for ( int i = 0; i < scai::common::context::MaxContext; ++i )                                                  \
        {                                                                                                              \
            scai::common::context::ContextType ctx = static_cast<scai::common::context::ContextType>( i + 1 );         \
            if ( Context::hasContext( ctx ) )                                                                          \
            {                                                                                                          \
                listofcontexts.push_back( ctx );                                                                       \
                SCAI_LOG_DEBUG( logger, "Context " << ctx << " is available");                                           \
            }                                                                                                          \
            else                                                                                                       \
                SCAI_LOG_INFO( logger, "The following context will be skipped, because it is not available: " << ctx );\
        }                                                                                                              \
    } else {                                                                                                           \
        listofcontexts.push_back( mapEnvContexttoContextType( contexttype ) );                                         \
        SCAI_LOG_INFO( logger, "Environment variable LAMA_TEST_CONTEXT contains context = " << getEnvContext() );      \
    }                                                                                                                  \
    for ( Iter = listofcontexts.begin(); Iter != listofcontexts.end(); Iter++ )

/*
 * @brief HelperMacro GETCONTEXT( loc )
 *
 * This macro creates the adressable context. It is only usable in a CONTEXTLOOP().
 *
 * @param loc          loc is the variable, which contains the context.
 */

#define GETCONTEXT( loc )                                                                                              \
    ContextPtr loc;                                                                                                    \
    loc = Context::getContextPtr( *Iter );

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
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context );                                                                                     \
            if ( loglevel_argument == "test_suite" )                                                                   \
                SCAI_LOG_INFO( logger, "    Entering context: " << context->getType() );                               \
                scai::lama::classname::name( context );                                                                \
        }                                                                                                              \
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

// we need the double trick, otherwise we just see ARITHMETIC_HOST_TYPE_##I in output

#define STR1( x ) #x
#define STR( x ) STR1( x )

/** This macro runs method<ValueType>( context ) where ValueType is given by the I-th arithmetic Host Type 
 *  (skips the run for long double types on CUDA as not supported there).
 *
 *  KernelRegistryException is caught with a correpsonding warn message on logger
 */

#define LAMA_RUN_TEST(z, I, method )                                                                            \
    try                                                                                                         \
    {                                                                                                           \
        if ( context->getType() == scai::common::context::CUDA )                                                \
        {                                                                                                       \
            switch( scai::common::getScalarType<ARITHMETIC_HOST_TYPE_##I>() )                                   \
            {                                                                                                   \
                case scai::common::scalar::LONG_DOUBLE:                                                         \
                case scai::common::scalar::LONG_DOUBLE_COMPLEX:                                                 \
                    continue;                                                                                   \
                default:                                                                                        \
                     ;                                                                                          \
             }                                                                                                  \
        }                                                                                                       \
        method<ARITHMETIC_HOST_TYPE_##I>( context );                                                            \
    }                                                                                                           \
    catch ( scai::kregistry::KernelRegistryException& )                                                         \
    {                                                                                                           \
        SCAI_LOG_WARN( logger, #method << "<" << STR( ARITHMETIC_HOST_TYPE_##I ) << "> cannot run on "          \
                       << context->getType() << ", corresponding function not implemented yet."        )        \
        return;                                                                                                 \
    }                                                                                                           \

/*
 * @brief HelperMacro LAMA_AUTO_TEST_CASE_CT( name, classname )
 *
 * This macro creates a boost test auto case, which uses all possible contexts.
 * The test case name is based on the name of the given test method.
 *
 * @param name          name of test method, which will invoke.
 * @param classname     name of the given test class.
*/
#define LAMA_AUTO_TEST_CASE_CT( name, classname, namespacename )                                     \
                                                                                                     \
    BOOST_AUTO_TEST_CASE( name )                                                                     \
    {                                                                                                \
        CONTEXTLOOP()                                                                                \
        {                                                                                            \
            GETCONTEXT( context )                                                                    \
            if ( loglevel_argument == "test_suite" )                                                 \
            {                                                                                        \
                SCAI_LOG_INFO( logger, "    Entering context: " << context->getType() );             \
            }                                                                                        \
            const std::string lama_name = #name;                                                     \
            const std::string lama_classname = #classname;                                           \
            BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_RUN_TEST, namespacename::classname::name )    \
        }                                                                                            \
    }

#define LAMA_RUN_TESTL(z, I, method )                                                                  \
    try                                                                                                \
    {                                                                                                  \
        method<ARITHMETIC_HOST_TYPE_##I>( context, logger );                                           \
    }                                                                                                  \
    catch ( scai::kregistry::KernelRegistryException& )                                                \
    {                                                                                                  \
        SCAI_LOG_WARN( logger, #method << "<" << STR( ARITHMETIC_HOST_TYPE_##I ) << "> cannot run on " \
                       << context->getType() << ", corresponding function not implemented yet." );     \
        return;                                                                                        \
    }                                                                                                  \

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
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context )                                                                                      \
            if ( loglevel_argument == "test_suite" )                                                                   \
            {                                                                                                          \
                SCAI_LOG_INFO( logger, "    Entering context: " << context->getType() );                               \
            }                                                                                                          \
            BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE__CNT, LAMA_RUN_TESTL, scai::lama::classname::name )                 \
        }                                                                                                              \
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
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context );                                                                                     \
            if ( loglevel_argument == "test_suite" )                                                                   \
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
        }                                                                                                              \
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
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context );                                                                                     \
            if ( loglevel_argument == "test_suite" )                                                                   \
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
        }                                                                                                              \
    }

#define LAMA_AUTO_TEST_CASE_CDUMMY( name, classname )                                                                 \
    BOOST_AUTO_TEST_CASE( name )                                                                                       \
    {                                                                                                                  \
        CONTEXTLOOP()                                                                                                  \
        {                                                                                                              \
            GETCONTEXT( context );                                                                                     \
            if ( loglevel_argument == "test_suite" )                                                                   \
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
        }                                                                                                              \
    }

    /*
     * @brief HelperMacro COMMONTESTCASEINVOKER( object_name, method_name )
     *
     * This macro represents the invoke of a common used test method.
     *
     * @param object_name     name of the object, which consists the test method.
     * @param method_name     name of method that will be called.
     */

#define COMMONTESTCASEINVOKER( object_name, method_name )                                                              \
    { if ( testcase == #method_name ) try {object_name.method_name();} catch(scai::common::UnsupportedException const& ex) { std::cout << "failed on " #method_name << std::endl; } }

    /*
     * @brief HelperMacro COMMONTESTCASEINVOKER_TEMPLATE( object_name, method_name, ValueType )
     *
     * This macro represents the invoke of a common used test method.
     *
     * @param object_name     name of the object, which consists the test method.
     * @param method_name     name of method that will be called.
     */

#define COMMONTESTCASEINVOKER_TEMPLATE( object_name, method_name, ValueType )                                          \
    { if ( testcase == #method_name ) object_name.method_name<ValueType>(); }

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
        if ( loglevel_argument == "test_suite" )                                                                       \
            BOOST_TEST_MESSAGE( "    Entering common test case \"" + lama_common_testcase_method + "\" " );            \

        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE_TM_END()
         *
         * This macro closes this testmethod.
         *
         */

#define LAMA_COMMON_TEST_CASE_TM_END();                                                                                \
    if ( loglevel_argument == "test_suite" )                                                                           \
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
        if ( loglevel_argument == "test_suite" )                                                                       \
            BOOST_TEST_MESSAGE( "    Entering common test case \"" + lama_common_testcase_method + "\" " );            \

        /*
         * @brief HelperMacro LAMA_COMMON_TEST_CASE_END()
         *
         * This macro closes this test method.
         *
         */

#define LAMA_COMMON_TEST_CASE_END();                                                                                   \
    if ( loglevel_argument == "test_suite" )                                                                           \
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
    extern std::string loglevel_argument;                                                                              \
    template<typename templatename>                                                                                    \
    void classname<templatename>::methodname()                                                                         \
    {                                                                                                                  \
        std::string lama_common_testcase_method = #methodname;                                                         \
        std::ostringstream omsg;                                                                                       \
        omsg << scai::common::getScalarType<templatename>();                                                                       \
        std::string lama_common_testcase_template = omsg.str();                                                        \
        if ( loglevel_argument == "test_suite" )                                                                       \
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

#define LAMA_COMMON_TEST_CASE_TEMPLATE_END();                                                                          \
    if ( loglevel_argument == "test_suite" )                                                                           \
        SCAI_LOG_INFO( logger, "    Leaving common test case \"" + lama_common_testcase_method + "\" " );                \
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

