/**
 * @file lamaTest.cpp
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
 * @brief Contains the implementation of the class lamaTest.
 * @author: Alexander Büchel
 * @date 31.01.2012
 * @since 1.0.0
 **/

#ifndef BOOST_TEST_DYN_LINK
    #define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE lamaTest

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <boost/regex.hpp>

/* Including all headerfiles of testcases */
#include <scai/lama/test/NormTest.hpp>
#include <scai/lama/test/MatrixStorageTest.hpp>
#include <scai/lama/test/SparseMatrixTest.hpp>

#include <scai/common/Settings.hpp>
#include <scai/hmemo/test/ContextFix.hpp>

#include <list>
#include <string>
#include <map>

#include <cstdio>

/* Use global fixture for context and define the static variable.
 * Avoids expensive calls of reserve/release routines of the Context for each test.
 */

BOOST_GLOBAL_FIXTURE( ContextFix )

scai::hmemo::ContextPtr ContextFix::testContext;

/** The init function returns true if it can get the specified context. */

bool init_function()
{
    try
    {
        scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();
        return true;
    }
    catch ( scai::common::Exception& ex )
    {
        std::cerr << "Could not get context for test: " << ex.what() << std::endl;
        return false;
    }
}

bool base_test_case = false;
std::string testcase;
//std::string loglevel_argument = "";

/*
 * @brief HelperMakro TC_REG( test_classes, test_methods )
 *
 * This makro creates the entries for the registration of base test cases.
 *
 * @param test_classes     String-array of testclasses
 * @param test_methods     String-array of testmethods
 */

#define TC_REG( test_classes, test_methods )                                                         \
    { for ( int i = 0; i < (int)( sizeof(test_classes) / sizeof(test_classes[0]) ); i++ )            \
            for ( int j = 0; j < (int)( sizeof(test_methods) / sizeof(test_methods[0]) ); j++ )      \
                inherited_test_cases.push_back( test_classes[i]+"/"+test_methods[j] ); }

int main( int argc, char* argv[] )
{
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    if ( argc > 1 )
    {
        std::string runtest_argument;
        std::string testsuite;
        std::list<std::string>::iterator Iter;
        std::list<std::string> inherited_test_cases;
        std::map<std::string, std::string> runtime_arguments;
        /* Registration of common Testclasses and their testmethods */
        TC_REG( normtestclasses, normtestmethods );
        TC_REG( sparsematrixtestclasses, sparsematrixtestmethods );
        TC_REG( storagetestclasses, storagetestmethods );

        /* Parse command line arguments */
        for ( int i = 0; i < argc; i++ )
        {
            std::string s = argv[i];
            int pos = s.find( '=' );

            if ( pos != -1 )
                runtime_arguments.insert(
                    std::make_pair( s.substr( 0, pos ),
                                    s.substr( pos + 1, s.length() - s.substr( 0, pos ).length() - 1 ) ) );
        }

        /* Find run_test in runtime parameters */
        std::map<std::string, std::string>::iterator iterator = runtime_arguments.find( "--run_test" );

        if ( iterator != runtime_arguments.end() )
        {
            runtest_argument = iterator->second;
        }

        //if a regular expression is used in a base test:

        int posasterix = runtest_argument.find( '*' );

        if ( posasterix != -1 )
        {
            std::list<std::string> matchedtestcases;
            std::string runtest_argument_regex = boost::regex_replace( runtest_argument, boost::regex( "[*]" ),
                                                 "[a-z]*" );

            for ( Iter = inherited_test_cases.begin(); Iter != inherited_test_cases.end(); Iter++ )
            {
                boost::regex re;
                re.assign( runtest_argument_regex, boost::regex_constants::icase );

                if ( boost::regex_match( *Iter, re ) )
                {
                    matchedtestcases.push_back( *Iter );
                }
            }

            //if there is no matched test case, we run main from boost automatically
            if ( matchedtestcases.size() == 0 )
            {
                // do nothing here, because a base test case will be not detected,
                // so main from boost will be invoked with the given expression automatically
            }
            //if there is just 1 matched test case, we map this test case to the full name
            else if ( matchedtestcases.size() == 1 )
            {
                for ( Iter = matchedtestcases.begin(); Iter != matchedtestcases.end(); Iter++ )
                {
                    runtest_argument = *Iter;
                }
            }
            else
                // if there are 2 or more matched test cases, we will invoke all test cases from the given testsuite.
                // We assume here, that there is just one individual test suite in the matched test cases.
            {
                for ( Iter = matchedtestcases.begin(); Iter != matchedtestcases.end(); Iter++ )
                {
                    runtest_argument = *Iter;
                }

                int pos = runtest_argument.find( '/' );
                testsuite = runtest_argument.substr( 0, pos );
                runtest_argument = testsuite + "/*";

                for ( int i = 0; i < argc; i++ )
                {
                    std::string valueofarg = argv[i];
                    int pos = valueofarg.find( '=' );
                    std::string category = valueofarg.substr( 0, pos );

                    if ( category == "--run_test" )
                    {
                        std::string newarg = "";
                        newarg = category + "=" + testsuite + "/*";
                        char* help = new char[newarg.size() + 1];
                        newarg.copy( help, newarg.size() );
                        argv[i] = help;
                        break;
                    }
                }
            }

            //TODO: Think about, what will happen, if there are different test suites ...
        }

        /* If this argument is directly a base test case or the runtime parameter was mapped to a full name,
         * the next step is to save this to change the run_test parameter to "commonTestCase" and to adjust
         * argv. The variable testcase contains the concrete testcasename.
         */

        /* Check, if argument of run_test is a base test case */
        for ( Iter = inherited_test_cases.begin(); Iter != inherited_test_cases.end(); Iter++ )
            if ( runtest_argument == *Iter )
            {
                base_test_case = true;
            }

        if ( base_test_case )
        {
            int pos = runtest_argument.find( '/' );
            testsuite = runtest_argument.substr( 0, pos );
            testcase = runtest_argument.substr(
                           pos + 1, runtest_argument.length() - runtest_argument.substr( 0, pos ).length() - 1 );
            runtime_arguments["--run_test"] = testsuite + "/commonTestCases*";

            for ( int i = 0; i < argc; i++ )
            {
                std::string value = argv[i];
                int pos = value.find( '=' );
                std::string category = value.substr( 0, pos );

                if ( category == "--run_test" )
                {
                    std::string string = category + "=" + runtime_arguments["--run_test"];
                    char* help = new char[string.size() + 1];
                    string.copy( help, string.size() + 1 );
                    argv[i] = help;
                    break;
                }
            }
        }
    }

    /* Call main() from boost with new arguments */
    /* it is just allowed to invoke this method once */

    return boost::unit_test::unit_test_main( &init_function, argc, argv );
}
