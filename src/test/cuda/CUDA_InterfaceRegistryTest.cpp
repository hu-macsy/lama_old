/**
 * @file CUDA_InterfaceRegistry.cpp
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
 * @brief Test of access to the LAMAInterface for CUDA device.
 * @author: Thomas Brandes
 * @date 13.04.2013
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/ContextFactory.hpp>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<double,float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDA_InterfaceRegistry );

LAMA_LOG_DEF_LOGGER( logger, "Test.CUDA_InterfaceRegistry" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( hasInterfaceTest )
{
    // Test will take the default CUDA device

    BOOST_CHECK( LAMAInterfaceRegistry::getRegistry().hasInterface( Context::CUDA ) );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getInterfaceTest, ValueType, test_types )
{
    // This test checks that CUDA routines have been registered correctly and are accessible

    ContextPtr context = ContextFactory::getContext( Context::CUDA );

    LAMA_INTERFACE_FN_T( scale, context, Utils, Transform, ValueType );

    BOOST_CHECK( scale ); 

    // normalGEMV is matrix-vector multiplication and is supported for each sparse matrix type

    {
        LAMA_INTERFACE_FN_T( normalGEMV, context, CSRUtils, Mult, ValueType );
        BOOST_CHECK( normalGEMV );
    }

    {
        LAMA_INTERFACE_FN_T( normalGEMV, context, ELLUtils, Mult, ValueType );
        BOOST_CHECK( normalGEMV );
    }

    {
        LAMA_INTERFACE_FN_T( normalGEMV, context, JDSUtils, Mult, ValueType );
        BOOST_CHECK( normalGEMV );
    }

    {
        LAMA_INTERFACE_FN_T( normalGEMV, context, DIAUtils, Mult, ValueType );
        BOOST_CHECK( normalGEMV );
    }

    {
        LAMA_INTERFACE_FN_T( normalGEMV, context, COOUtils, Mult, ValueType );
        BOOST_CHECK( normalGEMV );
    }

    // gemv is matrix-vector multiplication from BLAS2 used for dense matrix times vector

    {
        LAMA_INTERFACE_FN_T( gemv, context, BLAS, BLAS2, ValueType );
        BOOST_CHECK( gemv );
    }

    // normalGEVM is vector-matrix multiplication and is supported for each sparse matrix type

    {
        LAMA_INTERFACE_FN_T( normalGEVM, context, CSRUtils, Mult, ValueType );
        BOOST_CHECK( normalGEVM );
    }

    {
        LAMA_INTERFACE_FN_T( normalGEVM, context, ELLUtils, Mult, ValueType );
        BOOST_CHECK( normalGEVM );
    }

    {
        LAMA_INTERFACE_FN_T( normalGEVM, context, JDSUtils, Mult, ValueType );
        BOOST_CHECK( normalGEVM );
    }

    {
        LAMA_INTERFACE_FN_T( normalGEVM, context, DIAUtils, Mult, ValueType );
        BOOST_CHECK( normalGEVM );
    }

    {
        LAMA_INTERFACE_FN_T( normalGEVM, context, COOUtils, Mult, ValueType );
        BOOST_CHECK( normalGEVM );
    }

    // BLAS.LAPACK<ValueType>.getrf is not availabe for CUDA

    BOOST_CHECK_THROW( { LAMA_INTERFACE_FN_T( getrf, context, BLAS, LAPACK, ValueType ) }, Exception )

    // BLAS.LAPACK<ValueType>.getrf is not availabe for CUDA but for Host

    LAMA_INTERFACE_FN_DEFAULT_T( getrf, context, BLAS, LAPACK, ValueType )

    BOOST_CHECK( getrf );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

