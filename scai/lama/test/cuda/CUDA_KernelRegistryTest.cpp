/**
 * @file CUDA_KernelRegistryTest.cpp
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
 * @brief Test to verify that 'relevant' CUDA kernel implementations have been registered
 * @author: Thomas Brandes
 * @date 27.10.2015
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/LAMAKernel.hpp>
#include <scai/lama/UtilKernelTrait.hpp>

#include <scai/hmemo/Context.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDA_KernelRegistry );

SCAI_LOG_DEF_LOGGER( logger, "Test.CUDA_KernelRegistry" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getKernelTest, ValueType, test_types )
{
    // This test checks that CUDA routines have been registered correctly and are accessible

    ContextPtr context = Context::getContextPtr( scai::common::context::CUDA );

    {
        LAMAKernel<UtilKernelTrait::scale<ValueType> > scale;
        SCAI_LOG_INFO( logger, "scale = " << scale.printIt() )
        BOOST_CHECK( scale[context] != NULL );
    }
 
    {
        LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;
        SCAI_LOG_INFO( logger, "CSR::normalGEMV = " << normalGEMV.printIt() )
        BOOST_CHECK( normalGEMV[context] != NULL );
    }
 
    {
        LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;
        SCAI_LOG_INFO( logger, "ELL::normalGEMV = " << normalGEMV.printIt() )
        BOOST_CHECK( normalGEMV[context] != NULL );
    }
 
    {
        LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;
        SCAI_LOG_INFO( logger, "JDS::normalGEMV = " << normalGEMV.printIt() )
        BOOST_CHECK( normalGEMV[context] != NULL );
    }
 
    {
        LAMAKernel<COOKernelTrait::normalGEMV<ValueType> > normalGEMV;
        SCAI_LOG_INFO( logger, "COO::normalGEMV = " << normalGEMV.printIt() )
        BOOST_CHECK( normalGEMV[context] != NULL );
    }
 
    {
        LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;
        BOOST_CHECK( normalGEMV[context] != NULL );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

