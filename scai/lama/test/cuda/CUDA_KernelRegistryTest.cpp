/**
 * @file CUDA_KernelRegistryTest.cpp
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
 * @brief Test to verify that 'relevant' CUDA kernel implementations have been registered
 * @author Thomas Brandes
 * @date 27.10.2015
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/ELLKernelTrait.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>
#include <scai/sparsekernel/DIAKernelTrait.hpp>
#include <scai/sparsekernel/JDSKernelTrait.hpp>

#include <scai/hmemo/Context.hpp>

//using namespace scai::lama;
using namespace scai::sparsekernel;
using namespace scai::utilskernel;
using namespace scai::hmemo;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDA_KernelRegistry );

SCAI_LOG_DEF_LOGGER( logger, "Test.CUDA_KernelRegistry" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getKernelTest, ValueType, test_types )
{
    // This test checks that CUDA routines have been registered correctly and are accessible

    ContextPtr context = Context::getContextPtr( Context::CUDA );

    {
        LAMAKernel<UtilKernelTrait::setVal<ValueType, ValueType> > setVal;
        SCAI_LOG_INFO( logger, "setVal = " << setVal.printIt() )
        BOOST_CHECK( setVal[context] != NULL );
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

