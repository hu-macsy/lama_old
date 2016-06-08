/**
 * @file test/matrix/CSRSparseMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Contains only test specific for the CSR Sparse matrix
 * @author Thomas Brandes
 * @date 24.03.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/test/TestMacros.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CSRSparseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRSparseMatrixTest" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( defaultConstructorTest, ValueType, scai_arithmetic_test_types )
{
    CSRSparseMatrix<ValueType> matrix;

    // check zero sizes

    BOOST_CHECK_EQUAL( 0, matrix.getNumRows() );
    BOOST_CHECK_EQUAL( 0, matrix.getNumColumns() );

    // check correct format / type

    BOOST_CHECK_EQUAL( common::TypeTraits<ValueType>::stype, matrix.getValueType() );
    BOOST_CHECK_EQUAL( Matrix::CSR, matrix.getFormat() );

    const CSRStorage<ValueType>& local = matrix.getLocalStorage();
    const CSRStorage<ValueType>& halo = matrix.getHaloStorage();

    BOOST_CHECK_EQUAL( local.getNumRows(), halo.getNumRows() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
