/**
 * @file lama/test/DerivedSparseMatrixTest.cpp
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
 * @brief Contains the implementation of the class CSRSparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>

#include <scai/lama/test/SparseMatrixTest.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DerivedSparseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest.DerivedSparseMatrixTest" )

typedef SCAI_TEST_TYPE ValueType;

typedef boost::mpl::list<CSRSparseMatrix<ValueType>,
        COOSparseMatrix<ValueType>,
        DIASparseMatrix<ValueType>,
        DIASparseMatrix<ValueType>,
        ELLSparseMatrix<ValueType>,
        JDSSparseMatrix<ValueType> > test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( clear, MatrixType, test_types )
{
    MatrixType matrix;
    SparseMatrixTest<MatrixType> sparseMatrixtest( matrix );
    sparseMatrixtest.clearTest();
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( conversions, MatrixType, test_types )
{
    MatrixType matrix;
    SparseMatrixTest<MatrixType> sparseMatrixtest( matrix );
    sparseMatrixtest.testConversions();
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
