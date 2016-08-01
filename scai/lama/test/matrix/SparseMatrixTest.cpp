/**
 * @file test/matrix/SparseMatrixTest.cpp
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
 * @brief Test routines for specific methods/constructors of SparseMatrix
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/test/storage/TestStorages.hpp>

#include <scai/lama/matrix/SparseMatrix.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SparseMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef RealType ValueType;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( LocalConstructorTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    TypedStorages<ValueType> storages( ctx );

    for ( size_t i = 0; i < storages.size(); ++i )
    {
        common::shared_ptr<MatrixStorage<ValueType> > storage = storages[i];

        setDenseData( *storage );

        dmemo::DistributionPtr rowDist( new dmemo::BlockDistribution( storage->getNumRows(), comm ) );

        SparseMatrix<ValueType> matrix( storage, rowDist );

        BOOST_CHECK_EQUAL( matrix.getContextPtr(), storage->getContextPtr() );

        BOOST_CHECK( matrix.isConsistent() );

        BOOST_CHECK_EQUAL( matrix.getNumRows(), rowDist->getGlobalSize() );
        BOOST_CHECK_EQUAL( matrix.getNumColumns(), storage->getNumColumns() );
        BOOST_CHECK_EQUAL( matrix.getMatrixKind(), Matrix::SPARSE );

        // Note: storage has been modified (as shared pointer)

        BOOST_CHECK_EQUAL( storage->getNumRows(), rowDist->getLocalSize() );
    }
}

/* ------------------------------------------------------------------------- */

/*
BOOST_AUTO_TEST_CASE( FullConstructorTest )
{
    TypedStorages<ValueType> storages;

    for ( size_t i = 0; i < storages.size(); ++i )
    {
        storage = storages[i];

        matrix1.assign( storage );

        SparseMatrix matrix( local, halo, halo, rowDist, colDist );

    }
}
*/

/* ------------------------------------------------------------------------- */

/*
BOOST_AUTO_TEST_CASE( TransposeConstructorTest )
{
    TypedStorages<ValueType> storages;

    for ( size_t i = 0; i < storages.size(); ++i )
    {
        storage = storages[i];

        matrix1.assign( storage );

        SparseMatrix matrix( matrix1, transposeFlag );

        SparseMatrix matrix( matrix1, dist, dist );

    }
}
*/

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
