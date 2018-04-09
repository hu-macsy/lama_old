/**
 * @file test/matrix/StencilMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Test routines for specific methods/constructors of StencilMatrix
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/StencilMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StencilMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef DefaultReal ValueType;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Stencil1D3PTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType N1 = 98;

    common::Stencil1D<DefaultReal> stencil( 3 );
    common::Grid1D grid( N1 );

    StencilMatrix<DefaultReal> stencilMatrix( grid, stencil );

    CSRSparseMatrix<DefaultReal> csrMatrixPoisson;
    MatrixCreator::buildPoisson( csrMatrixPoisson, 1, 3, N1, 1, 1 );

    CSRSparseMatrix<DefaultReal> csrMatrixStencil( stencilMatrix );

    BOOST_REQUIRE_EQUAL( csrMatrixStencil.getNumRows(), csrMatrixStencil.getNumRows() );

    DefaultReal diff = csrMatrixStencil.maxDiffNorm( csrMatrixPoisson );

    BOOST_CHECK( diff < 1e-5 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Stencil2D5PTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType N1 = 10;
    const IndexType N2 = 4;

    common::Stencil2D<DefaultReal> stencil( 5 );
    common::Grid2D grid( N1, N2 );

    StencilMatrix<DefaultReal> stencilMatrix( grid, stencil );

    CSRSparseMatrix<DefaultReal> csrMatrixPoisson;
    MatrixCreator::buildPoisson( csrMatrixPoisson, 2, 5, N1, N2, 1 );

    CSRSparseMatrix<DefaultReal> csrMatrixStencil( stencilMatrix );

    BOOST_REQUIRE_EQUAL( csrMatrixStencil.getNumRows(), csrMatrixStencil.getNumRows() );

    DefaultReal diff = csrMatrixStencil.maxDiffNorm( csrMatrixPoisson );

    BOOST_CHECK( diff < 1e-5 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( Stencil3D27PTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType N1 = 10;
    const IndexType N2 = 4;
    const IndexType N3 = 7;

    common::Stencil3D<DefaultReal> stencil( 27 );
    common::Grid3D grid( N1, N2, N3 );

    StencilMatrix<DefaultReal> stencilMatrix( grid, stencil );

    CSRSparseMatrix<DefaultReal> csrMatrixPoisson;
    MatrixCreator::buildPoisson( csrMatrixPoisson, 3, 27, N1, N2, N3 );

    CSRSparseMatrix<DefaultReal> csrMatrixStencil( stencilMatrix );

    BOOST_REQUIRE_EQUAL( csrMatrixStencil.getNumRows(), csrMatrixStencil.getNumRows() );

    DefaultReal diff = csrMatrixStencil.maxDiffNorm( csrMatrixPoisson );

    BOOST_CHECK( diff < 1e-5 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CopyConstructorTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType N1 = 5;
    const IndexType N2 = 3;
    const IndexType N3 = 4;

    common::Stencil3D<DefaultReal> stencil( 27 );
    common::Grid3D grid( N1, N2, N3 );

    StencilMatrix<DefaultReal> stencilMatrix1( grid, stencil );
    SCAI_LOG_INFO( logger, "copy this stencil matrix: " << stencilMatrix1 )
    StencilMatrix<DefaultReal> stencilMatrix2( stencilMatrix1 );
    SCAI_LOG_INFO( logger, "copied stencil matrix = " << stencilMatrix2 )

    CSRSparseMatrix<DefaultReal> csrMatrixStencil1( stencilMatrix1 );
    SCAI_LOG_INFO( logger, "csr matrix 1 = " << csrMatrixStencil1 )
    CSRSparseMatrix<DefaultReal> csrMatrixStencil2( stencilMatrix2 );
    SCAI_LOG_INFO( logger, "csr matrix 2 = " << csrMatrixStencil2 )

    DefaultReal diff = csrMatrixStencil1.maxDiffNorm( csrMatrixStencil2 );

    BOOST_CHECK( diff < 1e-5 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( AssignTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType N1 = 5;
    const IndexType N2 = 3;
    const IndexType N3 = 4;

    common::Stencil3D<DefaultReal> stencil( 27 );
    common::Grid3D grid( N1, N2, N3 );

    StencilMatrix<DefaultReal> stencilMatrix1( grid, stencil );
    StencilMatrix<DefaultReal> stencilMatrix2;

    stencilMatrix2 = stencilMatrix1;

    CSRSparseMatrix<DefaultReal> csrMatrixStencil1( stencilMatrix1 );
    CSRSparseMatrix<DefaultReal> csrMatrixStencil2( stencilMatrix2 );

    DefaultReal diff = csrMatrixStencil1.maxDiffNorm( csrMatrixStencil2 );

    BOOST_CHECK( diff < 1e-5 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
