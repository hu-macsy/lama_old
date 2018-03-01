/**
 * @file ParMetisTest.cpp
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
 * @brief Tests for working ParMetis and sound interface to LAMA
 * @author Thomas Brandes
 * @date 24.08.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <parmetis.h>

#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

using namespace scai;
using namespace dmemo;

BOOST_AUTO_TEST_SUITE( ParMetisTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.ParMetisTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( simpleTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    int nPart = comm->getSize();
    int rank  = comm->getRank();

    // Example Graph 
    //
    //       0 --- 1 --- 2 --- 3 --- 4 
    //       |     |     |     |     |
    //       5 --- 6 --- 7 --- 8 --- 9 
    //       |     |     |     |     |
    //      10 -- 11 -- 12 -- 13 -- 14 

    int ia[] = { 0, 2, 5, 8, 11, 13, 16, 20, 24, 28, 31, 33, 36, 39, 42, 44 };
    int ja[] = { 1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2,
                 6, 8, 12, 3, 7, 9, 13, 4, 8, 14, 5, 11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};

    int allNodes = sizeof( ia ) / sizeof( int ) - 1;
    int numValues = sizeof( ja ) / sizeof( int );

    BOOST_CHECK_EQUAL( allNodes[ia], numValues );

    // set the vertex distribution

    std::unique_ptr<int[]> vtxdist( new int[ nPart + 1 ] );

    // scan of nNodes among all processors

    vtxdist[0] = 0;

    IndexType lb;
    IndexType ub;

    for ( int ip = 0; ip < nPart; ++ip )
    {
        BlockDistribution::getLocalRange( lb, ub, allNodes, ip, nPart );
        BOOST_CHECK_EQUAL( lb, vtxdist[ip] );
        vtxdist[ip+1] = ub;
    }

    BOOST_CHECK_EQUAL( vtxdist[nPart], allNodes );

    // now build the local CSR graph of it as required by ParMetis

    int nNodes = vtxdist[rank+1] - vtxdist[rank];

    common::scoped_array<int> adjIA( new int[ nNodes + 1 ] );
    common::scoped_array<int> vwgt( new int[ nNodes ] );

    // build local IA array

    int nEdges = 0;

    adjIA[ 0 ] = 0;

    for ( int i = 0; i < nNodes; ++i )
    {
        nEdges += ia[ vtxdist[rank] + i + 1] - ia[vtxdist[rank] + i];
        adjIA[i+1] = nEdges;
        vwgt[i] = nEdges;
    }

    SCAI_LOG_INFO( logger, *comm << " has " << nNodes << " nodes and " << nEdges << " edges"
                             << " of total " << allNodes << " nodes and " << numValues << " edges" )

    common::scoped_array<int> adjJA( new int[ nEdges ] );

    for ( int i = 0; i < nNodes; ++i )
    {
        for ( int jj = adjIA[i]; jj < adjIA[i+1]; ++jj )
        {
            adjJA[jj] = ja[ ia[ vtxdist[rank] ] + jj ];
        }
    }
    idx_t ncon = 1;

    idx_t options[4];
    options[0] = 0;

    common::scoped_array<idx_t> partition( new idx_t [nNodes] );

    for ( int i = 0; i < nNodes; ++i )
    {
        partition[i] = -1;
    }

    idx_t wgtflag = 2; // weight on the vertices only
    idx_t numflag = 0; // where indexing starts

    MPI_Comm mpiComm = MPI_COMM_WORLD;  

    real_t itr = 1000.0;
    idx_t edgecut = 0;

    common::scoped_array<real_t> tpwgts( new real_t[ nPart ] );

    for ( int ip = 0; ip < nPart; ++ip )
    {
        tpwgts[ip] = real_t( 1 ) / real_t( nPart );
    }

    common::scoped_array<real_t> ubvec( new real_t[ ncon ] );

    for ( int i = 0; i < ncon; ++i )
    {
        ubvec[i] = 2.0;
    }

    int res = ParMETIS_V3_AdaptiveRepart (
        vtxdist.get(),
        adjIA.get(),
        adjJA.get(),
        vwgt.get(),
        vwgt.get(),
        NULL, // no weight for edges
        &wgtflag,
        &numflag,
        &ncon,
        &nPart,
        tpwgts.get(),
        ubvec.get(),
        &itr,
        options,
        &edgecut,
        partition.get(),
        &mpiComm );

    SCAI_LOG_INFO( logger, *comm << " has finished, res = " << res << ", edgecut = " << edgecut )
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
