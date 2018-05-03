/**
 * @file CSRGraph.hpp
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
 * @brief Class that represents adjacency matrix as CSR graph.
 * @author Thomas Brandes
 * @date 02.01.2018
 */

#pragma once

#include <scai/logging.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/dmemo/Distribution.hpp>
#include <scai/logging.hpp>

namespace scai
{

namespace partitioning
{

/**
 *  Distributed CSR format as used by partitioning tools like ParMetis
 *
 *  @tparam IdxType data type used for index representaton (might be different to scai::IndexType)
 */
template<typename IdxType>
class COMMON_DLL_IMPORTEXPORT CSRGraph
{

public:

    /** Constructor for a partitioning.  */

    CSRGraph( const lama::_Matrix& matrix, bool isSingle = false );

    /** Destructor of partitionng */

    virtual ~CSRGraph();

    const hmemo::HArray<IdxType>& distOffsets() const 
    { 
        return mDistOffsets; 
    } 

    const hmemo::HArray<IdxType>& ia() const 
    { 
        return mIA; 
    } 

    const hmemo::HArray<IdxType>& ja() const 
    { 
        return mJA;
    } 

    const hmemo::HArray<IdxType>& weights() const
    {
        return mVertexWeight;
    }

private:

    hmemo::HArray<IdxType> mDistOffsets;    // vertex distribution, scan of local sizes of vertices
    hmemo::HArray<IdxType> mIA;             // offset array
    hmemo::HArray<IdxType> mJA;             // column indexes
    hmemo::HArray<IdxType> mVertexWeight;   // weights of local vertices

    /** build offset array for local sizes of processors, same as built for general block distribution */

    static void getDistributionOffsets( hmemo::HArray<IdxType>& offsets, const dmemo::Distribution& dist );

    void joinCSR(
        hmemo::HArray<IdxType>& outIA,
        hmemo::HArray<IdxType>& outJA,
        const hmemo::HArray<IndexType>& localIA,
        const hmemo::HArray<IndexType>& localJA,
        const hmemo::HArray<IndexType>& haloIA,
        const hmemo::HArray<IndexType>& haloJA );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename ValueType>
    void buildByCSRSparseMatrix( const lama::CSRSparseMatrix<ValueType>& matrix, bool isSingle );
};

} /* end namespace partitioning */

} /* end namespace scai */
