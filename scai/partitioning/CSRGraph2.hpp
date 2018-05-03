/**
 * @file CSRGraph2.hpp
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
 * @brief Class that represents adjacency matrix of a bipartite CSR graph, i.e. vertices for each row and column
 * @author Thomas Brandes
 * @date 02.01.2018
 */

#pragma once

#include <scai/logging.hpp>

#include <scai/lama/matrix/_Matrix.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/dmemo/Distribution.hpp>
#include <scai/logging.hpp>

namespace scai
{

namespace partitioning
{

/**
 *  Bipartite CSR graph, vertices for each row and each column 
 *
 *  @tparam IdxType data type used for index representaton (might be different to scai::IndexType)
 */
template<typename IdxType>
class COMMON_DLL_IMPORTEXPORT CSRGraph2
{

public:

    /** Constructor for a partitioning.  */

    CSRGraph2( const lama::_Matrix& matrix );

    /** Destructor of partitionng */

    virtual ~CSRGraph2();

    const hmemo::HArray<IdxType>& ia() const 
    { 
        return mIA;
    } 

    const hmemo::HArray<IdxType>& ja() const 
    { 
        return mJA; 
    } 

    const hmemo::HArray<IdxType>& vertexWeights() const 
    { 
        return mVertexWeights;
    } 

private:

    IdxType mNumConstraints = 3;    // number of weights used for each vertex

    hmemo::HArray<IdxType> mIA;             // offset array
    hmemo::HArray<IdxType> mJA;             // column indexes
    hmemo::HArray<IdxType> mVertexWeights;   // weights of local vertices

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace partitioning */

} /* end namespace scai */
