/**
 * @file Matrices.hpp
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
 * @brief Vector with all matrices, one for each supported matrix storage format/type
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/common/shared_ptr.hpp>
#include <scai/common/TypeTraits.hpp>

#include <vector>

/** Class for a list of matrix pointers, one for each supported 
 *  matrix storage format and each supported arithmetic type.
 */

class Matrices : public std::vector<scai::lama::MatrixPtr> 
{

public:

    /** Constructor creates already the list with all matrix pointers. */

    Matrices( scai::hmemo::ContextPtr ctx = scai::hmemo::ContextPtr() ) 
    {
        using namespace scai::lama;

        std::vector<MatrixCreateKeyType> values;  //  all create values

        Matrix::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            scai::lama::MatrixPtr matrixPtr( scai::lama::Matrix::create( values[i] ) );

            if ( ctx )
            {
                matrixPtr->setContextPtr( ctx );
            }

            push_back( matrixPtr );
        }
    }

    Matrices( scai::common::scalar::ScalarType stype, scai::hmemo::ContextPtr ctx = scai::hmemo::ContextPtr() ) 
    {
        using namespace scai::lama;

        std::vector<MatrixCreateKeyType> values;  //  all create values

        Matrix::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            scai::lama::MatrixPtr matrixPtr( scai::lama::Matrix::create( values[i] ) );

            if ( values[i].second != stype )
            {
                continue;
            }

            if ( ctx )
            {
                matrixPtr->setContextPtr( ctx );
            }

            push_back( matrixPtr );
        }
    }

    // Destructor will free all matrices due to use of shared pointers
};
