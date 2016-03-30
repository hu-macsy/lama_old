/**
 * @file Matrices.hpp
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
 * @brief Vector with all matrices, one for each supported matrix storage format/type
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
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
