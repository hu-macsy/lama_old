/**
 * @file Matrices.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Vector with all matrices, one for each supported matrix storage format/type
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/common/TypeTraits.hpp>

#include <vector>

/** Class for a standard vector of matrix pointers, one for each supported
 *  matrix storage format and each supported arithmetic type.
 */
class _Matrices : public std::vector<std::unique_ptr<scai::lama::_Matrix>>
{
public:

    /** Constructor creates already the list with all matrix pointers. */

    _Matrices( scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr() )
    {
        using namespace scai::lama;
        std::vector<MatrixCreateKeyType> values;  //  all create values
        _Matrix::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            std::unique_ptr<scai::lama::_Matrix> matrixPtr( scai::lama::_Matrix::create( values[i] ) );

            matrixPtr->setContextPtr( ctx );

            push_back( std::move( matrixPtr ) );
        }
    }

    _Matrices( scai::common::ScalarType stype, scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr() )
    {
        using namespace scai::lama;
        std::vector<MatrixCreateKeyType> values;  //  all create values
        _Matrix::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            if ( values[i].second != stype )
            {
                continue;
            }

            std::unique_ptr<scai::lama::_Matrix> matrixPtr( scai::lama::_Matrix::create( values[i] ) );

            matrixPtr->setContextPtr( ctx );

            push_back( std::move( matrixPtr ) );
        }
    }

    // Destructor will free all matrices due to use of unique pointers
};

/**
 *  @brief Derived class for a vector of matrix pointers
 *
 *  The constructor of this class generates a vector/set of matrices of the
 *  given value type. These matrices are all initialized as zero matrices.
 */
template<typename ValueType>
class Matrices : public std::vector<typename scai::lama::MatrixPtr<ValueType> >
{

public:
    
    /** Constructor allocates a set of typed matrices, one for each supported format */
    
    Matrices( scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr() )
    {   
        using namespace scai;
        using namespace lama;

        common::ScalarType stype = common::TypeTraits<ValueType>::stype;

        std::vector<MatrixCreateKeyType> values;  //  all create values

        _Matrix::getCreateValues( values );
        
        for ( size_t i = 0; i < values.size(); ++i )
        {   
            if ( values[i].second != stype )
            {
                continue;
            }
            
            MatrixPtr<ValueType> matrixPtr( Matrix<ValueType>::getMatrix( values[i].first ) );

            if ( ctx )
            {
                matrixPtr->setContextPtr( ctx );
            }

            this->push_back( matrixPtr );
        }
    }
};
