/**
 * @file Storages.hpp
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
 * @brief Vector with all storages, one for each supported matrix storage format/type
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <scai/lama/storage/MatrixStorage.hpp>

#include <scai/common/shared_ptr.hpp>
#include <scai/common/TypeTraits.hpp>

#include <vector>

typedef scai::common::shared_ptr<scai::lama::_MatrixStorage> StoragePtr;

/** Class for a list of matrix storage pointers, one for each supported
 *  matrix storage format and each supported arithmetic type.
 */

class Storages : public std::vector<StoragePtr>
{

public:

    /** Constructor creates already the list with all storage pointers. */

    Storages( scai::hmemo::ContextPtr ctx = scai::hmemo::ContextPtr() )
    {
        using namespace scai::lama;
        std::vector<MatrixStorageCreateKeyType> values;  //  all create values
        _MatrixStorage::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            StoragePtr storagePtr( scai::lama::_MatrixStorage::create( values[i] ) );

            if ( ctx )
            {
                storagePtr->setContextPtr( ctx );
            }

            push_back( storagePtr );
        }
    }

    // Destructor will free all matrix storages due to use of shared pointers
};

/** Class for a list of typed matrix storage pointers, one for each supported
 *  matrix storage format
 */

template<typename ValueType>
class TypedStorages : public std::vector<scai::common::shared_ptr<scai::lama::MatrixStorage<ValueType> > >
{
public:

    /** Constructor creates already the list with all storage pointers. */

    TypedStorages( scai::hmemo::ContextPtr ctx = scai::hmemo::ContextPtr() )
    {
        using namespace scai::lama;
        using namespace scai::common;
        std::vector<MatrixStorageCreateKeyType> values;  //  all create values
        _MatrixStorage::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            scalar::ScalarType stype = values[i].second;

            if ( stype == TypeTraits<ValueType>::stype )
            {
                _MatrixStorage* storage = _MatrixStorage::create( values[i] );
                MatrixStorage<ValueType>* typedStorage = dynamic_cast<MatrixStorage<ValueType>*>( storage );
                SCAI_ASSERT( typedStorage, "dynamic cast failed" )
                shared_ptr<MatrixStorage<ValueType> > typedStoragePtr( typedStorage );

                if ( ctx )
                {
                    typedStoragePtr->setContextPtr( ctx );
                }

                this->push_back( typedStoragePtr );
            }
        }
    }

    // Destructor will free all matrix storages due to use of shared pointers
};
