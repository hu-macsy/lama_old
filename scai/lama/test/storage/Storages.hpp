/**
 * @file Storages.hpp
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
 * @brief Vector with all storages, one for each supported matrix storage format/type
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
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

                this->push_back( typedStoragePtr );
            }
        }

        if ( ctx )
        {
            for ( size_t i = 0; i < values.size(); ++i )
            {
                (*this)[i]->setContextPtr( ctx );
            }
        }
    }

    // Destructor will free all matrix storages due to use of shared pointers
};
