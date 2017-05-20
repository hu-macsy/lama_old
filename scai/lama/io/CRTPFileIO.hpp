/**
 * @file CRTPFileIO.hpp
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
 * @brief CRTP wrapper class for all FileIO classes
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#pragma once

#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/IOWrapper.hpp>
#include <scai/lama/storage/MatrixStorage.hpp>

#include <cstdio>

namespace scai
{

namespace lama
{

/* --------------------------------------------------------------------------------- */

/** This help class uses the CRTP pattern to provide for all derived FileIO classes
 *  methods that are implemented in the same way.
 *
 *  - call typed version of writeStorage and readStorage
 *  - get default routines for file suffixes
 *
 *  @tparam Derived must be a class that is derived from this CRTP class.
 *
 *  By this way, implementation of new FileIO classes is very simple.
 */
template<class Derived>
class CRTPFileIO : public FileIO
{

public:

    /** Implementation of pure virtual method FileIO::writeStorage for all derived classes */

    void writeStorage( const _MatrixStorage& storage, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::readStorage  */

    void readStorage( _MatrixStorage& storage, const std::string& fileName, const IndexType offsetRow, const IndexType nRows );

    /** Implementation of pure virtual method FileIO::writeArray  */

    virtual void writeArray( const hmemo::_HArray& array, const std::string& fileName );

    /** Implementation of pure virtual method FileIO::writeSparse  */

    virtual void writeSparse( 
        const IndexType size, 
        const hmemo::HArray<IndexType>& indexes, 
        const hmemo::_HArray& array, 
        const std::string& fileName );

    /** Implementation of pure virtual method FileIO::readArray using same defaults */

    virtual void readArray(
        hmemo::_HArray& array,
        const std::string& fileName,
        const IndexType offset = 0,
        const IndexType n = nIndex );

    /** Implementation of pure virtual method FileIO::readSparse 
     *
     *  This CRTP class calls Derived::readSparseImpl with a typed value array.
     */

    virtual void readSparse(
        IndexType& size,
        hmemo::HArray<IndexType>& indexes,
        hmemo::_HArray& values,
        const std::string& fileName );

    /** Default implementation for query matrix file suffix, is createValue of derived class */

    virtual std::string getMatrixFileSuffix() const;

    /** Default implementation for query vector file suffix, is createValue of derived class */

    virtual std::string getVectorFileSuffix() const;
};

/* --------------------------------------------------------------------------------- */
/*    Implementation of methods for CRTPFileIO                                       */
/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::writeStorage( const _MatrixStorage& storage, const std::string& fileName )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( FileIO::hasSuffix( fileName, this->getMatrixFileSuffix() ),
                 fileName << " illegal file name for storage, must have suffix " << getMatrixFileSuffix() )

    IOWrapper<Derived, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorageImpl( ( Derived& ) *this, storage, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::readStorage(
    _MatrixStorage& storage,
    const std::string& fileName,
    const IndexType offsetRow,
    const IndexType nRows )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( FileIO::hasSuffix( fileName, getMatrixFileSuffix() ),
                 fileName << " illegal, must have suffix " << getMatrixFileSuffix() )

    // just call the corresponding typed routine

    IOWrapper<Derived, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorageImpl( ( Derived& ) *this, storage, fileName, offsetRow, nRows );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::writeArray( const hmemo::_HArray& array, const std::string& fileName )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( FileIO::hasSuffix( fileName, this->getVectorFileSuffix() ),
                 fileName << " illegal file name for array, must have suffix " << this->getVectorFileSuffix() )

    // now call the corresponding typed routine, use meta-programming to get the correct type

    IOWrapper<Derived, SCAI_ARRAY_TYPES_HOST_LIST>::writeArrayImpl( ( Derived& ) *this, array, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::writeSparse( const IndexType n, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values, const std::string& fileName )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( FileIO::hasSuffix( fileName, this->getVectorFileSuffix() ),
                 fileName << " illegal file name for array, must have suffix " << this->getVectorFileSuffix() )

    // now call the corresponding typed routine, use meta-programming to get the correct type

    IOWrapper<Derived, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparseImpl( ( Derived& ) *this, n, indexes, values, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::readArray( hmemo::_HArray& array, const std::string& fileName, const IndexType offset, const IndexType n )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( FileIO::hasSuffix( fileName, this->getVectorFileSuffix() ),
                 fileName << " illegal file name for array, must have suffix " << this->getVectorFileSuffix() )

    // just call the corresponding typed routine

    IOWrapper<Derived, SCAI_ARRAY_TYPES_HOST_LIST>::readArrayImpl( ( Derived& ) *this, array, fileName, offset, n );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
void CRTPFileIO<Derived>::readSparse( IndexType& size, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values, const std::string& fileName )
{
    SCAI_ASSERT( fileName.size() > 0 , "Error: fileName should not be empty" )

    SCAI_ASSERT( FileIO::hasSuffix( fileName, this->getVectorFileSuffix() ),
                 fileName << " illegal file name for array, must have suffix " << this->getVectorFileSuffix() )

    // just call the corresponding typed routine

    IOWrapper<Derived, SCAI_ARRAY_TYPES_HOST_LIST>::readSparseImpl( ( Derived& ) *this, size, indexes, values, fileName );
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
std::string CRTPFileIO<Derived>::getMatrixFileSuffix() const
{
    return Derived::createValue();
}

/* --------------------------------------------------------------------------------- */

template<class Derived>
std::string CRTPFileIO<Derived>::getVectorFileSuffix() const
{
    return Derived::createValue();
}

/* --------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai
