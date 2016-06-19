/**
 * @file PetSCIO.cpp
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
 * @brief Implementation of IO methods for PetSC format
 * @author Thomas Brandes
 * @date 10.06.2016
 */


#include "PetSCIO.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/lama/io/FileStream.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( PetSCIO::logger, "PetSCIO" )

/* --------------------------------------------------------------------------------- */

bool PetSCIO::isSupported( const bool binary ) const
{
    if ( binary )
    {
        return true; // binary is supported
    }
    else
    {
        return false; // formatted is unsupported
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::writeArrayFormatted(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    COMMON_THROWEXCEPTION( "writeArrayFormatted " << array << " to file " << fileName << " unsupported" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::writeArrayBinary(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName,
    const common::scalar::ScalarType valuesType) 
{
    COMMON_THROWEXCEPTION( "writeArrayBinary " << array << " to file "
                            << fileName << ", type = " << valuesType << " unsupported" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::readArrayTyped(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    COMMON_THROWEXCEPTION( "readArrayTyped " << array << " from file " << fileName << " unsupported" )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::writeStorageFormatted(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    COMMON_THROWEXCEPTION( "writeStorageFormatted< " << common::TypeTraits<ValueType>::id() << "> "
                           ", formatted not supported, storage = " << storage << ", filename = " << fileName )
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::writeStorageBinary(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName,
    const common::scalar::ScalarType iaType,
    const common::scalar::ScalarType jaType,
    const common::scalar::ScalarType valuesType ) 
{
    // int    MAT_FILE_CLASSID
    // int    number of rows
    // int    number of columns
    // int    total number of nonzeros
    // int    *number nonzeros in each row
    // int    *column indices of all nonzeros (starting index is zero)

    int nrows = storage.getNumRows();
    int ncols = storage.getNumColumns();

    HArray<IndexType> csrIA;
    HArray<IndexType> csrSizes;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    storage.buildCSRData( csrIA, csrJA, csrValues );

    int nnz = csrJA.size();
    int MAT_FILE_CLASSID = 0;

    // we need the CSR sizes, not the offsets

    {
        ContextPtr loc = hmemo::Context::getHostPtr();
        static utilskernel::LAMAKernel<sparsekernel::CSRKernelTrait::offsets2sizes > offsets2sizes;
        offsets2sizes.getSupportedContext( loc );
        ReadAccess<IndexType> rOffsets( csrIA, loc );
        WriteOnlyAccess<IndexType> wSizes( csrSizes, loc, nrows );
        SCAI_CONTEXT_ACCESS( loc )
        offsets2sizes[ loc ]( wSizes.get(), rOffsets.get(), nrows );
    }

    std::ios::openmode flags = std::ios::out | std::ios::trunc | std::ios::binary;

    FileStream outFile( fileName, flags );

    std::cout << "File " << fileName << " now open for binary write" << std::endl;

    // Note: PetSC starts indexing with 0

    std::fstream& f = outFile;   // avoid conflicts for method write on base class fstream

    f.write( reinterpret_cast<const char*>( &MAT_FILE_CLASSID ), sizeof( IndexType ) );
    f.write( reinterpret_cast<const char*>( &nrows ), sizeof( IndexType ) );
    f.write( reinterpret_cast<const char*>( &ncols ), sizeof( IndexType ) );
    f.write( reinterpret_cast<const char*>( &nnz ), sizeof( IndexType ) );

    outFile.write<IndexType>( csrSizes, 0, iaType, '\n' );
    outFile.write<IndexType>( csrJA , 0, jaType, '\n' ); 
    outFile.write<ValueType>( csrValues, 0, valuesType, '\n' );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::readStorageTyped(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName ) 
{
    // int    MAT_FILE_CLASSID
    // int    number of rows
    // int    number of columns
    // int    total number of nonzeros
    // int    *number nonzeros in each row
    // int    *column indices of all nonzeros (starting index is zero)

    std::ios::openmode flags = std::ios::in | std::ios::binary;

    std::cout << "Read from file " << fileName << std::endl;

    FileStream inFile( fileName, flags );

    std::fstream& f = inFile;   // avoid conflicts for method write on base class fstream

    int MAT_FILE_CLASSID;
    int nrows;
    int ncols;
    int nnz;

    f.read( reinterpret_cast<char*>( &MAT_FILE_CLASSID ), sizeof( IndexType ) );
    f.read( reinterpret_cast<char*>( &nrows ), sizeof( IndexType ) );
    f.read( reinterpret_cast<char*>( &ncols ), sizeof( IndexType ) );
    f.read( reinterpret_cast<char*>( &nnz ), sizeof( IndexType ) );

    std::cout << "Read: id = " << MAT_FILE_CLASSID << ", #rows = " << nrows 
              << ", #cols = " << ncols << ", #nnz = " << nnz << std::endl;

    HArray<IndexType> csrSizes;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    ValueType valAdjust = 0;

    inFile.read( csrSizes, nrows, 0, common::TypeTraits<IndexType>::stype, '\n' );
    inFile.read( csrJA, nnz, 0, common::TypeTraits<IndexType>::stype, '\n' );
    inFile.read( csrValues, nnz, valAdjust, common::TypeTraits<ValueType>::stype, '\n' );

    storage.setCSRData( nrows, ncols, nnz, csrSizes, csrJA, csrValues );
}

}  // lama

}  // scai
