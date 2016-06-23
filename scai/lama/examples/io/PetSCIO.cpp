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
#include "IOStream.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>

namespace scai
{

using namespace hmemo;

namespace lama
{

static int MAT_FILE_CLASSID = 1211216;  //<! internal id as specified by PetSC

static int VEC_FILE_CLASSID = 1211214;  //<! internal id as specified by PetSC

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( PetSCIO::logger, "FileIO.PetSCIO" )

/* --------------------------------------------------------------------------------- */

static std::string PETSC_SUFFIX   = ".psc";

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* PetSCIO::create()
{
    return new PetSCIO();
}   

std::string PetSCIO::createValue()
{
    return PETSC_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

void PetSCIO::writeAt( std::ostream& stream ) const
{
    stream << "PetSCIO ( ";
    stream << "suffix = " << PETSC_SUFFIX << ", ";
    writeMode( stream );
    stream << ", only binary, BIG Endian )";
}

/* --------------------------------------------------------------------------------- */

PetSCIO::PetSCIO() 
{
    if ( !mBinarySet )
    {
        mBinary = true;  // writes binary as default    
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_ASSERT( !mBinarySet || mBinary, "Formatted output not available for MatlabIO" )

    // int    VEC_FILE_CLASSID
    // int    number of rows
    // type   values

    int nrows = array.size();

    std::ios::openmode flags = std::ios::out | std::ios::binary;

    if ( mAppendMode )
    {
        flags |= std::ios::app;
    }
    else
    {
        flags |= std::ios::trunc;
    }

    IOStream outFile( fileName, flags, IOStream::BIG );

    SCAI_LOG_INFO( logger, "File " << fileName << " now open for binary write, append = " << mAppendMode )

    utilskernel::LArray<IndexType> headValues( 2 );

    headValues[0] = VEC_FILE_CLASSID;
    headValues[1] = nrows;

    outFile.writeBinary( headValues, mScalarTypeIndex );
    outFile.writeBinary( array, mScalarTypeData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    // int    VEC_FILE_CLASSID
    // int    number of rows
    // type   values

    std::ios::openmode flags = std::ios::in | std::ios::binary;

    SCAI_LOG_INFO( logger, "Read array<" << common::TypeTraits<ValueType>::id() 
                           << "> from file " << fileName << ", type = " << mScalarTypeData )

    IOStream inFile( fileName, flags, IOStream::BIG );

    utilskernel::LArray<int> headerVals;

    inFile.readBinary( headerVals, 2, common::scalar::INDEX_TYPE );

    int classid = headerVals[0];
    int n       = headerVals[1];

    SCAI_LOG_INFO( logger, "Read: id = " << classid << ", #n = " << n )

    SCAI_ASSERT_EQUAL( VEC_FILE_CLASSID, classid, "illegal VEC_FILE_CLASSID" )

    inFile.readBinary( array, n, mScalarTypeData );
 
    inFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    SCAI_ASSERT( !mBinarySet || mBinary, "Formatted output not available for MatlabIO" )

    // int    MAT_FILE_CLASSID
    // int    number of rows
    // int    number of columns
    // int    total number of nonzeros
    // int    *number nonzeros in each row
    // int    *column indices of all nonzeros (starting index is zero)

    int nrows = storage.getNumRows();
    int ncols = storage.getNumColumns();

    HArray<IndexType> csrIA;    // first offsets, later sizes
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    storage.buildCSRData( csrIA, csrJA, csrValues );

    int nnz = csrJA.size();

    // we need the CSR sizes, not the offsets

    utilskernel::HArrayUtils::unscan( csrIA );

    std::ios::openmode flags = std::ios::out | std::ios::binary;

    if ( mAppendMode )
    {
        flags |= std::ios::app;
    }
    else
    {
        flags |= std::ios::trunc;
    }

    IOStream outFile( fileName, flags, IOStream::BIG );

    SCAI_LOG_INFO( logger, "File " << fileName << " now open for binary write, append = " << mAppendMode )

    // Note: PetSC starts indexing with 0

    utilskernel::LArray<IndexType> headValues( 4 );

    headValues[0] = MAT_FILE_CLASSID;
    headValues[1] = nrows;
    headValues[2] = ncols;
    headValues[3] = nnz;

    // for binary output we make conversions to mScalarTypeData, mScalarTypeIndex

    outFile.writeBinary( headValues, mScalarTypeIndex );
    outFile.writeBinary( csrIA, mScalarTypeIndex );
    outFile.writeBinary( csrJA , mScalarTypeIndex ); 
    outFile.writeBinary( csrValues, mScalarTypeData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void PetSCIO::readStorageImpl(
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

    SCAI_LOG_INFO( logger, "Read storage<" << common::TypeTraits<ValueType>::id() << "> from file " << fileName )

    IOStream inFile( fileName, flags, IOStream::BIG );

    utilskernel::LArray<int> headerVals;

    inFile.readBinary( headerVals, 4, common::TypeTraits<int>::stype );

    int classid = headerVals[0];
    int nrows   = headerVals[1];
    int ncols   = headerVals[2];
    int nnz     = headerVals[3];

    SCAI_LOG_INFO( logger, "Read: id = " << MAT_FILE_CLASSID << ", #rows = " << nrows 
                           << ", #cols = " << ncols << ", #nnz = " << nnz )

    SCAI_ASSERT_EQUAL( MAT_FILE_CLASSID, classid, "illegal MAT_FILE_CLASSID" )

    HArray<IndexType> csrSizes;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    inFile.readBinary( csrSizes, nrows, mScalarTypeIndex );
    inFile.readBinary( csrJA, nnz, mScalarTypeIndex );
    inFile.readBinary( csrValues, nnz, mScalarTypeData );

    storage.setCSRData( nrows, ncols, nnz, csrSizes, csrJA, csrValues );
}

}  // lama

}  // scai
