/**
 * @file SAMGIO.cpp
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
 * @brief Implementation of IO methods for SAMG format
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#include "SAMGIO.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/lama/io/FileStream.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>

namespace scai
{

using namespace hmemo;

namespace lama
{

static int SAMG_VERSION_ID = 22;
static int SAMG_IVERSION   = 4;

/** SAMG file suffixes */

static std::string SAMG_MAT_HEADER_SUFFIX = ".srm";
static std::string SAMG_MAT_DATA_SUFFIX   = ".samg";
static std::string SAMG_VEC_HEADER_SUFFIX = ".srv";
static std::string SAMG_VEC_DATA_SUFFIX   = ".svec";

std::string SAMGIO::getVectorFileSuffix() const
{
    return SAMG_VEC_HEADER_SUFFIX;
}
    
std::string SAMGIO::getMatrixFileSuffix() const
{
    return SAMG_MAT_HEADER_SUFFIX;
}

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* SAMGIO::create()
{
    return new SAMGIO();
}

std::string SAMGIO::createValue()
{
    return SAMG_MAT_HEADER_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeAt( std::ostream& stream ) const
{
    stream << "SAMGIO (binary/formatted)";
}

/* --------------------------------------------------------------------------------- */

/** Help routine to get data file name by header file name 
 *
 *  @param[in] headerFileName is the file name of header file
 *  @return    name of the data file
 *
 *  Note: returns same name if no distinction between header and data 
 */
static std::string getDataFileName( const std::string& headerFileName )
{
    std::string result = headerFileName;

    if ( _StorageIO::hasSuffix( headerFileName, SAMG_MAT_HEADER_SUFFIX) )
    {
        size_t len = SAMG_MAT_HEADER_SUFFIX.length();
        result.replace( result.length() - len, len, SAMG_MAT_DATA_SUFFIX );
    }
    else if ( _StorageIO::hasSuffix( headerFileName, SAMG_VEC_HEADER_SUFFIX ) )
    {
        size_t len = SAMG_VEC_HEADER_SUFFIX.length();
        result.replace( result.length() - len, len, SAMG_VEC_DATA_SUFFIX );
    }

    return result;   // same name if no distinction between header and data
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( SAMGIO::logger, "SAMGIO" )

/* --------------------------------------------------------------------------------- */

SAMGIO::SAMGIO() 
{
    if ( common::Settings::getEnvironment( mAppendMode, "SCAI_IO_APPEND" ) )
    {
        if ( mAppendMode )
        {
            SCAI_LOG_WARN( logger, "SAMG format does not support append mode" ) 
        }
    }
}

/* --------------------------------------------------------------------------------- */

bool SAMGIO::isSupported( const bool binary ) const
{
    if ( binary )
    {
        return true; // binary is supported
    }
    else
    {
        return binary; // formatted is unsupported
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    // start by writing the SAMG Vector header file

    char fileType = mBinary ? 'b' : 'f';

    // ToDo: conversion not supported here

    int typeSize = sizeof( ValueType );

    if ( typeSize == 0 )
    {
        switch ( mDataType )
        {
            case common::scalar::INTERNAL:
                typeSize = sizeof( ValueType );
                break;

            default:
                SCAI_LOG_ERROR( logger, "Encountered invalid scalar type " << mDataType )
                break;
        }
    }

    FileStream outFile( fileName, std::ios::out );
    outFile << fileType << std::endl;
    outFile << array.size() << std::endl;
    outFile << typeSize;
    outFile.close();

    // write data into SAMG vector data file

    std::ios::openmode flags = std::ios::out | std::ios::trunc;

    if ( mBinary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );

    outFile.open( dataFileName, flags );
    outFile.write<ValueType>( array, 0, mDataType, '\n' );
    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    char fileType = ' ';
    int dataTypeSize = 0;
    IndexType numRows;
    // start with reading the *.frv header file
    FileStream inFile( fileName, std::ios::in );
    inFile >> fileType;
    inFile >> numRows;
    inFile >> dataTypeSize;
    inFile.close();
    // *.frv files can only store vectors

    if ( fileType != 'b' && fileType != 'f' )
    {
        COMMON_THROWEXCEPTION( "Invalid SAMG vector header file: " << fileName
                               << ", fileType = " << fileType << " illegal" )
    }

    // now read *.vec file in correct mode
    std::ios::openmode flags = std::ios::in;

    if ( fileType == 'b' )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );

    inFile.open( dataFileName, flags );
    common::scalar::ScalarType dataType;

    switch ( dataTypeSize )
    {
        // special cases for handling IndexType, int and long, as these are not properly supported yet
        case 4:
            dataType = common::scalar::FLOAT;
            break;

        case 8:
            dataType = common::scalar::DOUBLE;
            break;

        default:
            dataType = common::scalar::UNKNOWN;
            SCAI_LOG_ERROR( logger, "Encountered invalid type size " << dataTypeSize )
    }

    inFile.read( array, numRows, ValueType( 0 ), dataType, '\n' );
    inFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    storage.buildCSRData( csrIA, csrJA, csrValues );

    char fileType = mBinary ? 'b' : 'f';

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    FileStream outFile( fileName, std::ios::out | std::ios::trunc );

    IndexType size = 1;
    IndexType rank = 0;

    outFile << fileType << " \t" << SAMG_IVERSION << "\n";
    outFile << "\t\t" << numValues << "\t" << numRows << "\t" << SAMG_VERSION_ID << "\t" << size << "\t" << rank;
    outFile.close();
    SCAI_LOG_INFO( logger, "writeCSRToBinaryFile ( " << fileName << ")" << ", #rows = " << csrIA.size() - 1
                   << ", #values = " << csrJA.size() )
    std::ios::openmode flags = std::ios::out | std::ios::trunc;

    if ( mBinary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );

    outFile.open( dataFileName, flags );
    outFile.write<IndexType>( csrIA, 1, mIAType, '\n' );
    outFile.write<IndexType>( csrJA, 1, mJAType, '\n' ); 
    outFile.write<ValueType>( csrValues, 0, mDataType, '\n' );
    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName ) 
{
    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    // start with reading the header
    FileStream inFile( fileName, std::ios::in );
    int iversion; 
    char fileType = '!'; 
    IndexType numValues, numRows, id, size, rank;
    numValues = 0;
    inFile >> fileType >> iversion;

    if ( fileType != 'f' && fileType != 'b' )
    {   
        COMMON_THROWEXCEPTION( "Invalid file type" )
    }

    SCAI_ASSERT_EQUAL( SAMG_IVERSION, iversion, "SAMG version mismatch" )

    inFile >> numValues;
    inFile >> numRows;
    // AMG is always square?
    IndexType numColumns = numRows;
    inFile >> id;
    inFile >> size;
    inFile >> rank; 
    inFile.close(); // explicitly, otherwise done by destructor
    // now read *.amg file in correct mode
    std::ios::openmode flags = std::ios::in;

    if ( fileType == 'b' )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );
    inFile.open( dataFileName, flags );
    //TODO: allow different type to be read?
    inFile.read( csrIA, numRows + 1, -1, common::TypeTraits<IndexType>::stype, '\n' );
    inFile.read( csrJA, numValues, -1, common::TypeTraits<IndexType>::stype, '\n' );
    inFile.read( csrValues, numValues, ValueType( 0 ), common::TypeTraits<ValueType>::stype, '\n' );
    inFile.close();

    storage.setCSRData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
