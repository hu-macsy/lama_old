/**
 * @file MatlabIO.cpp
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
 * @brief Implementation of methods for FileIO class MatlabIO
 * @author Thomas Brandes
 * @date 10.06.2016
 */


#include "MatlabIO.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/io/IOStream.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/exception/IOException.hpp>

#include <sstream>

#define MATLAB_SUFFIX ".txt" 

using namespace std;

namespace scai
{

using namespace hmemo;
using namespace utilskernel;

namespace lama
{

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* MatlabIO::create()
{
    return new MatlabIO();
}

std::string MatlabIO::createValue()
{
    return MATLAB_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

bool MatlabIO::isSupportedMode( const FileMode mode ) const
{
    // binary is not supported

    if ( mode == BINARY )
    {
        return false;
    }

    return true;
}

/* --------------------------------------------------------------------------------- */

void MatlabIO::writeAt( std::ostream& stream ) const
{
    stream << "MatlabIO ( suffix = " << MATLAB_SUFFIX << ", ";
    writeMode( stream );
    stream << ", only formatted )";
}

/** Method to count number of lines of a text file and the maximal number of entries in one line 
 *
 *  @param[out]  number of lines the file has
 *  @param[out]  nEntries is maximal number of entries
 *  @param[in]   fileName is the name of the file
 *
 *  Note: it might be possible that one line contains less than 'nEntries' entries
 */
void MatlabIO::checkTextFile( IndexType& nLines, IndexType& nEntries, const char* fileName )
{
    nLines   = 0;
    nEntries = 0;

    std::ifstream infile( fileName, std::ios::in );

    if ( infile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    std::string line;
    std::vector<std::string> tokens;

    while ( std::getline( infile, line ) )
    {
        ++nLines;

        common::Settings::tokenize( tokens, line );

        IndexType nTokens = tokens.size();

        if ( nTokens > nEntries )
        {
            nEntries = nTokens;
            // LOG_DEBUG: cout << "max tokens = " << nEntries << " at line " << nLines << endl;
        }
    }

    SCAI_LOG_INFO( logger, "checkTextFile " << fileName << ": #lines = " << nLines << ", #entries = " << nEntries )
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MatlabIO::logger, "FileIO.MatlabIO" )

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_ASSERT( mFileMode != BINARY, "Binary mode not supported for " << *this )

    IOStream outFile( fileName, std::ios::out );

    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    outFile.writeFormatted( array, precData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readArrayImpl(
    hmemo::HArray<ValueType>& array,
    const std::string& fileName ) 
{
    int n;   // number of lines, size of array
    int k;   // numbe of entries in one line

    checkTextFile( n, k, fileName.c_str() );

    SCAI_LOG_INFO( logger, "File : " << fileName << ", #lines = " << n << ", #entries = " << k )

    SCAI_ASSERT_LE( k, 2, "#entries/row in file " << fileName << " must not excced 2" )

    // use local arrays instead of heteregeneous arrays as we want ops on them

    IOStream inFile( fileName, std::ios::in );

    inFile.readFormatted( array, n );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName ) 
{
    SCAI_ASSERT( mFileMode != BINARY, "Binary mode not supported for " << *this )

    COOStorage<ValueType> coo( storage );

    HArray<IndexType> cooIA = coo.getIA();
    HArray<IndexType> cooJA = coo.getJA();
    HArray<ValueType> cooValues = coo.getValues();

    IOStream outFile( fileName, std::ios::out );

    int precIndex = 0;
    int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

    outFile.writeFormatted( cooIA, precIndex, cooJA, precIndex, cooValues, precData );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    // binary mode does not matter as we have always formatted output

    // read coo entries lines by line, similiar to MatrixMarket
    // i , j, val     

    int nnz;
    int k;

    checkTextFile( nnz, k, fileName.c_str() );

    SCAI_LOG_INFO( logger, "File : " << fileName << ", #lines = " << nnz << ", #entries = " << k )

    if ( nnz == 0 )
    {
        storage.clear();
        return;
    }

    SCAI_ASSERT_GE( k, 3, "#entries/row in file " << fileName << " must be at least 3" )

    // use local arrays instead of heteregeneous arrays as we want ops on them

    LArray<double> dIA;
    LArray<double> dJA;
    LArray<ValueType> val;

    IOStream inFile( fileName, std::ios::in );

    inFile.readFormatted( dIA, dJA, val, nnz );

    LArray<IndexType> ia( dIA );
    LArray<IndexType> ja( dJA );

    SCAI_LOG_DEBUG( logger, "read ia  : " << ia  )
    SCAI_LOG_DEBUG( logger, "read ja  : " << ja  )
    SCAI_LOG_DEBUG( logger, "read val : " << val )

    int minRowIndex = ia.min();

    if ( minRowIndex == 0 )
    {
        // okay, seems that indexing start with 0
    }
    else if ( minRowIndex == 1 )
    {
        // offset base = 1, convert it to 0

        ia -= minRowIndex;
        ja -= minRowIndex;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Index base = " << minRowIndex << " is illegal" )
    }

    // we shape the matrix by maximal appearing indexes

    int nrows = ia.max() + 1;
    int ncols = ja.max() + 1;

    COOStorage<ValueType> coo( nrows, ncols, ia, ja, val );

    storage = coo;
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
