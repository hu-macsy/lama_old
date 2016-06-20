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
 * @brief Implementation of methods
 * @author Thomas Brandes
 * @date 10.06.2016
 */


#include "MatlabIO.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/lama/io/FileStream.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>

#include <sstream>

using namespace std;

template<typename ValueType>
bool readVal( ValueType& val, std::istringstream& input )
{
    input >> val;

    bool ok = !input.fail();

    return ok;
}

template<>
bool readVal( IndexType& val, std::istringstream& input )
{
    // int can als be represented as "1.4280000e+03" 

    double x;

    input >> x;

    bool ok = !input.fail();

    if ( ok )
    {
        val = static_cast<IndexType>( x );
        // do not accept e.g. 1.53
        ok  = static_cast<double>( val ) == x;
    }
    else
    {
        val = 0;
    }

    return ok;
}

namespace scai
{

using namespace hmemo;
using namespace utilskernel;

namespace lama
{

static std::string MATLAB_SUFFIX   = ".txt";

std::string MatlabIO::getVectorFileSuffix() const
{
    return MATLAB_SUFFIX;
}

std::string MatlabIO::getMatrixFileSuffix() const
{
    return MATLAB_SUFFIX;
}

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

void MatlabIO::writeAt( std::ostream& stream ) const
{
    stream << "MatlabIO (only formatted)";
}

/** Method to count number of lines of a text file and the maximal number of entries in one line 
 *
 *  @param[out]  number of lines the file has
 *  @param[out]  nEntries is maximal number of entries
 *  @param[in]   fileName is the name of the file
 *
 *  Note: it might be possible that one line contains less than 'nEntries' entries
 */
void checkTextFile( IndexType& nLines, IndexType& nEntries, const char* fileName )
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
}


/* --------------------------------------------------------------------------------- */

/** Read in a set of arrays from a file, one entry for each array in a line
 *
 *  @param[out] val1 val2 val3 are the arrays whose values are read, size will be nlines
 *  @param[in]  nlines is number of lines
 *  @param[in]  fileName is the name of the input file
 *
 *  Be careful: if ValueType is complex, two entries must be in the input file
 */

template<typename ValueType1, typename ValueType2, typename ValueType3>
void readTextFile( HArray<ValueType1>& val1,
                   HArray<ValueType2>& val2,
                   HArray<ValueType3>& val3,
                   const IndexType nlines,
                   const char* fileName )
{
    std::ifstream infile( fileName, std::ios::in );

    if ( infile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    ContextPtr ctx = Context::getHostPtr();

    // use two times the size if Value type is not complex

    std::cout << "Read " << fileName << ", nlines = " << nlines << std::endl;

    WriteOnlyAccess<ValueType1> wVal1( val1, ctx, nlines );
    WriteOnlyAccess<ValueType2> wVal2( val2, ctx, nlines );
    WriteOnlyAccess<ValueType3> wVal3( val3, ctx, nlines );

    bool error = false;

    std::string line;

    for ( IndexType k = 0; k < nlines; ++k )
    {
        std::getline( infile, line );

        if ( infile.fail() )
        {
            COMMON_THROWEXCEPTION( "Line mismatch at line " << k << ", expected " << nlines << " lines" )
        }

        std::istringstream iss( line );

        if ( readVal( wVal1[k], iss ) &&
             readVal( wVal2[k], iss ) &&
             readVal( wVal3[k], iss ) )
        {
            // okay
        }
        else
        {
            std::cout << "ERROR in " << fileName << ": Missed or wrong values in line " << k << std::endl;
            error = true;
        }
    }

    if ( error )
    {
        COMMON_THROWEXCEPTION( "File " << fileName << " contains errors" )
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void readTextFile( HArray<ValueType>& val,
                   const IndexType nlines,
                   const char* fileName )
{   
    std::ifstream infile( fileName, ios::in );
    
    if ( infile.fail() )
    {   
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }
    
    ContextPtr ctx = Context::getHostPtr();
    
    // use two times the size if Value type is not complex
    
    std::cout << "Read " << fileName << ", nlines = " << nlines << std::endl;
    
    WriteOnlyAccess<ValueType> wVal( val, ctx, nlines );
    
    std::string line;
    
    for ( IndexType k = 0; k < nlines; ++k )
    {   
        std::getline( infile, line );
        
        if ( infile.fail() )
        {   
            COMMON_THROWEXCEPTION( "Line mismatch at line " << k << ", expected " << nlines << " lines" )
        }
        
        std::istringstream iss( line );
        
        readVal( wVal[k], iss );
    }
}

/* --------------------------------------------------------------------------------- */


/** write a set of arrays from a file, one entry for each array in a line
 *
 *  Be careful: if ValueType is complex, two entries must be in the output file
 */

template<typename ValueType1, typename ValueType2, typename ValueType3>
void writeTextFile( const HArray<ValueType1>& val1,
                    const HArray<ValueType2>& val2,
                    const HArray<ValueType3>& val3,
                    const char* fileName )
{
    int n = val1.size();

    SCAI_ASSERT_EQUAL( n, val2.size(), "size mismatch" );
    SCAI_ASSERT_EQUAL( n, val3.size(), "size mismatch" );

    std::fstream outfile( fileName, std::ios::out );

    if ( outfile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    ContextPtr ctx = Context::getHostPtr();

    // use two times the size if Value type is not complex

    ReadAccess<ValueType1> rVal1( val1, ctx );
    ReadAccess<ValueType2> rVal2( val2, ctx );
    ReadAccess<ValueType3> rVal3( val3, ctx );

    bool error = false;

    for ( IndexType k = 0; k < n; ++k )
    {
        outfile << rVal1[k] << ' ' << rVal2[k] << ' ' << rVal3[k] << std::endl;

        if ( outfile.fail() )
        {
            error = true;
        }
    }

    if ( error )
    {
        COMMON_THROWEXCEPTION( "File " << fileName << " could not be written completely" )
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void writeTextFile( const HArray<ValueType>& val,
                    const char* fileName )
{
    int n = val.size();

    std::fstream outfile( fileName, std::ios::out );

    if ( outfile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open file '" << fileName << "'." )
    }

    ContextPtr ctx = Context::getHostPtr();

    // use two times the size if Value type is not complex

    ReadAccess<ValueType> rVal( val, ctx );

    bool error = false;

    for ( IndexType k = 0; k < n; ++k )
    {
        outfile << rVal[k] << std::endl;

        if ( outfile.fail() )
        {
            error = true;
        }
    }

    if ( error )
    {
        COMMON_THROWEXCEPTION( "File " << fileName << " could not be written completely" )
    }
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MatlabIO::logger, "MatlabIO" )

/* --------------------------------------------------------------------------------- */

bool MatlabIO::isSupported( const bool binary ) const
{
    if ( binary )
    {
        return false; // binary is not supported
    }
    else
    {
        return true;  // formatted supported
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    SCAI_ASSERT( !mBinary, "Binary mode not supported for MatlabIO" )

    writeTextFile( array, fileName.c_str() );
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

    readTextFile( array, n, fileName.c_str() );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void MatlabIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName ) 
{
    SCAI_ASSERT( !mBinary, "Binary mode not supported for MatlabIO" )

    COOStorage<ValueType> coo( storage );

    HArray<IndexType> cooIA = coo.getIA();
    HArray<IndexType> cooJA = coo.getJA();
    HArray<ValueType> cooValues = coo.getValues();

    writeTextFile( cooIA, cooJA, cooValues, fileName.c_str() );
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

    SCAI_ASSERT_GE( k, 3, "#entries/row in file " << fileName << " must be at least 3" )

    // use local arrays instead of heteregeneous arrays as we want ops on them

    LArray<IndexType> ia;
    LArray<IndexType> ja;
    LArray<ValueType> val;

    readTextFile( ia, ja, val, nnz, fileName.c_str() );

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
