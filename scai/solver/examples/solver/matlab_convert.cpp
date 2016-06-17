/**
 * @file solver/examples/solver/matlab_convert.cpp
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
 * @brief Conversion program for Matlab matrices / vectors
 * @author Thomas Brandes
 * @date 15.06.2016
 */

#include <scai/lama.hpp>

#include <scai/utilskernel/LArray.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>

using namespace std;

using namespace scai;
using namespace hmemo;
using namespace lama;
using namespace utilskernel;

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

    std::ifstream infile( fileName, ios::in );

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

template<typename ValueType>
bool readVal( ValueType& val, std::istringstream& input )
{
    return input >> val;
}

template<>
bool readVal( IndexType& val, std::istringstream& input )
{
    // int can als be represented as "1.4280000e+03" 

    double x;

    bool ok = input >> x;

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
    std::ifstream infile( fileName, ios::in );

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

template<typename ValueType1>
void readTextFile( HArray<ValueType1>& val1,
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

    WriteOnlyAccess<ValueType1> wVal1( val1, ctx, nlines );

    std::string line;

    for ( IndexType k = 0; k < nlines; ++k )
    {
        std::getline( infile, line );

        if ( infile.fail() )
        {
            COMMON_THROWEXCEPTION( "Line mismatch at line " << k << ", expected " << nlines << " lines" )
        }

        std::istringstream iss( line );

        readVal( wVal1[k], iss );
    }
}

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    File::FileType fileType = File::MATRIX_MARKET;

    common::scalar::ScalarType valueType = common::TypeTraits<ComplexDouble>::stype;
    common::scalar::ScalarType indexType = common::TypeTraits<IndexType>::stype;

    bool binary = false;

    if ( argc != 3 ) 
    {
        cerr << "Usage: matlab_convert infile_name outfile_name" << endl;
        return -1;
    }

    const char* in_filename = argv[1];
    const char* out_filename = argv[2];

    LArray<IndexType> ia;
    LArray<IndexType> ja;
    LArray<ComplexDouble> val;

    int nnz;
    int k;

    checkTextFile( nnz, k, in_filename );

    cout << "File : " << in_filename << ", #lines = " << nnz << ", #entries = " << k << endl;

    if ( k == 2 )
    {
        HArray<ComplexDouble> val;

        readTextFile( val, nnz, in_filename );
 
        cout << "File " << in_filename << " read, val = " << val << endl;

        StorageIO<ComplexDouble>::writeDenseToFile( val, 1, out_filename, fileType, valueType, binary );

        cout << "Output file written: " << out_filename << endl;

        return 0;
    }

    readTextFile( ia, ja, val, nnz, in_filename );

    cout << "Read ia : " << ia << endl;
    cout << "Read ja : " << ja << endl;
    cout << "Read vals : " << val << endl;

    int N = ia.min();

    if ( N == 0 )
    {
        // okay
    }
    else if ( N == 1 )
    {
        // offset base = 1, convert it to 0

        ia -= 1;
        ja -= 1;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Index base = " << N << " is illegal" )
    }

    N = ia.max() + 1;

    int M = ja.max() + 1;

    std::cout << "Matrix size = " << N << " x " << M << std::endl;

    COOStorage<ComplexDouble> coo( N, M, ia, ja, val );

    cout << "COO storage = " << coo << endl;

    CSRStorage<ComplexDouble> csr( coo );  // converts also COO to CSR

    cout << "CSR storage = " << csr << endl;

    csr.writeToFile( out_filename, fileType, valueType, indexType, indexType, binary );

    cout << "Output matrix file written: " << out_filename << endl;
}
