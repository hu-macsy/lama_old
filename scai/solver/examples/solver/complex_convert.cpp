/**
 * @file solver/examples/solver/complex_convert.cpp
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
 * @brief Conversion program for complex matrices/vectors to real matrices/vectors
 * @author Thomas Brandes
 * @date 15.06.2016
 */

#include <scai/lama.hpp>

#include <scai/utilskernel/LArray.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>

using namespace std;

using namespace scai;

using common::TypeTraits;

using namespace hmemo;
using namespace lama;
using namespace utilskernel;

template<typename RealValueType, typename ComplexValueType>
void convertCmplx2Real( HArray<IndexType>& outIA,
                        HArray<IndexType>& outJA,
                        HArray<RealValueType>& outValues,
                        const HArray<IndexType>& inIA,
                        const HArray<IndexType>& inJA,
                        const HArray<ComplexValueType>& inValues )
{
    IndexType nnz = inIA.size();

    SCAI_ASSERT_EQUAL( nnz, inJA.size(), "mismatch" )
    SCAI_ASSERT_EQUAL( nnz, inValues.size(), "mismatch" )

    // ComplexValueType must be complex, RealValueType must not

    common::ScalarType inType = TypeTraits<ComplexValueType>::stype;
    common::ScalarType outType = TypeTraits<RealValueType>::stype;

    SCAI_ASSERT( common::isComplex( inType ), inType << " illegal, must be complex" );
    SCAI_ASSERT( !common::isComplex( outType ), outType << " illegal, must not be complex" );

    // size changes from 2 x n to 2 x n, nnz increases by factor 4

    ContextPtr ctx = Context::getHostPtr();

    WriteOnlyAccess<IndexType> wia( outIA, ctx, 4 * nnz );
    WriteOnlyAccess<IndexType> wja( outJA, ctx, 4 * nnz );
    WriteOnlyAccess<RealValueType> wvals( outValues, ctx, 4 * nnz );

    ReadAccess<IndexType> ria( inIA, ctx );
    ReadAccess<IndexType> rja( inJA, ctx );
    ReadAccess<ComplexValueType> rvals( inValues, ctx );

    for ( IndexType k = 0; k < nnz; ++k )
    {
        RealValueType real = common::Math::real( rvals[ k ] );
        RealValueType imag = common::Math::imag( rvals[ k ] );

        IndexType i = ria[k];
        IndexType j = rja[k];

        wia[4 * k]     = 2 * i;
        wja[4 * k]     = 2 * j;
        wvals[4 * k]   = real;

        wia[4 * k + 1]     = 2 * i;
        wja[4 * k + 1]     = 2 * j + 1;
        wvals[4 * k + 1]   = -imag;

        wia[4 * k + 2]     = 2 * i + 1;
        wja[4 * k + 2]     = 2 * j;
        wvals[4 * k + 2]   = imag;

        wia[4 * k + 3]     = 2 * i + 1;
        wja[4 * k + 3]     = 2 * j + 1;
        wvals[4 * k + 3]   = real;
    }
}

template<typename RealValueType, typename ComplexValueType>
void convertCmplx2Real( HArray<RealValueType>& outValues,
                        const HArray<ComplexValueType>& inValues )
{
    IndexType n = inValues.size();

    // ComplexValueType must be complex, RealValueType must not

    SCAI_ASSERT( common::isComplex( TypeTraits<ComplexValueType>::stype ), "must be complex" );
    SCAI_ASSERT( !common::isComplex( TypeTraits<RealValueType>::stype ), "must not be complex" );

    // size changes from 2 x n to 2 x n, nnz increases by factor 4

    ContextPtr ctx = Context::getHostPtr();

    WriteOnlyAccess<RealValueType> wvals( outValues, ctx, 2 * n );
    ReadAccess<ComplexValueType> rvals( inValues, ctx );

    for ( IndexType k = 0; k < n; ++k )
    {
        RealValueType real = common::Math::real( rvals[ k ] );
        RealValueType imag = common::Math::imag( rvals[ k ] );

        wvals[2 * k]     = real;
        wvals[2 * k + 1] = imag;
    }
}

int main( int argc, const char* argv[] )
{
    typedef ComplexDouble ComplexValueType;

    // take the abs type as real value type

    typedef TypeTraits<ComplexValueType>::RealType RealValueType;

    common::Settings::parseArgs( argc, argv );

    File::FileType fileType = File::MATRIX_MARKET;

    bool binary = false;
    common::ScalarType valueType = TypeTraits<RealValueType>::stype;
    common::ScalarType indexType = TypeTraits<IndexType>::stype;

    if ( argc != 3 )
    {
        cerr << "Usage: complex_convert complex_infile_name real_outfile_name" << endl;
        return -1;
    }

    const char* in_filename = argv[1];

    const char* out_filename = argv[2];

    HArray<ComplexValueType> values;

    int numColumns;

    StorageIO<ComplexValueType>::readDenseFromFile( values, numColumns, in_filename );

    cout << in_filename << ", read values = " << values << ", numColumns = " << numColumns << endl;

    if ( numColumns == 1 )
    {
        cout << "convert vector<" << TypeTraits<ComplexValueType>::id() << ">"
             << " to vector<" << TypeTraits<RealValueType>::id() << ">" << endl;

        HArray<RealValueType> newValues;
        convertCmplx2Real( newValues, values );

        cout << "new values = " << newValues << endl;

        StorageIO<RealValueType>::writeDenseToFile( newValues, 1, out_filename, fileType, valueType, binary );

        cout << "written vector to " << out_filename << ".mtx ( MatrixMarket )" << endl;
    }
    else
    {
        COOStorage<ComplexValueType> coo;

        coo.readFromFile( in_filename );

        int numRows = coo.getNumRows();
        numColumns  = coo.getNumColumns();

        cout << "COO storage = " << coo << endl;

        const HArray<IndexType>& ia = coo.getIA();
        const HArray<IndexType>& ja = coo.getJA();
        const HArray<ComplexValueType> values = coo.getValues();

        cout << "coo::ia = " << ia << endl;
        cout << "coo::ja = " << ja << endl;
        cout << "coo::values = " << values << endl;

        HArray<IndexType> newIA;
        HArray<IndexType> newJA;
        HArray<RealValueType> newValues;

        convertCmplx2Real( newIA, newJA, newValues, ia, ja, values );

        cout << "new::ia = " << newIA << endl;
        cout << "new::ja = " << newJA << endl;
        cout << "new::values = " << newValues << endl;

        COOStorage<double> newCOO( 2 * numRows, 2 * numColumns, newIA, newJA, newValues );

        newCOO.writeToFile( out_filename, fileType, valueType, indexType, indexType, binary );

        cout << "Output matrix file written: " << out_filename << ".mtx ( _Matrix Market)" << endl;
    }
}

