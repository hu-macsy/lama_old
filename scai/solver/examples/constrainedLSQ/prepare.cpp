/**
 * @file prepare.cpp
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
 * @brief Example of least square problem with boundary conditions
 * @author Thomas Brandes, Andreas Borgen Longva
 * @date 21.07.2017
 */

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>

using namespace scai;
using namespace lama;

void setupSmoothMatrix( CSRSparseMatrix<double>& L, const IndexType ny, const IndexType nz, const double strength )
{
    // not neccesary, but to keep in accordance with the original

    common::Stencil2D<double> stencil( 5 );

    common::Grid2D grid( ny, nz );

    grid.setBorderType( 0, common::Grid::BORDER_REFLECTING );
    grid.setBorderType( 1, common::Grid::BORDER_REFLECTING );

    StencilMatrix<double> stencilMatrix( grid, stencil );

    L = stencilMatrix;

    L.scale( - strength / 4 );
}

void joinMatrix( CSRSparseMatrix<double>& result, const CSRSparseMatrix<double>& a, const CSRSparseMatrix<double>& b )
{
    SCAI_ASSERT_EQ_ERROR( a.getNumColumns(), b.getNumColumns(), "joined matrices must have same number of columns" );

    typedef scai::common::shared_ptr<scai::lama::_MatrixStorage> StoragePtr;

    StoragePtr shared_ptrA( a.getLocalStorage().copy() );
    StoragePtr shared_ptrB( b.getLocalStorage().copy() );

    std::vector<StoragePtr> bothMatrices;

    bothMatrices.push_back( shared_ptrA );
    bothMatrices.push_back( shared_ptrB );

    scai::lama::CSRStorage<double> joinedStorage;

    joinedStorage.rowCat( bothMatrices );

    result.assign( joinedStorage );
}

void zeroExtend( DenseVector<double>& T_ext,
                 const DenseVector<double>& T, const IndexType nZeros )
{
    T_ext.allocate( T.size() + nZeros );

    hmemo::WriteAccess<double> wT( T_ext.getLocalValues() );
    hmemo::ReadAccess<double> rT( T.getLocalValues() );

    for ( IndexType i = 0; i < rT.size(); ++i )
    {
        wT[i] = rT[i];
    }
}

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.driver" )
    
    common::Settings::parseArgs( argc, argv );

    if ( argc <= 8 )
    {
        std::cout << argv[0] << " <input_D> <input_T> <input_So> <input_hrz> " << std::endl;
        std::cout << "            <output_A> <output_b> <output_lb> <output_ub>" << std::endl;
        std::cout << "            [ --SCAI_STRENGTH=<strength> ]" << std::endl;
        std::cout << "            [ --SCAI_VARIATION=<variation> ]" << std::endl;
        return -1;
    }

    std::cout << "Read D from "  << argv[1] << ", T from " << argv[2] 
              << ", So from " << argv[3] << ", hrz from " << argv[4] << std::endl;

    CSRSparseMatrix<double> D( argv[1] );
    DenseVector<double> T( argv[2] );
    DenseVector<double> So( argv[3]  );
    DenseVector<double> hrz( argv[4] );

    std::cout << "D = " << D << std::endl;
    std::cout << "hrz = " << hrz << std::endl;
    std::cout << "So = " << So << std::endl;
    std::cout << "T = " << T << std::endl;

    IndexType ny = hrz.size();
    IndexType n  = So.size();
    IndexType nz = n / ny ;

    SCAI_ASSERT_EQ_ERROR( ny * nz, n , "Illegal factors ny = " << ny << ", nz = " << nz )

    IndexType nray = D.getNumRows();

    SCAI_ASSERT_EQ_ERROR( D.getNumColumns(), n, "D must have #colums equal to problem size " << ny << " x " << nz )
    SCAI_ASSERT_EQ_ERROR( T.size(), D.getNumRows(), "T cannot be rhs for D" )

    int strength = 10;
    int variation = 2;

    common::Settings::getEnvironment( strength, "SCAI_STRENGTH" );
    common::Settings::getEnvironment( variation, "SCAI_VARIATION" );

    std::cout << "Use strength = " << strength << ", variation = " << variation << std::endl;
    CSRSparseMatrix<double> A;
    CSRSparseMatrix<double> L;

    setupSmoothMatrix( L, ny, nz, double( strength ) );

    joinMatrix( A, D, L );

    DenseVector<double> T_ext;

    zeroExtend( T_ext, T, L.getNumRows() );

    DenseVector<double> lb( So );
    DenseVector<double> ub( So );

    lb *= ( 1. - variation / 100. ); //for testing: 0.01;
    ub  *= ( 1. + variation / 100. ); //for testing: 100.0;
 
    for ( IndexType i = 0; i < ny; ++i )
    {
        for ( IndexType j = 0; j < nz; ++j )
        {
            if ( hrz[i] > j + 1 )
            {
                IndexType pos = i * nz + j;
                lb[pos] = So[pos] * 0.9999;
                ub[pos] = So[pos] * 1.0001;
            }
        }
    }

    A.writeToFile( argv[5] );
    T_ext.writeToFile( argv[6] );
    lb.writeToFile( argv[7] );
    ub.writeToFile( argv[8] );

    std::cout << "Written A to " << argv[5] << ", b to " << argv[6]
              << ", lb to " << argv[7] << ", ub to " << argv[8] << std::endl;
}
