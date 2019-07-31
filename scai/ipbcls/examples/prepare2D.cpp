/**
 * @file prepare2D.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Example of least square problem with boundary conditions
 * @author Thomas Brandes, Andreas Borgen Langva
 * @date 21.07.2017
 */

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>

using namespace scai;
using namespace lama;

typedef DefaultReal ValueType;

void setupSmoothMatrix( CSRSparseMatrix<ValueType>& L, const IndexType ny, const IndexType nz, const ValueType strength )
{
    // not neccesary, but to keep in accordance with the original

    common::Stencil2D<ValueType> stencil( 5 );

    common::Grid2D grid( ny, nz );

    grid.setBorderType( 0, common::BorderType::ABSORBING );
    grid.setBorderType( 1, common::BorderType::ABSORBING );

    StencilMatrix<ValueType> stencilMatrix( grid, stencil );

    L = stencilMatrix;

    L.scale( - strength / 4 );
}

void zeroExtend( DenseVector<ValueType>& T_ext,
                 const DenseVector<ValueType>& T, const IndexType nZeros )
{
    T_ext.allocate( T.size() + nZeros );

    hmemo::WriteAccess<ValueType> wT( T_ext.getLocalValues() );
    hmemo::ReadAccess<ValueType> rT( T.getLocalValues() );

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

    auto D = read<CSRSparseMatrix<ValueType>>( argv[1] );
    auto T = read<DenseVector<ValueType>>( argv[2] );
    auto So = read<DenseVector<ValueType>>( argv[3]  );
    auto hrz = read<DenseVector<ValueType>>( argv[4] );

    std::cout << "D = " << D << std::endl;
    std::cout << "hrz = " << hrz << std::endl;
    std::cout << "So = " << So << std::endl;
    std::cout << "T = " << T << std::endl;

    IndexType ny = hrz.size();
    IndexType n  = So.size();
    IndexType nz = n / ny ;

    SCAI_ASSERT_EQ_ERROR( ny * nz, n , "Illegal factors ny = " << ny << ", nz = " << nz )

    IndexType nray = D.getNumRows();
    std::cout << "#ray = " << nray << std::endl;

    SCAI_ASSERT_EQ_ERROR( D.getNumColumns(), n, "D must have #colums equal to problem size " << ny << " x " << nz )
    SCAI_ASSERT_EQ_ERROR( T.size(), D.getNumRows(), "T cannot be rhs for D" )

    int strength = 10;
    int variation = 10;

    common::Settings::getEnvironment( strength, "SCAI_STRENGTH" );
    common::Settings::getEnvironment( variation, "SCAI_VARIATION" );

    std::cout << "Use strength = " << strength << ", variation = " << variation << std::endl;
    CSRSparseMatrix<ValueType> A;
    CSRSparseMatrix<ValueType> L;

    setupSmoothMatrix( L, ny, nz, ValueType( strength ) );

    A.hcat( D, L );

    std::cout << "A = " << A << std::endl;

    DenseVector<ValueType> T_ext;

    zeroExtend( T_ext, T, L.getNumRows() );

    DenseVector<ValueType> lb( So );
    DenseVector<ValueType> ub( So );

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
