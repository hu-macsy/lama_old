/**
 * @file BenchOps.cpp
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
 * @brief Benchmarking of binary and UnaryOp operations
 * @author Thomas Brandes
 * @date 21.10.2016
 */

#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

using namespace std;
using namespace scai;
using namespace hmemo;
using namespace utilskernel;
using common::BinaryOp;
using common::UnaryOp;

int main( int argc, const char* argv[] )
{
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    typedef DefaultReal ValueType;

    const IndexType N = 10 * 1000 * 1000;

    HArray<ValueType> values1( N, ctx );
    HArray<ValueType> values2( ctx );
    HArray<ValueType> values3( ctx );

    values1.resize( N );
    values2.resize( N );

    HArrayUtils::setRandom( values1, 2 );
    HArrayUtils::setRandom( values2, 2 );

    // move random values in range 1.0 - 3.0 instead of 0 .. 2 

    HArrayUtils::setScalar( values1, ValueType( 1 ), BinaryOp::ADD );
    HArrayUtils::setScalar( values2, ValueType( 1 ), BinaryOp::ADD );

    for ( int i = 0; i < static_cast<int>( UnaryOp::MAX_UNARY_OP ); ++i )
    {
        UnaryOp op = UnaryOp( i );

        double start = common::Walltime::get();

        for ( int iter = 0; iter < 10; ++iter )
        {
            HArrayUtils::unaryOp( values3, values1, op, ctx );
        }

        double time1 = common::Walltime::get() - start;

        cout << "Time for UnaryOp op = " << op << ": " << time1 << " seconds." << endl;
    }

    for ( int i = 0; i < static_cast<int>( BinaryOp::MAX_BINARY_OP ); ++i )
    {
        BinaryOp op = BinaryOp( i );

        double start = common::Walltime::get();

        for ( int iter = 0; iter < 10; ++iter )
        {
            HArrayUtils::binaryOp( values3, values1, values2, op, ctx );
        }

        double time1 = common::Walltime::get() - start;

        cout << "Time for binary op = " << op << ": " << time1 << " seconds." << endl;
    }

    for ( int i = 0; i < static_cast<int>( BinaryOp::MAX_BINARY_OP ); ++i )
    {
        ValueType sum = 0;

        BinaryOp op = BinaryOp( i );

        if ( op != BinaryOp::ADD && op != BinaryOp::MAX && op != BinaryOp::MIN && op != BinaryOp::ABS_MAX )
        {
            continue;
        }

        double start = common::Walltime::get();

        for ( int iter = 0; iter < 10; ++iter )
        {
            sum = HArrayUtils::reduce( values1, op, ctx );
        }

        double time1 = common::Walltime::get() - start;

        cout << "Time for reduce op = " << op << ": " << time1 << " seconds. Result = " << sum << endl;
    }
}
