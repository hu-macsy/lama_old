/**
 * @file EquationHelper.hpp
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
 * @brief Definition of class that provides some equation systems for tests.
 * @author Matthias Makulla
 * @date 09.02.2012
 */

#pragma once

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>

namespace scai
{

class EquationHelper
{
public:

    template<typename ValueType>
    struct EquationSystem
    {
        lama::CSRSparseMatrix<ValueType> coefficients;
        lama::DenseVector<ValueType> rhs;
        lama::DenseVector<ValueType> solution;
    };

    typedef EquationSystem<float> ( *ssystem )();
    typedef EquationSystem<double> ( *dsystem )();

    template<typename ValueType>
    static EquationSystem<ValueType> get3x3SystemA();
    template<typename ValueType>
    static EquationSystem<ValueType> get4x4SystemA();
    template<typename ValueType>
    static EquationSystem<ValueType> get8x8SystemA();
    template<typename ValueType>
    static EquationSystem<ValueType> get4x4SystemB();
    template<typename ValueType>
    static EquationSystem<ValueType> get8x8SystemB();
    template<typename ValueType>
    static EquationSystem<ValueType> get8x8SystemC();
    template<typename ValueType>
    static EquationSystem<ValueType> get8x8SystemD();
    template<typename ValueType>
    static EquationSystem<ValueType> get8x8EmptyDiagonal();

    template<typename T>
    static const std::vector<EquationSystem<T> (* )()>& getFunctions();
};

template<typename T>
const std::vector<EquationHelper::EquationSystem<T> (* )()>& EquationHelper::getFunctions()
{
    typedef EquationSystem<T> ( *t_system )();
    static const t_system t_arr[] =
    {
        get3x3SystemA<T>, get4x4SystemA<T>, get8x8SystemA<T>, get4x4SystemB<T>, get8x8SystemB<T>, get8x8SystemC<T>,
        get8x8SystemD<T>, get8x8EmptyDiagonal<T>
    };
    static std::vector<t_system> t_vec;

    if ( t_vec.empty() )
    {
        for ( std::size_t i = 0; i < sizeof( t_arr ) / sizeof( t_system ); ++i )
        {
            t_vec.push_back( t_arr[i] );
        }
    }

    return t_vec;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get3x3SystemA()
{
    IndexType dim = 3;

    ValueType coefficientValues[] = { 9.0f, 1.0f, 2.0f, -2.0f, 7.0f, 1.0f, -1.0f, 1.0f, 8.0f };
    ValueType rhsValues[]         = { 9.0f, 11.0f, -7.0f };
    ValueType solutionValues[]    = { 1.0f, 2.0f, -1.0f };

    EquationSystem<ValueType> equationSystem;

    equationSystem.coefficients.setRawDenseData( dim, dim, coefficientValues );
    equationSystem.rhs.setRawData( dim, rhsValues );
    equationSystem.solution.setRawData( dim, solutionValues );

    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get4x4SystemA()
{
    IndexType dim = 4;

    ValueType values[] =
    { 1.0f, 2.0f, 3.0f, 2.0f, 2.0f, 4.0f, -1.0f, -1.0f, 3.0f, -1.0f, 2.0f, 7.0f, 2.0f, -1.0f, 7.0f, -4.0f };
    ValueType rhsValues[] =
    { 8.0f, 13.0f, 12.0f, -16.0f };
    ValueType solutionValues[] =
    { 1.0f, 3.0f, -1.0f, 2.0f };

    SCAI_ASSERT_EQ_ERROR( dim * dim, sizeof( values ) / sizeof( ValueType ), "serious size mismatch" )

    EquationSystem<ValueType> equationSystem;

    equationSystem.coefficients.setRawDenseData( dim, dim, values );
    equationSystem.rhs.setRawData( dim, rhsValues );
    equationSystem.solution.setRawData( dim, solutionValues );

    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get8x8SystemA()
{
    IndexType dim = 8;

    ValueType values[] =
    {
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f
    };

    SCAI_ASSERT_EQ_ERROR( dim * dim, sizeof( values ) / sizeof( ValueType ), "serious size mismatch" )

    ValueType rhsValues[]      = { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };
    ValueType solutionValues[] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };

    EquationSystem<ValueType> equationSystem;

    equationSystem.coefficients.setRawDenseData( dim, dim, values );
    equationSystem.rhs.setRawData( dim, rhsValues );
    equationSystem.solution.setRawData( dim, solutionValues );

    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get8x8EmptyDiagonal()
{
    IndexType dim = 8;
    ValueType values[] =
    {
        0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 2.0f, 1.0f, 0.0f,
        1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f,
        2.0f, 3.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 7.0f,
        6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f
    };

    SCAI_ASSERT_EQ_ERROR( dim * dim, sizeof( values ) / sizeof( ValueType ), "serious size mismatch" )

    ValueType rhsValues[] = { 28.0f, 22.0f, 18.0f, 16.0f, 16.0f, 18.0f, 22.0f, 28.0f };
    ValueType solutionValues[] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };

    EquationSystem<ValueType> equationSystem;

    equationSystem.coefficients.setRawDenseData( dim, dim, values );
    equationSystem.rhs.setRawData( dim, rhsValues );
    equationSystem.solution.setRawData( dim, solutionValues );

    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get4x4SystemB()
{
    IndexType dim = 4;
    ValueType values[] =
    { 4.0f, 1.0f, 2.0f, 1.0f, 2.0f, 3.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 2.0f, 1.0f, 2.0f, 1.0f };
    ValueType rhsValues[] =
    { 8.0f, 8.0f, 6.0f, 6.0f };
    ValueType solutionValues[] =
    { 1.0f, 1.0f, 1.0f, 1.0f };

    EquationSystem<ValueType> equationSystem;

    equationSystem.coefficients.setRawDenseData( dim, dim, values );
    equationSystem.rhs.setRawData( dim, rhsValues );
    equationSystem.solution.setRawData( dim, solutionValues );

    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get8x8SystemB()
{
    IndexType dim = 8;
    ValueType values[] =
    {
        7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, // 0
        6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, // 1
        5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, // 2
        0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, // 3
        3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, // 6
        4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, // 7
        2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, // 5
        1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f // 4
    };

    SCAI_ASSERT_EQ_ERROR( dim * dim, sizeof( values ) / sizeof( ValueType ), "serious size mismatch" )

    ValueType rhsValues[] =
    { 28.0f, 22.0f, 18.0f, 28.0f, 16.0f, 16.0f, 18.0f, 22.0f };
    ValueType solutionValues[] =
    { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };

    EquationSystem<ValueType> equationSystem;

    equationSystem.coefficients.setRawDenseData( dim, dim, values );
    equationSystem.rhs.setRawData( dim, rhsValues );
    equationSystem.solution.setRawData( dim, solutionValues );

    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get8x8SystemC()
{
    IndexType dim = 8;
    ValueType values[] =
    {
        7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 1.0f, 0.0f, 1.0f,
        2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f,
        3.0f, 4.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 6.0f,
        5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f
    };
    // (0, 7, 2, 3, 4, 5, 6, 1)
    ValueType rhsValues[] =
    { 28.0f, 28.0f, 22.0f, 18.0f, 16.0f, 16.0f, 18.0f, 22.0f };
    ValueType solutionValues[] =
    { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };

    EquationSystem<ValueType> equationSystem;

    equationSystem.coefficients.setRawDenseData( dim, dim, values );
    equationSystem.rhs.setRawData( dim, rhsValues );
    equationSystem.solution.setRawData( dim, solutionValues );

    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get8x8SystemD()
{
    IndexType dim = 8;

    ValueType values[] =
    {
        7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, // 0
        6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, // 1
        5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, // 6
        0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, // 7
        3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, // 4
        4.0f, 3.0f, 2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, // 5
        2.0f, 1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, // 3
        1.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f // 2
    };

    ValueType rhsValues[] =      { 28.0f, 22.0f, 18.0f, 28.0f, 16.0f, 16.0f, 18.0f, 22.0f };
    ValueType solutionValues[] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };

    EquationSystem<ValueType> equationSystem;

    equationSystem.coefficients.setRawDenseData( dim, dim, values );
    equationSystem.rhs.setRawData( dim, rhsValues );
    equationSystem.solution.setRawData( dim, solutionValues );

    return equationSystem;
}

}
