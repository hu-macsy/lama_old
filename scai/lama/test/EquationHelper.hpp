/**
 * @file EquationHelper.h
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
 * @brief EquationHelper.h
 * @author Alexander Büchel, Matthias Makulla
 * @date 09.02.2012
 * @since 1.0.0
 */

#ifndef LAMA_EQUATIONHELPER_H_
#define LAMA_EQUATIONHELPER_H_

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>

using namespace lama;

class EquationHelper
{
public:

    template<typename ValueType>
    struct EquationSystem
    {
        CSRSparseMatrix<ValueType> coefficients;
        DenseVector<ValueType> rhs;
        DenseVector<ValueType> solution;
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
    static const std::vector<EquationSystem<T> ( * )()>& getFunctions();
};

template<typename T>
const std::vector<EquationHelper::EquationSystem<T> ( * )()>& EquationHelper::getFunctions()
{
    typedef EquationSystem<T> ( *t_system )();
    static const t_system t_arr[] =
    {
        get3x3SystemA<T>, get4x4SystemA<T>, get8x8SystemA<T>, get4x4SystemB<T>, get8x8SystemB<T>, get8x8SystemC<T>,
        get8x8SystemD<T>, get8x8EmptyDiagonal<T>
    };
    static std::vector<t_system> t_vec;

    if ( t_vec.empty() )
        for ( std::size_t i = 0; i < sizeof( t_arr ) / sizeof( t_system ); ++i )
        {
            t_vec.push_back( t_arr[i] );
        }

    return t_vec;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get3x3SystemA()
{
    IndexType dim = 3;
    ValueType coefficientValues[] =
    { 9.0f, 1.0f, 2.0f, -2.0f, 7.0f, 1.0f, -1.0f, 1.0f, 8.0f };
    ValueType rhsValues[] =
    { 9.0f, 11.0f, -7.0f };
    ValueType solutionValues[] =
    { 1.0f, 2.0f, -1.0f };
    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setRawDenseData( dim, dim, coefficientValues );
    DenseVector<ValueType> rhs( dim, rhsValues );
    DenseVector<ValueType> solution( dim, solutionValues );
    EquationSystem<ValueType> equationSystem =
    { coefficients, rhs, solution };
    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get4x4SystemA()
{
    IndexType dim = 4;
    ValueType coefficientValues[] =
    { 1.0f, 2.0f, 3.0f, 2.0f, 2.0f, 4.0f, -1.0f, -1.0f, 3.0f, -1.0f, 2.0f, 7.0f, 2.0f, -1.0f, 7.0f, -4.0f };
    ValueType rhsValues[] =
    { 8.0f, 13.0f, 12.0f, -16.0f };
    ValueType solutionValues[] =
    { 1.0f, 3.0f, -1.0f, 2.0f };
    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setRawDenseData( dim, dim, coefficientValues );
    DenseVector<ValueType> rhs( dim, rhsValues );
    DenseVector<ValueType> solution( dim, solutionValues );
    EquationSystem<ValueType> equationSystem =
    { coefficients, rhs, solution };
    return equationSystem;
}

template<typename ValueType>
EquationHelper::EquationSystem<ValueType> EquationHelper::get8x8SystemA()
{
    IndexType dim = 8;
    ValueType matrixValues[] =
    {
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        2.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f
    };
    ValueType rhsValues[] =
    { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };
    ValueType solutionValues[] =
    { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    CSRSparseMatrix<ValueType> coefficients;
    coefficients.setRawDenseData( dim, dim, matrixValues );
    DenseVector<ValueType> rhs( dim, rhsValues );
    DenseVector<ValueType> solution( dim, solutionValues );
    EquationSystem<ValueType> equationSystem =
    { coefficients, rhs, solution };
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
    ValueType rhsValues[] =
    { 28.0f, 22.0f, 18.0f, 16.0f, 16.0f, 18.0f, 22.0f, 28.0f };
    ValueType solutionValues[] =
    { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    CSRSparseMatrix<ValueType> coefficients( dim, dim, values );
    DenseVector<ValueType> rhs( dim, rhsValues );
    DenseVector<ValueType> solution( dim, solutionValues );
    EquationSystem<ValueType> equationSystem =
    { coefficients, rhs, solution };
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
    CSRSparseMatrix<ValueType> coefficients( dim, dim, values );
    DenseVector<ValueType> rhs( dim, rhsValues );
    DenseVector<ValueType> solution( dim, solutionValues );
    EquationSystem<ValueType> equationSystem =
    { coefficients, rhs, solution };
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
    ValueType rhsValues[] =
    { 28.0f, 22.0f, 18.0f, 28.0f, 16.0f, 16.0f, 18.0f, 22.0f };
    ValueType solutionValues[] =
    { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    CSRSparseMatrix<ValueType> coefficients( dim, dim, values );
    DenseVector<ValueType> rhs( dim, rhsValues );
    DenseVector<ValueType> solution( dim, solutionValues );
    EquationSystem<ValueType> equationSystem =
    { coefficients, rhs, solution };
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
    CSRSparseMatrix<ValueType> coefficients( dim, dim, values );
    DenseVector<ValueType> rhs( dim, rhsValues );
    DenseVector<ValueType> solution( dim, solutionValues );
    EquationSystem<ValueType> equationSystem =
    { coefficients, rhs, solution };
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
    // (0, 1, 6, 7, 4, 5, 3, 2)
    ValueType rhsValues[] =
    { 28.0f, 22.0f, 18.0f, 28.0f, 16.0f, 16.0f, 18.0f, 22.0f };
    ValueType solutionValues[] =
    { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
    CSRSparseMatrix<ValueType> coefficients( dim, dim, values );
    DenseVector<ValueType> rhs( dim, rhsValues );
    DenseVector<ValueType> solution( dim, solutionValues );
    EquationSystem<ValueType> equationSystem =
    { coefficients, rhs, solution };
    return equationSystem;
}

#endif // LAMA_EQUATIONHELPER_H_
