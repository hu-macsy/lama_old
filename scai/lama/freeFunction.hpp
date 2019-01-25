/**
 * @file lama/freeFunction.hpp
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
 * @brief Definiton of free functions to construct objects in lama namespace.
 * @author Thomas Brandes, Andreas Borgen Longva
 * @date 25.01.18
 */
#pragma once

#include <scai/hmemo/Context.hpp>
#include <scai/dmemo/Distribution.hpp>
#include <scai/lama/_Vector.hpp>
#include <scai/lama/matrix/_Matrix.hpp>
#include <scai/lama/storage/_MatrixStorage.hpp>

namespace scai
{

namespace lama
{

/** 
 *  @brief Function that returns a vector/matrix by evaluation of an expression
 * 
 *  @tparam    ObjectType can be any derived Vector or Matrix class
 *  @param[in] exp        is a vector or matrix expression
 *  @param[in] ctx        context for the new created object
 *  @returns              a new object, distribution is implicitly given by the expression
 *
 *  \code
 *     const auto vector = eval<DenseVector<double>>( a, b );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind that the function requires the result type as template argument.
 */
template<typename ObjectType, typename expression>
ObjectType evalDEP( expression exp, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    ObjectType obj( ctx );
    obj = exp;
    return obj;
}

/** 
 *  @brief Function that returns a zero vector or array of a given size 
 * 
 *  @tparam    ObjectType can be any derived MatrixStorage or Matrix class.
 *  @param[in] size       number of entries
 *  @param[in] ctx        context where operations on matrix will be execued
 *  @returns              a zero vector of given size
 *
 *  \code
 *     const auto v = zero<DenseVector<double>>( n );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind the function requires the result type as template argument.
 */
template<typename ObjectType>
ObjectType zeroDEP( const IndexType size, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    static_assert(std::is_base_of<_Vector, ObjectType>::value, "zero<ObjectType>( n [,ctx] ): ObjectType is not a vector class");
    ObjectType obj( ctx );
    typename ObjectType::ObjectValueType zeroVal = 0;
    obj.setSameValue( size, zeroVal );
    return obj;
}

template<typename ObjectType>
ObjectType zeroDEP( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    // ObjectType must be a vector class
    static_assert(std::is_base_of<_Vector, ObjectType>::value, "zero<ObjectType>( dist [,ctx] ): ObjectType is not a vector class");
    ObjectType obj( ctx );
    typename ObjectType::ObjectValueType zeroVal = 0;
    obj.setSameValue( dist, zeroVal );
    return obj;
}

template<typename ObjectType>
ObjectType undefinedDEP( const IndexType size, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    // ObjectType must be a vector class

    static_assert(std::is_base_of<_Vector, ObjectType>::value, "undefined<ObjectType>( n [,ctx] ): ObjectType is not a vector class");
    ObjectType obj( ctx );
    obj.allocate( size );
    return obj;
}

template<typename ObjectType>
ObjectType undefinedDEP( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    // ObjectType must be a vector class

    static_assert(std::is_base_of<_Vector, ObjectType>::value, "undefined<ObjectType>( dist [,ctx] ): ObjectType is not a vector class");
    ObjectType obj( ctx );
    obj.allocate( dist );
    return obj;
}

/** 
 *  @brief Function that returns a zero matrix or storage for a given size 
 * 
 *  @tparam    ObjectType can be any derived MatrixStorage or Matrix class.
 *  @param[in] numRows    number of rows
 *  @param[in] numColumns number of columns
 *  @param[in] ctx        context where operations on matrix will be execued
 *  @returns              a zero matrix of given size
 *
 *  \code
 *     const auto storage = zero<CSRStorage<double>>( m, n );
 *     const auto matrix  = zero<CSRSparseMatrix<double>>( m, n );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind the function requires the result type as template argument.
 */
template<typename ObjectType>
ObjectType zero( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    static_assert(std::is_base_of<_Matrix, ObjectType>::value || std::is_base_of<_MatrixStorage, ObjectType>::value, 
                  "zero<ObjectType>( numRows, numCols [, ctx] ): ObjectType is not a matrix/storage class");
    ObjectType obj( ctx );
    obj.allocate( numRows, numColumns );
    return obj;
}

/** 
 *  @brief Function that returns a zero matrix with a given distribution
 * 
 *  @tparam    ObjectType can be any derived Matrix class.
 *  @param[in] rowDist    distribution of rows
 *  @param[in] colDist    distribution of columns
 *  @param[in] ctx        context where operations on matrix will be execued
 *  @returns              a zero matrix of given size
 *
 *  \code
 *     const auto rowDist = std::make_shared<BlockDistribution>( m, comm );
 *     const auto colDist = std::make_shared<BlockDistribution>( n, comm );
 *     const auto storage = zero<CSRStorage<double>>( rowDist, colDist );
 *     const auto matrix  = zero<CSRSparseMatrix<double>>( rowDist, colDist );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind the function requires the result type as template argument.
 */
template<typename ObjectType>
ObjectType zero( 
    dmemo::DistributionPtr rowDist,
    dmemo::DistributionPtr colDist,
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr()  )
{
    static_assert(std::is_base_of<_Matrix, ObjectType>::value, "zero<ObjectType>( rowDist, colDist [, ctx] ): ObjectType is not a matrix class");
    ObjectType obj( ctx );
    obj.allocate( rowDist, colDist );
    return obj;
}

/** 
 *  @brief Function that returns a identy matrix or storage for a given size 
 * 
 *  @tparam    ObjectType can be any derived MatrixStorage or Matrix class.
 *  @param[in] size       number of rows and columns
 *  @param[in] ctx        context where operations on matrix will be execued
 *  @returns              a zero matrix of given size
 *
 *  \code
 *     const auto storage = identity<CSRStorage<double>>( n );
 *     const auto matrix  = identity<CSRSparseMatrix<double>>( n );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind the function requires the result type as template argument.
 */
template<typename ObjectType>
ObjectType identity( const IndexType size, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    ObjectType obj( ctx );
    obj.setIdentity( size );
    return obj;
}

/** 
 *  @brief Function that returns a identy matrix for a given distributon
 * 
 *  @tparam    ObjectType can be any derived MatrixStorage or Matrix class.
 *  @param[in] dist       row and column distribution (must be the same)
 *  @param[in] ctx        context where operations on matrix will be execued
 *  @returns              a zero matrix of given size
 *
 *  \code
 *     const auto dist   = std::make_shared<BlockDistribution>( n, comm );
 *     const auto matrix = identity<CSRSparseMatrix<double>>( dist );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind the function requires the result type as template argument.
 */
template<typename ObjectType>
ObjectType identity( dmemo::DistributionPtr dist, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    ObjectType obj( ctx );
    obj.setIdentity( dist );
    return obj;
}

/** 
 *  @brief Function that returns an object (storage, matrix, array) read in from a file.
 * 
 *  @tparam    ObjectType can be any derived MatrixStorage or Matrix class.
 *  @param[in] fileName   specifies the name of the file from which the object is read.
 *  @param[in] ctx        context where operations on matrix will be execued
 *  @returns              a zero matrix of given size
 *
 *  \code
 *     const auto matrix = read<CSRSparseMatrix<double>>( "Emilia_923.mtx" );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind the function requires the result type as template argument.
 */
template<typename ObjectType>
ObjectType read( const std::string& fileName, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    ObjectType obj( ctx );
    obj.readFromFile( fileName );
    return obj;
}

template<typename ObjectType>
ObjectType read( const std::string& fileName, 
                 dmemo::CommunicatorPtr comm,
                 hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    ObjectType obj( ctx );
    obj.readFromFile( fileName, comm );
    return obj;
}

/** 
 *  @brief Function that returns a vector of a given size initialized with the same value.
 * 
 *  @tparam    ObjectType can be any derived vector class.
 *  @param[in] n          specifies the size of the vector                             
 *  @param[in] value      is the vector assigned to all elements of the vector
 *  @param[in] ctx        Context that is used for the filling and the generated vector
 *  @returns              a new object of type ObjectType
 *
 *  \code
 *     const auto v1 = fill<SparseVector<double>>( n, 10 );
 *     const auto v2 = fill<DenseVector<double>>( n, 10 );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind the function requires the result type as template argument.
 *
 *  Note: DEPRECATED function, please use sparseVectorFill or denseVectorFill
 */
template<typename ObjectType>
ObjectType fillDEP( const IndexType n, typename ObjectType::ObjectValueType value, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    ObjectType obj( ctx );
    obj.setSameValue( n, value );
    return obj;
}

/** 
 *  @brief Function that returns a vector for a given distribution initialized with the same value.
 * 
 *  @tparam    ObjectType can be any derived vector class.
 *  @param[in] dist       specifies the distribution of the new vector                             
 *  @param[in] value      is the vector assigned to all elements of the vector
 *  @param[in] ctx        Context that is used for the filling and the generated vector
 *  @returns              a new object of type ObjectType
 *
 *  \code
 *     const auto dist = std::make_shared<BlockDistribution>( n );
 *     const auto v1 = fill<SparseVector<double>>( dist, 10 );
 *     const auto v2 = fill<DenseVector<double>>( dist, 10 );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  Please keep in mind the function requires the result type as template argument.
 */
template<typename ObjectType>
ObjectType fillDEP( const dmemo::DistributionPtr dist, typename ObjectType::ObjectValueType value, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    ObjectType obj( ctx );
    obj.setSameValue( dist, value );
    return obj;
}

template<typename TargetType, typename SourceType>
TargetType convert( const SourceType& input, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    TargetType obj( ctx );
    obj.assign( input );
    return obj;
}

template<typename TargetType, typename ValueType>
TargetType convertRawCSR(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType* const ia,
    const IndexType* const ja,
    const ValueType* const values,
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    // wrap the pointer data into LAMA arrays ( without copies )
    hmemo::HArrayRef<IndexType> csrIA( numRows + 1, ia );
    hmemo::HArrayRef<IndexType> csrJA( numValues, ja );
    hmemo::HArrayRef<ValueType> csrValues( numValues, values );
    // now set the data on the context of this storage via virtual method
    TargetType obj( ctx );
    obj.setCSRData( numRows, numColumns, csrIA, csrJA, csrValues );
    return obj;
}

/** 
 *  @brief Function that returns a copy of a (sparse, dense) vector with a new distribution
 * 
 *  @tparam    ObjectType can be any derived vector class.
 *  @param[in] input      is the vector that is copied
 *  @param[in] dist       is the distribution of the new vector
 *  @param[in] ctx        is the context of the new vector
 *  @returns              a new object of type ObjectType
 *
 *  \code
 *     auto repV  = read<DenseVector<double>>( "vector.txt" );
 *     auto distV = distribute( repV, std::make_shared<BlockDistribution>( repV.size() ), repV.getContextPtr() );
 *  \endcode
 *
 *  This template function can be considered as a syntactical help for code abbreviation.
 *  In contrary to the other routines the template argument and return type can be deduced from the input vector.
 */
template<typename DistObjType, typename GlobalObjType>
DistObjType distribute( const GlobalObjType& input, dmemo::DistributionPtr dist, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    DistObjType obj( ctx );
    obj.assignDistribute( input, dist );
    return obj;
}

/**
 *  @brief Resize is nearly the same as distribute but here the new distributon might have a different global size
 *
 *  Note: The input vector might be either truncated or filled up with zero.
 */
template<typename DistObjType, typename GlobalObjType>
DistObjType resize( 
    const GlobalObjType& input, 
    dmemo::DistributionPtr dist, 
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    DistObjType obj( ctx );
    obj.assign( input );
    obj.resize( dist );
    return obj;
}

/** 
 *  @brief Function that (re-)distributes a global matrix/storage 
 * 
 *  @tparam    DistObjType   type of the distributed object, should be a matrix class
 *  @tparam    GlobalObjType type of the input storage/matrix that will be distributed
 *  @param[in] input         is the matrix/storage that will be distributed
 *  @param[in] rowDist       is the row distribution of the new matrix
 *  @param[in] colDist       is the column distribution of the new matrix
 *  @param[in] ctx           is the context of the new matrix
 *  @returns                 a copy of the orginal matrix that is now distributed as specified
 *
 *  The following relations must hold:
 *
 *   - input.getNumRows() == rowDist->getGlboalSize()
 *   - input.getNumColumns() == colDist->getGlboalSize()
 *
 *  \code
 *      auto globalStorage = read<CSRStorage<double>>( "input.mtx" );
 *      auto rowDist = std::make_shared<BlockDistribution>( globalStorage.getNumRows() );
 *      auto colDist = std::make_shared<BlockDistribution>( globalStorage.getNumColumns() );
 *      auto matrix  = distribute<CSRSparseMatrix<double>>( globalStorage, rowDist, colDist );
 *      auto rowDist1 = std::make_shared<CyclicDistribution>( globalStorage.getNumRows() );
 *      auto colDist1 = std::make_shared<CyclicDistribution>( globalStorage.getNumColumns() );
 *      auto matrix1  = distribute<CSRSparseMatrix<double>>( matrix, rowDist1, colDist1 );
 *  \endcode
 *
 *  Please keep in mind that for a matrix the row distribution specifies the distribution
 *  of the target vector and the column distribution the distribution of the source vector.
 */
template<typename DistObjType, typename GlobalObjType>
DistObjType distribute( 
    const GlobalObjType& input, 
    dmemo::DistributionPtr rowDist, 
    dmemo::DistributionPtr colDist, 
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    DistObjType obj( ctx );
    obj.assignDistribute( input, rowDist, colDist );
    return obj;
}

/**
 *  @brief Make a copy of an existing storage/matrix and resize it with new distributions.
 */
template<typename DistObjType, typename GlobalObjType>
DistObjType resize( 
    const GlobalObjType& input, 
    dmemo::DistributionPtr rowDist, 
    dmemo::DistributionPtr colDist, 
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    DistObjType obj( ctx );
    obj.assign( input );
    obj.resize( rowDist, colDist );
    return obj;
}

/** 
 *  @brief Function that returns a matrix where only the diagonal is set.
 */
template<typename MatrixType, typename ArrayType>
MatrixType diagonal(
    const ArrayType& diag,
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() )
{
    MatrixType obj( ctx );
    obj.assignDiagonal( diag );
    return obj;
}

} /* end namespace lama */

} /* end namespace scai */

