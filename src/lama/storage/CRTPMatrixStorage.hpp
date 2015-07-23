/**
 * @file CRTPMatrixStorage.hpp
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
 * @brief Definition of a common base class for all matrix storage formats.
 * @author Thomas Brandes
 * @date 27.04.2011
 * @since 1.0.0
 */

#ifndef LAMA_CRTP_MATRIX_STORAGE_HPP_
#define LAMA_CRTP_MATRIX_STORAGE_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <lama/storage/MatrixStorage.hpp>

namespace lama
{

/** This template class supports static polymorphism to define
 *  common routines for derived classes of MatrixStorage.
 *
 *  Therefore it uses the Curiously Recurring Template Pattern (CRTP)
 *  as a C++ idiom where a class X derived from this template is
 *  itself a template argument.
 *
 *  @tparam Derived is the derived class of MatrixStorage<ValueType>
 *  @tparam ValueType specifies the type of the matrix valuess
 *
 *  @todo: create, copy should also be defined here
 */

template<class Derived,typename ValueType>
class COMMON_DLL_IMPORTEXPORT CRTPMatrixStorage: public MatrixStorage<ValueType>
{
public:

    /** Constructor for a zero-matrix of a certain size. */

    CRTPMatrixStorage( const IndexType numRows, const IndexType numColumns )

        : MatrixStorage<ValueType>( numRows, numColumns )

    {
    }

    /** Constructor for an empty matrix. */

    CRTPMatrixStorage()

        : MatrixStorage<ValueType>( 0, 0 )

    {
    }

    /** Implementation of _MatrixStorage::setCSRData for derived class.
     *
     *  This routine requires that the derived class provides a corresponding
     *  routine setCSRDataImpl that can deal with a typed version of values.
     */

    void setCSRData(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValues,
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const _LAMAArray& values )
    {
        memory::ScalarType arrayType = values.getValueType();

        switch( arrayType )
        {

#define LAMA_SET_CSR_CALL( z, I, _ )                                                            \
case common::scalar::SCALAR_ARITHMETIC_TYPE##I:                                                 \
{                                                                                       \
    const LAMAArray<ARITHMETIC_TYPE##I>& typedValues =                                  \
            dynamic_cast<const LAMAArray<ARITHMETIC_TYPE##I>&>( values );                   \
    static_cast<Derived*>( this )->setCSRDataImpl(                                      \
            numRows, numColumns, numValues, ia, ja, typedValues, this->getContextPtr() );   \
    break;                                                                              \
}                                                                                       \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_SET_CSR_CALL, _ )

#undef LAMA_SET_CSR_CALL

default            :
{
    COMMON_THROWEXCEPTION( *this << ": setCSRData with value type " << arrayType << " not supported" )
}
}
}

/** Implementation for _MatrixStorage::buildCSRSizes */

void buildCSRSizes( LAMAArray<IndexType>& ia ) const
{
// The sizes will be available via buildCSR with NULL for ja, values

LAMAArray<IndexType>* ja = NULL;
LAMAArray<ValueType>* values = NULL;

static_cast<const Derived*>( this )->buildCSR( ia, ja, values, this->getContextPtr() );
}

void buildCSRData( LAMAArray<IndexType>& csrIA, LAMAArray<IndexType>& csrJA, _LAMAArray& csrValues ) const
{
memory::ScalarType arrayType = csrValues.getValueType();

switch ( arrayType )
{
#define LAMA_BUILD_CSR_CALL( z, I, _ )                                        \
case common::scalar::SCALAR_ARITHMETIC_TYPE##I:                               \
{                                                                     \
    LAMAArray<ARITHMETIC_TYPE##I>& typedValues =                      \
            dynamic_cast<LAMAArray<ARITHMETIC_TYPE##I>&>( csrValues );    \
    static_cast<const Derived*>( this )->buildCSR(                    \
            csrIA, &csrJA, &typedValues, this->getContextPtr() );         \
    break;                                                            \
}                                                                     \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_BUILD_CSR_CALL, _ )

#undef LAMA_BUILD_CSR_CALL

default:
{
    COMMON_THROWEXCEPTION( *this << ": build CSR with value type " << arrayType << " not supported" )
}

}

}

/** Get the i-th row of a storage as LAMA array. */

void getRow( _LAMAArray& row, const IndexType irow ) const
{
common::scalar::ScalarType arrayType = row.getValueType();

switch ( arrayType )
{
#define LAMA_GET_ROW_CALL( z, I, _ )                                           \
case common::scalar::SCALAR_ARITHMETIC_TYPE##I:                                \
{                                                                      \
    LAMAArray<ARITHMETIC_TYPE##I>& typedRow =                          \
            dynamic_cast<LAMAArray<ARITHMETIC_TYPE##I>&>( row );           \
    static_cast<const Derived*>( this )->getRowImpl( typedRow, irow ); \
    break;                                                             \
}                                                                      \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_GET_ROW_CALL, _ )

#undef LAMA_GET_ROW_CALL

default:
{
    COMMON_THROWEXCEPTION( "getRow for array of type " << arrayType << " not supported" )
}

}
}

void getDiagonal( _LAMAArray& diagonal ) const
{
if ( !this->hasDiagonalProperty() )
{
COMMON_THROWEXCEPTION( *this << ": has not diagonal property, cannot get diagonal" )
}

memory::ScalarType arrayType = diagonal.getValueType();

switch ( arrayType )
{
#define LAMA_GET_DIAGONAL_CALL( z, I, _ )                                           \
case common::scalar::SCALAR_ARITHMETIC_TYPE##I:                                     \
{                                                                           \
    LAMAArray<ARITHMETIC_TYPE##I>& typedDiagonal =                          \
            dynamic_cast<LAMAArray<ARITHMETIC_TYPE##I>&>( diagonal );           \
    static_cast<const Derived*>( this )->getDiagonalImpl( typedDiagonal );  \
    break;                                                             \
}                                                                           \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_GET_DIAGONAL_CALL, _ )

#undef LAMA_GET_DIAGONAL_CALL

default:
    COMMON_THROWEXCEPTION( "getDiagonal for array of type " << arrayType << " not supported" )

}
}

void setDiagonal( const Scalar value )
{
if ( !this->hasDiagonalProperty() )
{
COMMON_THROWEXCEPTION( *this << ": has not diagonal property, cannot set diagonal" )
}

static_cast<Derived*>( this )->setDiagonalImpl( value );
}

void setDiagonal( const _LAMAArray& diagonal )
{
IndexType numDiagonalElements = diagonal.size();

if ( numDiagonalElements > this->getNumRows() || numDiagonalElements > this->getNumColumns() )
{
COMMON_THROWEXCEPTION( "Diagonal of size " << numDiagonalElements << " too large for matrix: " << *this )
}

if ( !this->hasDiagonalProperty() )
{
COMMON_THROWEXCEPTION( *this << ": has not diagonal property, cannot set diagonal" )
}

memory::ScalarType arrayType = diagonal.getValueType();

switch ( arrayType )
{
#define LAMA_SET_DIAGONAL_CALL( z, I, _ )                                           \
case common::scalar::SCALAR_ARITHMETIC_TYPE##I:                                     \
{                                                                           \
    const LAMAArray<ARITHMETIC_TYPE##I>& typedDiagonal =                    \
            dynamic_cast<const LAMAArray<ARITHMETIC_TYPE##I>&>( diagonal );     \
    static_cast<Derived*>( this )->setDiagonalImpl( typedDiagonal );        \
    break;                                                                  \
}                                                                           \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_SET_DIAGONAL_CALL, _ )

#undef LAMA_SET_DIAGONAL_CALL

default:
    COMMON_THROWEXCEPTION( "setDiagonal to array of type " << arrayType << " not supported" )
}
}

void scale( const Scalar value )
{
static_cast<Derived*>( this )->scaleImpl( value );
}

/** Polymorph implementation for MatrixStorage<ValueType>::scale */

void scale( const _LAMAArray& diagonal )
{
LAMA_ASSERT_EQUAL_ERROR( this->getNumRows(), diagonal.size() )

memory::ScalarType arrayType = diagonal.getValueType();

switch ( arrayType )
{
#define LAMA_SCALE_CALL( z, I, _ )                                                  \
case common::scalar::SCALAR_ARITHMETIC_TYPE##I:                                     \
{                                                                           \
    const LAMAArray<ARITHMETIC_TYPE##I>& typedDiagonal =                    \
            dynamic_cast<const LAMAArray<ARITHMETIC_TYPE##I>&>( diagonal );     \
    static_cast<Derived*>( this )->scaleImpl( typedDiagonal );              \
    break;                                                                  \
}                                                                           \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_SCALE_CALL, _ )

#undef LAMA_SCALE_CALL

default:
    COMMON_THROWEXCEPTION( "scale of type " << arrayType << " not supported" )
}
}

/** Implementation of MatrixStorage::getTypeName for derived class. */

virtual const char* getTypeName() const
{
// each derived class provides static method to get the type name.

return Derived::typeName();
}
};

}
// namespace lama

#endif // LAMA_CRTP_MATRIX_STORAGE_HPP_
