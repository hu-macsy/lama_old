/**
 * @file CRTPMatrixStorage.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

#ifndef LAMA_CRTP_MATRIX_STORAGE_HPP_
#define LAMA_CRTP_MATRIX_STORAGE_HPP_

// for dll_import
#include <lama/config.hpp>

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
class LAMA_DLL_IMPORTEXPORT CRTPMatrixStorage: public MatrixStorage<ValueType>
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
        Scalar::ScalarType arrayType = values.getValueType();

        if( arrayType == Scalar::DOUBLE )
        {
            const LAMAArray<double>& typedValues = dynamic_cast<const LAMAArray<double>&>( values );
            static_cast<Derived*>( this )->setCSRDataImpl( numRows, numColumns, numValues, ia, ja, typedValues,
                    this->getContextPtr() );
        }
        else if( arrayType == Scalar::FLOAT )
        {
            const LAMAArray<float>& typedValues = dynamic_cast<const LAMAArray<float>&>( values );
            static_cast<Derived*>( this )->setCSRDataImpl( numRows, numColumns, numValues, ia, ja, typedValues,
                    this->getContextPtr() );
        }
        else
        {
            LAMA_THROWEXCEPTION( *this << ": setCSRData with value type " << arrayType << " not supported" );
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

    void buildCSRData( LAMAArray<IndexType>& ia, LAMAArray<IndexType>& ja, _LAMAArray& values ) const
    {
        Scalar::ScalarType arrayType = values.getValueType();

        if( arrayType == Scalar::DOUBLE )
        {
            LAMAArray<double>& typedValues = dynamic_cast<LAMAArray<double>&>( values );
            static_cast<const Derived*>( this )->buildCSR( ia, &ja, &typedValues, this->getContextPtr() );
        }
        else if( arrayType == Scalar::FLOAT )
        {
            LAMAArray<float>& typedValues = dynamic_cast<LAMAArray<float>&>( values );
            static_cast<const Derived*>( this )->buildCSR( ia, &ja, &typedValues, this->getContextPtr() );
        }
        else
        {
            LAMA_THROWEXCEPTION( *this << ": build CSR with value type " << arrayType << " not supported" );
        }
    }

    /** Get the i-th row of a storage as LAMA array. */

    void getRow( _LAMAArray& row, const IndexType i ) const
    {
        Scalar::ScalarType arrayType = row.getValueType();

        if( arrayType == Scalar::DOUBLE )
        {
            LAMAArray<double>& typedRow = dynamic_cast<LAMAArray<double>&>( row );
            static_cast<const Derived*>( this )->getRowImpl( typedRow, i );
        }
        else if( arrayType == Scalar::FLOAT )
        {
            LAMAArray<float>& typedRow = dynamic_cast<LAMAArray<float>&>( row );
            static_cast<const Derived*>( this )->getRowImpl( typedRow, i );
        }
        else
        {
            LAMA_THROWEXCEPTION( "getRow for array of type " << arrayType << " not supported" );
        }
    }

    void getDiagonal( _LAMAArray& diagonal ) const
    {
        if( !this->hasDiagonalProperty() )
        {
            LAMA_THROWEXCEPTION( *this << ": has not diagonal property, cannot set diagonal" );
        }

        Scalar::ScalarType arrayType = diagonal.getValueType();

        if( arrayType == Scalar::DOUBLE )
        {
            LAMAArray<double>& typedDiagonal = dynamic_cast<LAMAArray<double>&>( diagonal );
            static_cast<const Derived*>( this )->getDiagonalImpl( typedDiagonal );
        }
        else if( arrayType == Scalar::FLOAT )
        {
            LAMAArray<float>& typedDiagonal = dynamic_cast<LAMAArray<float>&>( diagonal );
            static_cast<const Derived*>( this )->getDiagonalImpl( typedDiagonal );
        }
        else
        {
            LAMA_THROWEXCEPTION( "getDiagonal for array of type " << arrayType << " not supported" );
        }
    }

    void setDiagonal( const Scalar value )
    {
        static_cast<Derived*>( this )->setDiagonalImpl( value );
    }

    void setDiagonal( const _LAMAArray& diagonal )
    {
        IndexType numDiagonalElements = diagonal.size();

        if( numDiagonalElements > this->getNumRows() || numDiagonalElements > this->getNumColumns() )
        {
            LAMA_THROWEXCEPTION( "Diagonal of size " << numDiagonalElements << " too large for matrix: " << *this );
        }

        if( !this->hasDiagonalProperty() )
        {
            LAMA_THROWEXCEPTION( *this << ": has not diagonal property, cannot set diagonal" );
        }

        Scalar::ScalarType arrayType = diagonal.getValueType();

        if( arrayType == Scalar::DOUBLE )
        {
            const LAMAArray<double>& typedDiagonal = dynamic_cast<const LAMAArray<double>&>( diagonal );
            static_cast<Derived*>( this )->setDiagonalImpl( typedDiagonal );
        }
        else if( arrayType == Scalar::FLOAT )
        {
            const LAMAArray<float>& typedDiagonal = dynamic_cast<const LAMAArray<float>&>( diagonal );
            static_cast<Derived*>( this )->setDiagonalImpl( typedDiagonal );
        }
        else
        {
            LAMA_THROWEXCEPTION( "setDiagonal to array of type " << arrayType << " not supported" );
        }
    }

    void scale( const Scalar value )
    {
        static_cast<Derived*>( this )->scaleImpl( value );
    }

    /** Polymorph implementation for MatrixStorage<T>::scale */

    void scale( const _LAMAArray& diagonal )
    {
        LAMA_ASSERT_EQUAL_ERROR( this->getNumRows(), diagonal.size() );

        Scalar::ScalarType arrayType = diagonal.getValueType();

        if( arrayType == Scalar::DOUBLE )
        {
            const LAMAArray<double>& typedDiagonal = dynamic_cast<const LAMAArray<double>&>( diagonal );
            static_cast<Derived*>( this )->scaleImpl( typedDiagonal );
        }
        else if( arrayType == Scalar::FLOAT )
        {
            const LAMAArray<float>& typedDiagonal = dynamic_cast<const LAMAArray<float>&>( diagonal );
            static_cast<Derived*>( this )->scaleImpl( typedDiagonal );
        }
        else
        {
            LAMA_THROWEXCEPTION( "scale of type " << arrayType << " not supported" );
        }
    }

    /** Implementation of MatrixStorage::getTypeName for derived class. */

    virtual const char* getTypeName() const
    {
        // each derived class provides static method to get the type name.

        return Derived::typeName();
    }
};

} // namespace lama

#endif // LAMA_CRTP_MATRIX_STORAGE_HPP_
