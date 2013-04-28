/**
 * @file SparseAssemblyStorage.cpp
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
 * @brief SparseAssemblyStorage.cpp
 * @author Jiri Kraus, Thomas Brandes
 * @date 07.11.2011
 * $Id$
 */

// hpp
#include <lama/storage/SparseAssemblyStorage.hpp>

#include <lama/openmp/OpenMPUtils.hpp>
#include <lama/openmp/OpenMPCSRUtils.hpp>

#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>

#include <cmath>

namespace lama
{

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, SparseAssemblyStorage<ValueType>::logger,
                              "MatrixStorage.SparseAssemblyStorage" )

/* --------------------------------------------------------------------------- */

template<typename T>
SparseAssemblyStorage<T>::SparseAssemblyStorage()
    : CRTPMatrixStorage<SparseAssemblyStorage<ValueType>, ValueType> (0, 0),
      mRows(0), mNumValues(0)
{
}

/* --------------------------------------------------------------------------- */

template<typename T>
SparseAssemblyStorage<T>::SparseAssemblyStorage( const SparseAssemblyStorage<ValueType>& other )
    : CRTPMatrixStorage<SparseAssemblyStorage<ValueType>,ValueType>( other.getNumRows(),
            other.getNumColumns() ), mRows(
                other.mRows ), mNumValues( other.mNumValues )
{
}

/* --------------------------------------------------------------------------- */

template<typename T>
SparseAssemblyStorage<T>::SparseAssemblyStorage( const _MatrixStorage& other )
    : CRTPMatrixStorage<SparseAssemblyStorage<ValueType>,ValueType>( other.getNumRows(),
            other.getNumColumns() )
{
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename T>
SparseAssemblyStorage<T>::SparseAssemblyStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow /* = 10*/)
    : CRTPMatrixStorage<SparseAssemblyStorage<ValueType>,ValueType>( numRows, numColumns ), mRows(
        numRows ), mNumValues( 0 )
{
    LAMA_LOG_INFO( logger,
                   "Creating with " << mNumRows <<" rows, " << mNumColumns << " columns, " << numValuesPerRow << " values per row." )
    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        LAMA_LOG_TRACE( logger, "Reserving storage for row " << i )
        mRows[i].reserve( numValuesPerRow );
    }
    LAMA_LOG_DEBUG( logger, "Created." )
}

template<typename T>
SparseAssemblyStorage<T>::~SparseAssemblyStorage()
{
}

/* --------------------------------------------------------------------------- */

template<typename T>
SparseAssemblyStorage<T>& SparseAssemblyStorage<T>::operator=( const SparseAssemblyStorage<T>& other )
{
    mNumRows = other.mNumRows;
    mNumColumns = other.mNumColumns;
    mRows = other.mRows;
    mNumValues = other.mNumValues;
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::swap( SparseAssemblyStorage<ValueType>& other )
{
    std::swap( mNumValues, other.mNumValues );
    mRows.swap( other.mRows );

    MatrixStorage<ValueType>::swap( other );
}

/* --------------------------------------------------------------------------- */

template<typename T>
SparseAssemblyStorage<T>& SparseAssemblyStorage<T>::operator=( const _MatrixStorage& other )
{
    assign( other );
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::allocate( const IndexType numRows, const IndexType numColumns )
{
    LAMA_LOG_INFO( logger, "allocate sparse assembly storage " << numRows << " x " << numColumns )

    _MatrixStorage::init( numRows, numColumns );

    mNumValues = 0;

    mRows.clear(), mRows.resize( numRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorageFormat SparseAssemblyStorage<ValueType>::getFormat() const
{
    return ASSEMBLY;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::prefetch( const ContextPtr /* location */) const
{
    // not supported, no implementation required
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::wait() const
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::clear()
{
    // clear is here like a purge due to destructor of mRows

    mNumRows = 0;
    mNumColumns = 0;
    mNumValues = 0;

    mRows.clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::purge()
{
    mNumRows = 0;
    mNumColumns = 0;
    mNumValues = 0;

    mRows.clear();
}

/* --------------------------------------------------------------------------- */

#ifdef LAMA_ASSERT_LEVEL_OFF
template<typename ValueType>
void SparseAssemblyStorage<ValueType>::check( const char* ) const
{}
#else
template<typename ValueType>
void SparseAssemblyStorage<ValueType>::check( const char* msg ) const
{
    if ( mNumRows != static_cast<IndexType>( mRows.size() ) )
    {
        LAMA_THROWEXCEPTION(
            msg << ": SparseAssemblyStorage: mNumRows = " << mNumRows << " does not match size of mRows = " << mRows.size() );
    }

    // TODO check that numValues is equal sum( mRows[i].ja.size() ), 0 <= i < mNumRows
    // TODO check for valid column indexes
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType SparseAssemblyStorage<ValueType>::maxNorm() const
{
    ValueType maxval = static_cast<ValueType>( 0.0 );

    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        const std::vector<ValueType>& values = mRows[i].values;

        for ( size_t jj = 0; jj < values.size(); ++jj )
        {
            const ValueType val = std::abs( mRows[i].values[jj] );

            if ( val > maxval )
            {
                maxval = val;
            }
        }
    }

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t SparseAssemblyStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof(IndexType) * mNumValues;
    memoryUsage += sizeof(ValueType) * mNumValues;
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool SparseAssemblyStorage<ValueType>::checkDiagonalProperty() const
{
    bool diagonalProperty = true;

    const IndexType numDiagonals = std::min( mNumRows, mNumColumns );

    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        const Row& row = mRows[i];
        if ( row.ja.size() == 0 )
        {
            diagonalProperty = false;
            break;
        }
        if ( row.ja[0] != i )
        {
            diagonalProperty = false;
            break;
        }
    }
    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename T>
SparseAssemblyStorage<T>::Row::Row()
{
}

template<typename T>
SparseAssemblyStorage<T>::Row::Row( const IndexType numValuesPerRow )
{
    ja.reserve( numValuesPerRow );
    values.reserve( numValuesPerRow );
}

template<typename T>
SparseAssemblyStorage<T>::Row::Row( const typename SparseAssemblyStorage<T>::Row& other )
    : ja( other.ja ), values( other.values )
{
}

template<typename T>
typename SparseAssemblyStorage<T>::Row& SparseAssemblyStorage<T>::Row::operator=(
    const typename SparseAssemblyStorage<T>::Row& other )
{
    ja = other.ja;
    values = other.values;
    return *this;
}

template<typename T>
void SparseAssemblyStorage<T>::Row::reserve( const IndexType numValuesPerRow )
{
    ja.reserve( numValuesPerRow );
    values.reserve( numValuesPerRow );
}

template<typename T>
void SparseAssemblyStorage<T>::Row::scale( const ValueType val )
{
    for ( size_t i = 0; i < values.size(); i++ )
    {
        values[i] *= val;
    }
}

template<typename T>
typename SparseAssemblyStorage<T>::ValueType SparseAssemblyStorage<T>::operator()(
    const IndexType i,
    const IndexType j ) const
{
    if ( j >= mNumColumns )
    {
        LAMA_THROWEXCEPTION( "Passed column Index " << j << " exceeds column count " << mNumColumns << "." )
    }

    const std::vector<IndexType>& rJA = mRows[i].ja;

    for ( size_t k = 0; k < mRows[i].ja.size(); ++k )
    {
        if ( j == rJA[k] )
        {
            const std::vector<ValueType>& rValues = mRows[i].values;
            return rValues[k];
        }
    }
    return 0.0;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::print() const
{
    using std::cout;
    using std::endl;

    cout << "AssemblyStorage " << mNumRows << " x " << mNumColumns << ", #values = " << mNumValues << endl;

    for ( IndexType i = 0; i < mNumRows; i++ )
    {
        const Row& row = mRows[i];

        cout << "Row " << i << " ( " << row.ja.size() << " ) :";

        for ( size_t k = 0; k < row.ja.size(); ++k )
        {
            cout << " " << row.ja[k] << ":" << row.values[k];
        }
        cout << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename T>
void SparseAssemblyStorage<T>::set( const IndexType i, const IndexType j, const ValueType value )
{
    if ( i >= mNumRows )
    {
        LAMA_THROWEXCEPTION( "Passed row Index " << i << " exceeds row count " << mNumRows << "." )
    }

    if ( j >= mNumColumns )
    {
        LAMA_THROWEXCEPTION( "Passed column Index " << j << " exceeds column count " << mNumColumns << "." )
    }

    {
        const std::vector<IndexType>& rJA = mRows[i].ja;

        for ( size_t k = 0; k < mRows[i].ja.size(); ++k )
        {
            if ( j == rJA[k] )
            {
                std::vector<ValueType>& wValues = mRows[i].values;

                LAMA_LOG_TRACE( logger,
                                "set( " << i << ", " << j << ", " << value << ") : override existing value " << wValues[k] )

                wValues[k] = value;
                return;
            }
        }
    }

    LAMA_LOG_TRACE( logger, "set( " << i << ", " << j << ", " << value << ") : new entry " )

    std::vector<IndexType>& wJA = mRows[i].ja;
    std::vector<ValueType>& wValues = mRows[i].values;

    wValues.push_back( value );
    wJA.push_back( j );
    ++mNumValues;

    // if we have a diagonal element and the row was not empty before we need to
    // fix the diagonal property

    if ( i == j && ( wJA.size() - 1 > 0 ) )
    {
        LAMA_LOG_TRACE( logger, "diagonal element swapped to first element of row" )

        std::swap( wValues[0], wValues[wValues.size() - 1] );
        std::swap( wJA[0], wJA[wJA.size() - 1] );
    }
}

template<typename T>
const std::vector<IndexType>&
SparseAssemblyStorage<T>::getJa( const IndexType i ) const
{
    return mRows[i].ja;
}

template<typename T>
std::vector<IndexType>& SparseAssemblyStorage<T>::getJa( const IndexType i )
{
    return mRows[i].ja;
}

template<typename T>
const std::vector<typename SparseAssemblyStorage<T>::ValueType>&
SparseAssemblyStorage<T>::getValues( const IndexType i ) const
{
    return mRows[i].values;
}

template<typename T>
IndexType SparseAssemblyStorage<T>::getNumValues() const
{
    return mNumValues;
}

template<typename T>
void SparseAssemblyStorage<T>::setRow(
    const IndexType i,
    const LAMAArrayConstView<IndexType>& ja,
    const LAMAArrayConstView<ValueType>& values )
{
    //LAMA_ASSERT_EQUAL_ERROR( ja.size(), values.size() )

    #pragma omp atomic
    mNumValues -= mRows[i].ja.size();

    mRows[i].ja.resize( ja.size() );
    mRows[i].values.resize( values.size() );

    HostReadAccess<IndexType> rJA( ja );
    HostReadAccess<ValueType> rValues( values );

    for ( IndexType k = 0; k < ja.size(); ++k )
    {
        mRows[i].ja[k] = rJA[k];
        mRows[i].values[k] = rValues[k];
    }
    const IndexType rowSize = mRows[i].ja.size();
    #pragma omp atomic
    mNumValues += rowSize;
}

/* --------------------------------------------------------------------------- */

template<typename T>
void SparseAssemblyStorage<T>::fixDiagonalProperty( const IndexType i )

{
    // fix diagonal property if necessary

    if ( i >= mNumColumns )
    {
        return;
    }

    if ( mRows[i].ja.size() == 0 )
    {
        #pragma omp atomic
        ++mNumValues;
        mRows[i].ja.push_back( i );
        mRows[i].values.push_back( 0.0 );
        return;
    }

    if ( mRows[i].ja[0] == i )
    {
        return;
    }

    // try to find diagonal element

    std::vector<IndexType>& wJA = mRows[i].ja;
    std::vector<ValueType>& wValues = mRows[i].values;

    for ( size_t k = 0; k < wJA.size(); ++k )
    {
        if ( i == wJA[k] )
        {
            std::swap( wValues[0], wValues[k] );
            std::swap( wJA[0], wJA[k] );
            break;
        }
    }

    if ( wJA[0] != i )
    {
        #pragma omp atomic
        ++mNumValues;
        wJA.push_back( i );
        wValues.push_back( 0.0 );
        std::swap( wValues[0], wValues[wValues.size() - 1] );
        std::swap( wJA[0], wJA[wValues.size() - 1] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename T>
void SparseAssemblyStorage<T>::setNumColumns( const IndexType numColumns )
{
    mNumColumns = numColumns;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::setIdentity( const IndexType n )
{
    allocate( n, n );

    ValueType one = static_cast<ValueType>( 1.0 );

    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        set( i, i, one );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void SparseAssemblyStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<OtherValueType>& values,
    const ContextPtr /* loc */)
{
    // no more error checks here on the sizes, but on the content

    HostReadAccess<IndexType> csrIA( ia );
    HostReadAccess<IndexType> csrJA( ja );
    HostReadAccess<OtherValueType> csrValues( values );

    if ( !OpenMPCSRUtils::validOffsets( csrIA.get(), numRows, numValues ) )
    {
        LAMA_THROWEXCEPTION( "invalid offset array" )
    }

    if ( !OpenMPUtils::validIndexes( csrJA.get(), numValues, numColumns ) )
    {
        LAMA_THROWEXCEPTION( "invalid column indexes in ja = " << ja << ", #columns = " << numColumns )
    }

    mNumRows = numRows;
    mNumColumns = numColumns;
    mNumValues = numValues;

    mRows.resize( mNumRows );

    LAMA_ASSERT_EQUAL_ERROR( csrIA[numRows], numValues )

    LAMA_LOG_DEBUG( logger, "fill " << *this << " with csr data, " << numValues << " non-zero values" )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        const IndexType n = csrIA[i + 1] - csrIA[i];

        Row& row = mRows[i];

        row.ja.resize( n );
        row.values.resize( n );

        IndexType offset = 0;

        for ( IndexType k = csrIA[i]; k < csrIA[i + 1]; ++k )
        {
            row.ja[offset] = csrJA[k];
            row.values[offset] = static_cast<ValueType>( csrValues[k] );
            ++offset;
        }
    }

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void SparseAssemblyStorage<ValueType>::buildCSR(
    LAMAArray<IndexType>& ia,
    LAMAArray<IndexType>* ja,
    LAMAArray<OtherValueType>* values,
    const ContextPtr /* loc */) const
{
    // TODO all done on host, so loc is unused

    LAMA_LOG_INFO( logger, *this << ": build CSR data from it" )

    HostWriteOnlyAccess<IndexType> csrIA( ia, mNumRows + 1 );

    // build csrSizes in ia

    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        csrIA[i] = mRows[i].ja.size();
        LAMA_ASSERT_EQUAL_DEBUG( mRows[i].ja.size(), mRows[i].values.size() )
    }

    if ( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    // build csrOffset in ia from the sizes

    OpenMPCSRUtils::sizes2offsets( csrIA.get(), mNumRows );

    // copy ja, values

    HostWriteOnlyAccess<IndexType> csrJA( *ja, mNumValues );
    HostWriteOnlyAccess<OtherValueType> csrValues( *values, mNumValues );

    #pragma omp parallel for
    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        IndexType offset = 0;

        for ( IndexType k = csrIA[i]; k < csrIA[i + 1]; ++k )
        {
            csrJA[k] = mRows[i].ja[offset];
            csrValues[k] = static_cast<OtherValueType>( mRows[i].values[offset] );
            ++offset;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void SparseAssemblyStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherValueType>& diagonal )
{
    IndexType numDiagonalElements = diagonal.size();

    HostReadAccess<OtherValueType> rDiagonal( diagonal );

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        set( i, i, static_cast<ValueType>( rDiagonal[i] ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void SparseAssemblyStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherValueType>& diagonal ) const
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    HostWriteOnlyAccess<OtherValueType> wDiagonal( diagonal, numDiagonalElements );

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        wDiagonal[i] = static_cast<OtherValueType>( operator()( i, i ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void SparseAssemblyStorage<ValueType>::getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const
{
    LAMA_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    HostWriteOnlyAccess<OtherType> wRow( row, mNumColumns );

    for ( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = 0.0;
    }

    const std::vector<IndexType>& ja = mRows[i].ja;
    const std::vector<ValueType>& values = mRows[i].values;

    LAMA_ASSERT_EQUAL_DEBUG( ja.size(), values.size() )

    for ( size_t k = 0; k < ja.size(); ++k )
    {
        const IndexType j = ja[k];
        LAMA_ASSERT_DEBUG( j >= 0 && j < mNumColumns,
                           "col index " << j << " out of range" << " is at row " << i << ":" << k )

        wRow[j] = static_cast<OtherType>( values[k] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::setDiagonalImpl( const Scalar value )
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    ValueType val = value.getValue<ValueType>();

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        set( i, i, val );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    ValueType value = scalar.getValue<ValueType>();

    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        mRows[i].scale( value );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void SparseAssemblyStorage<ValueType>::scaleImpl( const LAMAArray<OtherValueType>& diagonal )
{
    IndexType n = std::min( mNumRows, diagonal.size() );

    HostReadAccess<OtherValueType> rDiagonal( diagonal );

    for ( IndexType i = 0; i < n; ++i )
    {
        mRows[i].scale( static_cast<ValueType>( rDiagonal[i] ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void SparseAssemblyStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "SparseAssemblyStorage: (" << mNumRows << " x " << mNumColumns << ", #values = " << mNumValues
           << ", diag = " << mDiagonalProperty << " )";
}

/* --------------------------------------------------------------------------- */

template<>
const char* SparseAssemblyStorage<float>::typeName()
{
    return "SparseAssemblyStorage<float>";
}

template class LAMA_DLL_IMPORTEXPORT SparseAssemblyStorage<float> ;

template<>
const char* SparseAssemblyStorage<double>::typeName()
{
    return "SparseAssemblyStorage<double>";
}

template class LAMA_DLL_IMPORTEXPORT SparseAssemblyStorage<double> ;

} //namespace lama

