/**
 * @file DIAStorage.cpp
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
 * @brief Instantiation for template class DIAStorage.
 * @author Thomas Brandes
 * @date 04.06.2011
 * $Id$
 */

// hpp
#include <lama/storage/DIAStorage.hpp>

// others
#include <lama/openmp/OpenMPUtils.hpp>
#include <lama/openmp/OpenMPCSRUtils.hpp>
#include <lama/openmp/OpenMPDIAUtils.hpp>

#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/ContextAccess.hpp>

// macros
#include <lama/macros/unused.hpp>

// boost
#include <boost/scoped_array.hpp>

namespace lama
{

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DIAStorage<ValueType>::logger, "MatrixStorage.DIAStorage" );

/* --------------------------------------------------------------------------- */

template<>
const char* DIAStorage<float>::typeName()
{
    return "DIAStorage<float>";
}

template<>
const char* DIAStorage<double>::typeName()
{
    return "DIAStorage<double>";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( const IndexType numRows, const IndexType numColumns )

    : CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>( numRows, numColumns ), mNumDiagonals( 0 )
{
    LAMA_LOG_DEBUG( logger, "DIAStorage for matrix " << mNumRows << " x " << mNumColumns << ", no non-zero elements" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage()
    : CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>( 0, 0 ), mNumDiagonals( 0 )
{
    LAMA_LOG_DEBUG( logger, "DIAStorage, matrix is 0 x 0." );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const LAMAArray<IndexType>& offsets,
    const LAMAArray<ValueType>& values )

    : CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>( numRows, numColumns ), mNumDiagonals(
        numDiagonals ), mOffset( offsets ), mValues( values )
{
    // set diagonal property inherited as given

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( const DIAStorage<ValueType>& other )

    : CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>( 0, 0 )
{
    // @todo: copy of same storage format should be implemented more efficiently

    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::print() const
{
    using std::cout;
    using std::endl;

    cout << "DIAStorage " << mNumRows << " x " << mNumColumns << ", #diags = " << mNumDiagonals << endl;

    HostReadAccess<IndexType> offset( mOffset );
    HostReadAccess<ValueType> values( mValues );

    cout << "Diagonal offsets:";
    for ( IndexType d = 0; d < mNumDiagonals; d++ )
    {
        cout << " " << offset[d];
    }
    cout << endl;

    for ( IndexType i = 0; i < mNumRows; i++ )
    {
        cout << "Row " << i << " :";

        for ( IndexType ii = 0; ii < mNumDiagonals; ++ii )
        {
            const IndexType j = i + offset[ii];

            if ( j < 0 )
            {
                continue;
            }

            if ( j >= mNumColumns )
            {
                break;
            }

            cout << " " << j << ":" << values[i + ii * mNumRows];
        }
        cout << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::clear()
{
    mNumRows = 0;
    mNumColumns = 0;
    mNumDiagonals = 0;

    mOffset.clear();
    mValues.clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorageFormat DIAStorage<ValueType>::getFormat() const
{
    return DIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setDiagonalImpl( const Scalar scalar )
{
    LAMA_ASSERT_ERROR( mDiagonalProperty, *this << ": has not diagonal property, cannot set diagonal" );

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    {
        HostWriteAccess<ValueType> wValues( mValues );
        HostReadAccess<IndexType> rOffset( mOffset );

        ValueType value = scalar.getValue<ValueType>();

        OpenMPUtils::setVal( wValues.get(), numDiagonalElements, value );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const
{
    LAMA_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" );

    HostWriteOnlyAccess<OtherType> wRow( row, mNumColumns );

    const HostReadAccess<IndexType> offset( mOffset );
    const HostReadAccess<ValueType> values( mValues );

    for ( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = 0.0;
    }

    for ( IndexType d = 0; d < mNumDiagonals; ++d )
    {
        IndexType j = i + offset[d];

        if ( j < 0 || j >= mNumColumns )
        {
            continue;
        }

        wRow[j] = static_cast<OtherType>( values[index( i, d )] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherType>& diagonal ) const
{
    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    HostWriteOnlyAccess<OtherType> wDiagonal( diagonal, numDiagonalElements );
    HostReadAccess<ValueType> rValues( mValues );

    // Diagonal is first column

    OpenMPUtils::set( wDiagonal.get(), rValues.get(), numDiagonalElements );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherType>& diagonal )
{
    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );
    numDiagonalElements = std::min( numDiagonalElements, diagonal.size() );

    {
        HostReadAccess<OtherType> rDiagonal( diagonal );
        HostWriteAccess<ValueType> wValues( mValues );

        OpenMPUtils::set( wValues.get(), rDiagonal.get(), numDiagonalElements );
    }

    if ( LAMA_LOG_TRACE_ON( logger ) )
    {
        LAMA_LOG_TRACE( logger, "DIA after setDiagonal diagonal" );
        print();
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    HostWriteAccess<ValueType> wValues( mValues );

    ValueType value = scalar.getValue<ValueType>();

    for ( IndexType i = 0; i < mValues.size(); i++ )
    {
        wValues[i] *= value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::scaleImpl( const LAMAArray<OtherType>& diagonal )
{
    {
        HostReadAccess<OtherType> rDiagonal( diagonal );
        HostWriteAccess<ValueType> wValues( mValues );
        HostReadAccess<IndexType> rOffset( mOffset );

        for ( IndexType i = 0; i < mNumRows; i++ )
        {
            for ( IndexType ii = 0; ii < mNumDiagonals; ++ii )
            {
                const IndexType j = i + rOffset[ii];
                if ( j >= mNumColumns )
                {
                    break;
                }
                if ( j >= 0 )
                {
                    wValues[ii * mNumRows + i] *= static_cast<ValueType>( rDiagonal[j] );
                }
            }
        }
    }

    if ( LAMA_LOG_TRACE_ON( logger ) )
    {
        LAMA_LOG_TRACE( logger, "DIA after scale diagonal" );
        print();
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool DIAStorage<ValueType>::checkDiagonalProperty() const
{
    if ( mNumRows != mNumColumns )
    {
        return false;
    }

    if ( mOffset.size() == 0 )
    {
        return false;
    }

    // diagonal property is given if first diagonal is the main one

    HostReadAccess<IndexType> offset( mOffset );
    return ( offset[0] == 0 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::check( const char* /* msg */) const
{
    LAMA_ASSERT_EQUAL_ERROR( mNumDiagonals, mOffset.size() );

    if ( mNumDiagonals == 0 )
    {
        LAMA_ASSERT_EQUAL_ERROR( false, mDiagonalProperty );
    }

    LAMA_ASSERT_EQUAL_ERROR( mNumDiagonals * mNumRows, mValues.size() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setIdentity( const IndexType size )
{
    LAMA_LOG_DEBUG( logger, "set identity, size = " << size );

    mNumRows = size;
    mNumColumns = size;

    mNumDiagonals = 1; // identity has exactly one diagonal

    HostWriteOnlyAccess<IndexType> offset( mOffset, mNumDiagonals );
    HostWriteOnlyAccess<ValueType> values( mValues, mNumRows );

    offset[0] = 0;

    ValueType one = static_cast<ValueType>( 1.0 );

    OpenMPUtils::setVal( values.get(), mNumRows, one );

    mDiagonalProperty = true;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setUsedDiagonal(
    bool upperDiagonalUsed[],
    bool lowerDiagonalUsed[],
    IndexType i,
    IndexType j )
{
    if ( j >= i )
    {
        bool& flag = upperDiagonalUsed[j - i];

        // set flag only if not already set, improves cache usage

        if ( !flag )
        {
            flag = true; // write only if not true,
        }
    }
    else
    {
        bool& flag = lowerDiagonalUsed[i - j];

        if ( !flag )
        {
            flag = true; // write only if not true,
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::buildCSR(
    LAMAArray<IndexType>& ia,
    LAMAArray<IndexType>* ja,
    LAMAArray<OtherValueType>* values,
    const ContextPtr /* loc */) const
{
    // TODO all done on host, so loc is unused

    LAMA_LOG_INFO( logger,
                   "buildTypedCSRData<" << Scalar::getType<OtherValueType>() << ">" << " from DIA<" << Scalar::getType<ValueType>() << "> = " << *this << ", diagonal property = " << mDiagonalProperty );

    HostReadAccess<IndexType> diaOffsets( mOffset );
    HostReadAccess<ValueType> diaValues( mValues );

    ValueType eps = 0.0;

    HostWriteOnlyAccess<IndexType> csrIA( ia, mNumRows + 1 );

    OpenMPDIAUtils::getCSRSizes( csrIA.get(), mDiagonalProperty, mNumRows, mNumColumns, mNumDiagonals, diaOffsets.get(),
                                 diaValues.get(), eps );

    if ( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = OpenMPCSRUtils::sizes2offsets( csrIA.get(), mNumRows );

    LAMA_LOG_INFO( logger, "CSR: #non-zero values = " << numValues );

    HostWriteOnlyAccess<IndexType> csrJA( *ja, numValues );
    HostWriteOnlyAccess<OtherValueType> csrValues( *values, numValues );

    OpenMPDIAUtils::getCSRValues( csrJA.get(), csrValues.get(), csrIA.get(), mDiagonalProperty, mNumRows, mNumColumns,
                                  mNumDiagonals, diaOffsets.get(), diaValues.get(), eps );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType UNUSED( numValues ),
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<OtherValueType>& values,
    ContextPtr UNUSED( loc ) )
{
    // loc is ignored, we do it on the Host

    HostReadAccess<IndexType> csrIA( ia );
    HostReadAccess<IndexType> csrJA( ja );
    HostReadAccess<OtherValueType> csrValues( values );

    _MatrixStorage::init( numRows, numColumns );

    LAMA_LOG_DEBUG( logger, "fill DIA sparse maxtrix " << mNumRows << " x " << mNumColumns << " from csr data" );

    // build a set of all used lower and upper diagonals

    IndexType maxNumDiagonals = std::max( mNumRows, mNumColumns );

    boost::scoped_array<bool> upperDiagonalUsed( new bool[maxNumDiagonals] );
    boost::scoped_array<bool> lowerDiagonalUsed( new bool[maxNumDiagonals] );

    for ( IndexType i = 0; i < maxNumDiagonals; i++ )
    {
        upperDiagonalUsed[i] = false;
        lowerDiagonalUsed[i] = false;
    }

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; jj++ )
        {
            IndexType j = csrJA[jj]; // column used
            setUsedDiagonal( upperDiagonalUsed.get(), lowerDiagonalUsed.get(), i, j );
        }
    }

    mDiagonalProperty = OpenMPCSRUtils::hasDiagonalProperty( numRows, csrIA.get(), csrJA.get() );

    // mDiagonalProperty forces upper diagonal to be the first one

    setOffsets( maxNumDiagonals, upperDiagonalUsed.get(), lowerDiagonalUsed.get() );

    // now we can allocate and set the values
    {
        HostReadAccess<IndexType> offset( mOffset );
        HostWriteOnlyAccess<ValueType> myValues( mValues, mNumDiagonals * mNumRows );

        #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
        for ( int i = 0; i < mNumRows; i++ )
        {
            for ( IndexType d = 0; d < mNumDiagonals; d++ )
            {
                IndexType j = i + offset[d];

                if ( j < 0 || j >= mNumColumns )
                {
                    continue;
                }

                ValueType& addrValue = myValues[index( i, d )];

                // check for j >= 0 and j < mNumColumns not needed here

                addrValue = 0.0;

                for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    if ( csrJA[jj] == j )
                    {
                        addrValue = static_cast<ValueType>( csrValues[jj] );
                        break;
                    }
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setOffsets(
    const IndexType maxNumDiagonals,
    const bool upperDiagonalUsed[],
    const bool lowerDiagonalUsed[] )
{
    // Help routine to set offsets from used diagonals

    mNumDiagonals = 0;

    IndexType firstIndex = 0;

    if ( mDiagonalProperty )
    {
        firstIndex = 1;
        mNumDiagonals = 1;
    }

    for ( IndexType i = firstIndex; i < maxNumDiagonals; i++ )
    {
        if ( upperDiagonalUsed[i] )
        {
            mNumDiagonals++;
        }
        if ( lowerDiagonalUsed[i] )
        {
            mNumDiagonals++;
        }
    }

    LAMA_LOG_INFO( logger, "storage data requires " << mNumDiagonals << " diagonals a " << mNumRows << " values" );

    HostWriteOnlyAccess<IndexType> wOffset( mOffset, mNumDiagonals );

    if ( mNumDiagonals > 0 )
    {
        mNumDiagonals = 0;
        firstIndex = 0;

        if ( mDiagonalProperty )
        {
            wOffset[mNumDiagonals++] = 0;
            firstIndex = 1;
        }

        for ( IndexType i = maxNumDiagonals - 1; i > 0; i-- )
        {
            if ( lowerDiagonalUsed[i] )
            {
                wOffset[mNumDiagonals++] = -i;
            }
        }

        LAMA_LOG_INFO( logger, "lower diagonals = " << mNumDiagonals );

        for ( IndexType i = firstIndex; i < maxNumDiagonals; i++ )
        {
            if ( upperDiagonalUsed[i] )
            {
                wOffset[mNumDiagonals++] = i;
            }
        }
    }

    LAMA_ASSERT_EQUAL_DEBUG( mNumDiagonals, wOffset.size() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::~DIAStorage()
{
    LAMA_LOG_DEBUG( logger,
                    "~DIAStorage for matrix " << mNumRows << " x " << mNumColumns << ", # diags = " << mNumDiagonals );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::purge()
{
    mNumColumns = 0;
    mNumRows = 0;
    mNumDiagonals = 0;

    mOffset.purge();
    mValues.purge();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    LAMA_LOG_INFO( logger, "allocate DIA sparse matrix of size " << numRows << " x " << numColumns );

    mNumRows = numRows;
    mNumColumns = numColumns;

    // clear arrays to avoid unnecessary data transfer (write-only)

    mValues.clear();
    mOffset.clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "DIA(rows=" << mNumRows << ",cols=" << mNumColumns << ",nd=" << mNumDiagonals << ")";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::maxNorm() const
{
    LAMA_LOG_INFO( logger, *this << ": maxNorm()" );

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( absMaxVal, loc, DIAUtils, Reductions, ValueType );

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc );

    ValueType maxval = absMaxVal( mNumRows, mNumColumns, mNumDiagonals, diaOffsets.get(), diaValues.get() );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::getValue( IndexType i, IndexType j ) const
{
    LAMA_LOG_DEBUG( logger, "get value (" << i << ", " << j << ") from " << *this );

    const HostReadAccess<IndexType> offset( mOffset );

    ValueType myValue = 0.0;

    // check for a matching diagonal element in the row i

    for ( IndexType d = 0; d < mNumDiagonals; ++d )
    {
        if ( i + offset[d] == j )
        {
            const HostReadAccess<ValueType> values( mValues );
            LAMA_LOG_DEBUG( logger,
                            "get value (" << i << ", " << j << ") is diag = " << d << ", offset = " << offset[d] << ", index = " << index(i,d) );
            myValue = values[index( i, d )];
            break;
        }
    }

    return myValue;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mOffset.prefetch( location );
    mValues.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DIAStorage<ValueType>::getNumDiagonals() const
{
    return mNumDiagonals;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<IndexType>& DIAStorage<ValueType>::getOffsets() const
{
    return mOffset;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<ValueType>& DIAStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::wait() const
{
    mOffset.wait();
    mValues.wait();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::swap( DIAStorage<ValueType>& other )
{
    std::swap( mNumDiagonals, other.mNumDiagonals );
    mOffset.swap( other.mOffset );
    mValues.swap( other.mValues );

    MatrixStorage<ValueType>::swap( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t DIAStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof(IndexType);
    memoryUsage += sizeof(IndexType) * mOffset.size();
    memoryUsage += sizeof(ValueType) * mValues.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::matrixTimesVector(

    LAMAArrayView<ValueType> result,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const

{
    LAMA_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y );

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumColumns );
    LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumRows );

    // @todo support on GPU for matrixTimesVector with DIA format

    ContextPtr loc = ContextFactory::getContext( Context::Host );

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc );

    LAMA_INTERFACE_FN_T( normalGEMV, loc, DIAUtils, Mult, ValueType );

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        LAMA_LOG_INFO( logger, "result == y" );
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        LAMA_CONTEXT_ACCESS( loc );

        normalGEMV( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets.get(), diaValues.get(), NULL );
    }
    else
    {
        LAMA_LOG_INFO( logger, "result != y" );
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        LAMA_CONTEXT_ACCESS( loc );

        normalGEMV( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets.get(), diaValues.get(), NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::jacobiIterate(
    LAMAArrayView<ValueType> solution,
    const LAMAArrayConstView<ValueType> oldSolution,
    const LAMAArrayConstView<ValueType> rhs,
    const ValueType omega ) const
{
    LAMA_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." );

    LAMA_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" );

    if ( solution == oldSolution )
    {
        LAMA_THROWEXCEPTION( "alias of solution and oldSolution unsupported" );
    }

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() );
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() );
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns );
    // matrix must be square

    HostWriteAccess<ValueType> wSolution( solution );
    HostReadAccess<IndexType> diaOffset( mOffset );
    HostReadAccess<ValueType> diaValues( mValues );
    HostReadAccess<ValueType> rOldSolution( oldSolution );
    HostReadAccess<ValueType> rRhs( rhs );

    OpenMPDIAUtils::jacobi( wSolution.get(), mNumColumns, mNumDiagonals, diaOffset.get(), diaValues.get(),
                            rOldSolution.get(), rRhs.get(), omega, mNumRows, NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::create() const
{
    return new DIAStorage<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::copy() const
{
    return new DIAStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template class LAMA_DLL_IMPORTEXPORT DIAStorage<float> ;
template class LAMA_DLL_IMPORTEXPORT DIAStorage<double> ;

}
