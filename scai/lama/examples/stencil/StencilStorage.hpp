/**
 * @file StencilStorage.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Definition of a structure for a (non-distributed) sparse matrix based on a grid and a stencil on it
 * @author Thomas Brandes
 * @date 13.04.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/storage/CRTPMatrixStorage.hpp>
#include <scai/common/Grid.hpp>

#include "Stencil.hpp"

namespace scai
{

namespace lama
{

/** Storage format for a stencil matrix
 *
 *  The Stencil storage contains the following data:
 *
 *  - grid specifies the grid elements that will be connected
 *  - stencil specifies the pattern how elements are connected
 *
 *  Keep in mind that the elements of the grid are linearized.
 *
 *  @tparam ValueType is the value type of the matrix values.
 *
 *  Note: StencilStorage is not registed at factory at can only be
 *        created with a grid and a stencil
 *
 *  In contrary to the other matrix storage formats many operations are
 *  not supported:
 * 
 *  - copy constrcuctor with any other storage, e.g. StencilStorage s( csrStorage );
 *    but StencilStrage csr( stencilStorage ) is allowed
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT StencilStorage: public MatrixStorage<ValueType>
{
public:

    typedef ValueType StorageValueType;

    /** get typename of the matrix storage format. */

    static const char* typeName();

    /** Constructor of a stencil stroage
     *
     */
    StencilStorage( const common::Grid& grid, const Stencil<ValueType>&  stencil );

    /** Override default copy constructor to guarantee same behavior */

    StencilStorage( const StencilStorage<ValueType>& other )

        : MatrixStorage<ValueType>(),
          mGrid( other.mGrid ),
          mStencil( other.mStencil )
    {
        _MatrixStorage::setContextPtr( other.getContextPtr() );
    }
 
    /** 
     * @brief Implementation of pure method _MatrixStorage::allocate
     */
    void allocate( const IndexType /* numRows */, const IndexType /* numColumns */ )
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /**
     * @brief Implementation of pure method for _MatrixStorage.
     */
    virtual void clear()
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge()
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual StencilStorage* copy() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual StencilStorage* newMatrixStorage() const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Destructor of Stencil sparse storage. */

    virtual ~StencilStorage();

    /** Test the storage data for inconsistencies.
     *
     *  @param[in] msg is additional message string that should be used to identify calling routine
     *  @throws Exception in case of any inconsistency.
     */
    void check( const char* msg ) const
    {
        COMMON_THROWEXCEPTION( "check: " << msg );
    }

    /** Getter routine for the enum value that stands for this format. */

    virtual Format::MatrixStorageFormat getFormat() const
    {
        COMMON_THROWEXCEPTION( "getDiagonal unsuported" );
    }

    /** Implementation of pure method.  */

    virtual void setIdentity( const IndexType )
    {
        COMMON_THROWEXCEPTION( "getDiagonal unsuported" )
    }

    /** Implementation of pure method _MatrixStorage::checkDiagonalProperty */

    virtual bool checkDiagonalProperty() const
    {
        COMMON_THROWEXCEPTION( "getDiagonal unsuported" )
    }

    /** Implementation of pure method _MatrixStorage::swap */

    virtual void swap( _MatrixStorage& )
    {
        COMMON_THROWEXCEPTION( "getDiagonal unsuported" )
    }

    /** Implementation of pure method _MatrixStorage::getTypeName */

    virtual const char* getTypeName() const
    {
        COMMON_THROWEXCEPTION( "getTypeName unsuported" )
    }

    virtual scai::lama::MatrixStorageCreateKeyType getCreateValue() const;

    /** _MatrixStorage */

    virtual void getRow(scai::hmemo::_HArray&, IndexType) const
    {
        COMMON_THROWEXCEPTION( "getDiagonal unsuported" )
    }

    /** _MatrixStorage */

    virtual void setRow(const scai::hmemo::_HArray&, IndexType, scai::common::binary::BinaryOp)
    {
        COMMON_THROWEXCEPTION( "getDiagonal unsuported" )
    }

    /** _MatrixStorage */

 	virtual void setColumn(const scai::hmemo::_HArray&, IndexType, scai::common::binary::BinaryOp)
    {
        COMMON_THROWEXCEPTION( "getDiagonal unsuported" )
    }

    /** _MatrixStorage */

 	virtual void getDiagonal(scai::hmemo::_HArray&) const
    {
        COMMON_THROWEXCEPTION( "getDiagonal unsuported" )
    }

    /** _MatrixStorage */

 	virtual void setDiagonalV(const scai::hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "setDiagonal unsuported" )
    }

    /** _MatrixStorage */

 	virtual void scaleRows(const scai::hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "scaleRows unsuported" )
    }

    /** _MatrixStorage */

 	virtual void buildCSRSizes(scai::hmemo::HArray<int>&) const
    {
        COMMON_THROWEXCEPTION( "buildCSRSizes unsuported" )
    }

    /** _MatrixStorage */

 	virtual void buildCSRData(scai::hmemo::HArray<int>&, scai::hmemo::HArray<int>&, scai::hmemo::_HArray&) const
    {
        COMMON_THROWEXCEPTION( "buildCSRData unsuported" )
    }

    /** _MatrixStorage */

 	virtual void setCSRData(IndexType, IndexType, IndexType, const scai::hmemo::HArray<int>&, const scai::hmemo::HArray<int>&, const scai::hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "setcSRData unsuported" )
    }

    /** _MatrixStorage */

 	virtual void setDIAData(IndexType, IndexType, IndexType, const scai::hmemo::HArray<int>&, const scai::hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "setDIAData unsuported" )
    }

    /** MatrixStorage<ValueType> */

 	void scale(ValueType)
    {
        COMMON_THROWEXCEPTION( "scale unsuported" )
    }

    /** MatrixStorage<ValueType> */

 	void setDiagonal(ValueType)
    {
        COMMON_THROWEXCEPTION( "setDiagonal unsuported" )
    }

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "StencilStorage";
    }

    /** Getter routine for the number of stored values. */

    IndexType getNumValues() const
    {
        return 1;
    }

    /** Implementation of pure method MatrixStorage::getSparseRow */

    virtual void getSparseRow( hmemo::HArray<IndexType>& /* jA */, hmemo::_HArray& /* values */, const IndexType /* i */ ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Implementation of pure method MatrixStorage::getColumn */

    void getColumn( hmemo::_HArray& /* column */, const IndexType /* j */ ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Implementation of pure method MatrixStorage::getSparseColumn */

    virtual void getSparseColumn( hmemo::HArray<IndexType>&, hmemo::_HArray&, const IndexType ) const
    {
        COMMON_THROWEXCEPTION( "unsupported" )
    }

    /** Implementation of pure method.  */

    void conj();

    /** Get a value of the matrix.
     *
     * @param[in] i is the row index, 0 <= i < mNumRows
     * @param[in] j is the colum index, 0 <= j < mNumRows
     *
     * Out-of-range check is enabled for DEBUG version.
     */

    ValueType getValue( const IndexType, const IndexType ) const
    {
        COMMON_THROWEXCEPTION( "getValue unsupported" )
    }

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for Stencil storage */

    void setValue( const IndexType, const IndexType, const ValueType,
                   const common::binary::BinaryOp = common::binary::COPY )
    {
        COMMON_THROWEXCEPTION( "setValue unsupported" )
    }

    /** Initiate an asynchronous data transfer to a specified location. */

    virtual void prefetch( const hmemo::ContextPtr ) const
    {
        COMMON_THROWEXCEPTION( "print unsupported" )
    }

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const
    {
        COMMON_THROWEXCEPTION( "print unsupported" )
    }

    /** Implementation for pure method is provided. */

    virtual size_t getMemoryUsageImpl() const
    {
        COMMON_THROWEXCEPTION( "print unsupported" )
    }

    /** Override splitHalo with version that avoids unnecessary conversions.
     *
     *  This routine expects localData to be a StencilStorage and haloData to be
     *  a CSRStroage. Furthermore, colDist must be a grid distribution that fits 
     *  to the grid of this stencil storage.
     */

    virtual void splitHalo(
        MatrixStorage<ValueType>& /* localData */,
        MatrixStorage<ValueType>& /* haloData */,
        dmemo::Halo& /* halo */,
        const dmemo::Distribution& /* colDist */,
        const dmemo::Distribution* /* rowDist */ ) const
    {
        COMMON_THROWEXCEPTION( "print unsupported" )
    }
    

    /** General routine to build a CSR version of this stroage
     *
     *  @param[out] ia is the CSR offset array
     *  @param[out] ja is the array with the column indexes (optional)
     *  @param[out] values is the array with the non-zero matrix values (optional)
     *  @param[in]  loc is the Context where conversion should be done
     */

    template<typename OtherValueType>
    void buildCSR(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>* ja,
        hmemo::HArray<OtherValueType>* values,
        const hmemo::ContextPtr loc ) const;

    /** Implementation for MatrixStorage::l1Norm */

    virtual ValueType l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual ValueType l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual ValueType maxNorm() const;

    /** Implemenation of pure method of class MatrixStorage. */

    virtual void print( std::ostream& ) const
    {
        COMMON_THROWEXCEPTION( "print unsupported" )
    }
   

    using MatrixStorage<ValueType>::prefetch;
    using MatrixStorage<ValueType>::getContextPtr;

protected:

    using MatrixStorage<ValueType>::mNumRows;
    using MatrixStorage<ValueType>::mNumColumns;
    using MatrixStorage<ValueType>::mDiagonalProperty;

    common::Grid mGrid;    //! grid for which this matrix storage stands
    Stencil<ValueType> mStencil;      //! stencil that specifies connections between elements

private:

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
