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
#include <scai/lama/storage/MatrixStorage.hpp>

#include <scai/common/Grid.hpp>
#include <scai/common/Stencil.hpp>

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

    /* ==================================================================== */
    /*  static getter methods and corresponding pure methods                */
    /* ==================================================================== */

    /** Static method that returns a unique name for this storage class */

    static const char* typeName();

    /** Implementation of pure method _MatrixStorage:getTypeName    */

    virtual const char* getTypeName() const;

    /** Implementation of pure method _MatrixStorage:getCreateValue    */

    virtual MatrixStorageCreateKeyType getCreateValue() const;

    /* ==================================================================== */
    /*  Constructor / Destructor                                            */
    /* ==================================================================== */

    StencilStorage();

    /** Constructor of a stencil stroage
     *
     */
    StencilStorage( const common::Grid& grid, const common::Stencil<ValueType>&  stencil );

    /** Override default copy constructor to guarantee same behavior */

    StencilStorage( const StencilStorage<ValueType>& other ) : 

        MatrixStorage<ValueType>( other ),
        mGrid( other.mGrid ),
        mStencil( other.mStencil )

    {
    }

    /** 
     * @brief Implementation of pure method _MatrixStorage::allocate
     */
    void allocate( const IndexType numRows, const IndexType numColumns )
    {
        COMMON_THROWEXCEPTION( "unsupported: StencilStorage cannot be allocated by " << numRows << " x " << numColumns )
    }

    /**
     * @brief Implementation of pure method for _MatrixStorage.
     */
    virtual void clear();

    /** Implementation of pure method of class MatrixStorage. */

    virtual void purge();

    /** Implementation of MatrixStorage::copy for derived class. */

    virtual StencilStorage* copy() const;

    /** Implementation of MatrixStorage::newMatrixStorage for derived class. */

    virtual StencilStorage* newMatrixStorage( const IndexType, const IndexType ) const
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

    virtual Format getFormat() const
    {
        return Format::STENCIL;
    }

    /** Getter routine for the stencil used to construct this stencil storage. */

    const common::Stencil<ValueType>& getStencil() const
    {
        return mStencil;
    }

    /** Getter routine for the grid used to construct this stencil storage. */

    const common::Grid& getGrid() const
    {
        return mGrid;
    }

    /** Set identity uses trivial stencil for the given grid. */

    void setIdentity( const common::Grid& grid );

    /** Implementation of pure method MatrixStorage::setIdentity */

    virtual void setIdentity( const IndexType n )
    {
        setIdentity( common::Grid1D( n ) );
    }

    /******************************************************************/
    /*  set / get diagonal                                            */
    /******************************************************************/

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::getDiagonal
     */
    virtual void getDiagonal( hmemo::HArray<ValueType>& diagonal ) const;

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::setDiagonalV
     */
    virtual void setDiagonalV( const hmemo::HArray<ValueType>& )
    {
        COMMON_THROWEXCEPTION( "setDiagonal unsuported for stencil storage" )
    }

    /** 
     * Implementation of pure method MatrixStorage<ValueType>::setDiagonal
     */
    virtual void setDiagonal( const ValueType )
    {
        COMMON_THROWEXCEPTION( "setDiagonal unsuported for stencil storage" )
    }

    /** _MatrixStorage */

    virtual void scaleRows( const scai::hmemo::HArray<ValueType>& )
    {
        COMMON_THROWEXCEPTION( "scaleRows cannot be applied for stencil storage" )
    }

    virtual void scaleColumns( const scai::hmemo::HArray<ValueType>& )
    {
        COMMON_THROWEXCEPTION( "scaleColumns cannot be applied for stencil storage" )
    }

    /** Implementation of _MatrixStorage::buildCSRSizes */

 	virtual void buildCSRSizes( scai::hmemo::HArray<IndexType>& sizeIA ) const;

    /** Implementation of _MatrixStorage::buildCSRData */

 	virtual void buildCSRData(
       scai::hmemo::HArray<IndexType>& csrIA, 
       scai::hmemo::HArray<IndexType>& csrJA, 
       scai::hmemo::_HArray& csrValues ) const;

    /** Implementation of _MatrixStorage::setCSRData i
     *
     *  This method throws an exception as arbitrary CSR data cannot be converted to a stencil storage
     */
 	virtual void setCSRData( 
        IndexType, 
        IndexType, 
        const scai::hmemo::HArray<IndexType>&, 
        const scai::hmemo::HArray<IndexType>&, 
        const scai::hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "setcSRData unsuported" )
    }

    /** _MatrixStorage */

 	virtual void setDIAData( IndexType, IndexType, IndexType, const scai::hmemo::HArray<IndexType>&, const scai::hmemo::_HArray& )
    {
        COMMON_THROWEXCEPTION( "setDIAData cannot be applied for stencil storage" )
    }

    /** Implementation of pure method MatrixStorage<ValueType>::scale
     *
     *  This operation is supported as scale is available for a stencil.
     */
 	void scale( const ValueType factor );

    /* Print relevant information about matrix storage format. */

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "StencilStorage( " << this->getNumRows() << " x " << this->getNumColumns() 
               << ", grid = " << mGrid << ", stencil = " << mStencil << ")";
    }

    /** Getter routine for the number of stored values. */

    IndexType getNumValues() const
    {
        return 1;
    }

    /******************************************************************/
    /*  set - get  row - column                                       */
    /******************************************************************/

    /** Implementation of pure method MatrixStorage<ValueType>::getRow */

    virtual void getRow( hmemo::HArray<ValueType>& row, const IndexType i ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::getColumn */

    virtual void getColumn( hmemo::HArray<ValueType>& column, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::getSparseRow */

    virtual void getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const;

    /** Implementation of pure method MatrixStorage::getSparseColumn */

    virtual void getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setRow */

    virtual void setRow( const hmemo::HArray<ValueType>& row, const IndexType i, const common::BinaryOp op );

    /** Implementation of pure method MatrixStorage<ValueType>::setColumn */

    virtual void setColumn( const hmemo::HArray<ValueType>& column, const IndexType j, const common::BinaryOp op );

    /** Implementation of pure method.  */

    void conj();

    /** Implementation of pure method MatrixStorage::getValue */

    ValueType getValue( const IndexType i, const IndexType j ) const;

    /** Implementation of pure method MatrixStorage<ValueType>::setValue for Stencil storage */

    void setValue( const IndexType, const IndexType, const ValueType,
                   const common::BinaryOp = common::BinaryOp::COPY )
    {
        COMMON_THROWEXCEPTION( "setValue cannot be applied for stencil storage" )
    }

    /** Initiate an asynchronous data transfer to a specified location. */

    virtual void prefetch( const hmemo::ContextPtr ) const
    {
        // nothing to do here
    }

    /** Will wait for all outstanding asynchronous data transfers. */

    void wait() const
    {
    }

    /** Implementation for pure method is provided. */

    virtual size_t getMemoryUsageImpl() const;

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

    virtual RealType<ValueType> l1Norm() const;

    /** Implementation for MatrixStorage::l2Norm */

    virtual RealType<ValueType> l2Norm() const;

    /** Implementation for MatrixStorage::maxNorm */

    virtual RealType<ValueType> maxNorm() const;

    /** Implemenation of pure method of class MatrixStorage. */

    virtual void print( std::ostream& ) const
    {
        COMMON_THROWEXCEPTION( "print unsupported" )
    }
   
    /** Override the default implementation MatrixStorage<ValueType>::assign
     *
     *  Assign is only supported if other matrix is also a stencil storage.
     *  Otherwise an exception is thrown.
     */
    virtual void assign( const _MatrixStorage& other );

    /** Implementation of pure method MatrixStorage<ValueType>::assignDiagonal */

    virtual void assignDiagonal( const hmemo::HArray<ValueType>& )
    {
        COMMON_THROWEXCEPTION( "assignDiagonal unsupported" )
    }

    /** Override the default implementation MatrixStorage<ValueType>::assignTranspose 
     *
     *  Assign transpose is only supported if other matrix is also a stencil storage.
     *  Otherwise an exception is thrown.
     */
    virtual void assignTranspose( const MatrixStorage<ValueType>& other );

    /** Implementation of MatrixStorage::matrixTimesVector for stencil storage */

    virtual void matrixTimesVector(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    /** Implementation of MatrixStorage::jacobiIterate for stencil storage */

    virtual void jacobiIterate(
        hmemo::HArray<ValueType>& solution,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const ValueType omega ) const;

    /** Implementation of MatrixStorage::matrixTimesVectorAsync for stencil storage */

    virtual tasking::SyncToken* matrixTimesVectorAsync(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const common::MatrixOp op ) const;

    using _MatrixStorage::getNumRows;
    using _MatrixStorage::getNumColumns;
    using _MatrixStorage::prefetch;
    using _MatrixStorage::getContextPtr;

protected:

    common::Grid mGrid;                    //! grid for which this matrix storage stands
    common::Stencil<ValueType> mStencil;   //! stencil that specifies the linear mapping  with involved neighbors

private:

    static std::string initTypeName();

    tasking::SyncToken* incGEMV( hmemo::HArray<ValueType>& result, const ValueType alpha, const hmemo::HArray<ValueType>& x, bool async ) const;

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void StencilStorage<ValueType>::setRow( const scai::hmemo::HArray<ValueType>&, IndexType, scai::common::BinaryOp )
{
    COMMON_THROWEXCEPTION( "setRow unsuported" )
}

template<typename ValueType>
void StencilStorage<ValueType>::setColumn( const scai::hmemo::HArray<ValueType>&, IndexType, scai::common::BinaryOp )
{
    COMMON_THROWEXCEPTION( "setColumn unsuported" )
}

template<typename ValueType>
void StencilStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>&, hmemo::HArray<ValueType>&, const IndexType ) const
{
    COMMON_THROWEXCEPTION( "get(Sparse)Column unsupported for stencil storage" )
}

template<typename ValueType>
void StencilStorage<ValueType>::getColumn( hmemo::HArray<ValueType>&, const IndexType ) const
{
    COMMON_THROWEXCEPTION( "get(Sparse)Column unsupported for stencil storage" )
}

} /* end namespace lama */

} /* end namespace scai */
