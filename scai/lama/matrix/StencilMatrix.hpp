/**
 * @file lama/matrix/StencilMatrix.hpp
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
 * @brief Definition of matrix class for stencils on distributed grids
 * @author Thomas Brandes
 * @date 13.04.2017
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/SparseMatrix.hpp>

// local library
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/StencilStorage.hpp>
#include <scai/common/Stencil.hpp>

#include <scai/dmemo/GridDistribution.hpp>

namespace scai
{

namespace lama
{

/** Definition of a derived class for SparseMatrix that uses the stencil storage 
 *  for the local part and the CSR format for the halo part.
 *
 *  This matrix class will not register at the factory.
 *
 *  In contrary to the other sparse/dense matrix format, the following operations
 *  are unsupported:
 *
 *  - cannot convert other matrix formats to a stencil matrix, but a stencil matrix
 *    can be easily converted to all other formats
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT StencilMatrix:

    public SparseMatrix<ValueType>

{

public:

    /** Static method that returns the name of the matrix class. */

    static const char* typeName();

    /** Default constructor, creates a replicated matrix of size 0 x 0 */

    StencilMatrix();

    void define( const common::Grid& grid, const common::Stencil<ValueType>& stencil );

    void define( dmemo::DistributionPtr dist, const common::Stencil<ValueType>& stencil );

    /** Constructor, creates a replicated zero-matrix of size numRows x numColums */

    StencilMatrix( const common::Grid& grid, const common::Stencil<ValueType>& stencil );

    /** Constructor of a stencil matrix by a distributed grid and a stencil 
     *
     *  @param dist must be a GridDistribution 
     *  @param stencil specifies operation on the distributed grid points
     *
     *  The dimensions of the stencil must match the number of dimensions of the grid.
     *
     *  Note: the stencil matrix is like a sparse matrix where the local part is a stencil storage
     *        and the halo part is given by a CSR storage.
     */
    StencilMatrix( dmemo::DistributionPtr dist, const common::Stencil<ValueType>& stencil );

    /**
     * @brief Destructor. Releases all allocated resources.
     */
    ~StencilMatrix();

    /** Override the default assignment operator that would not make deep copies. */

    StencilMatrix& operator=( const StencilMatrix& matrix );

    /** Override the default move assignment operator */

    StencilMatrix& operator=( StencilMatrix&& matrix );

    /** Copy constructor */

    StencilMatrix( const StencilMatrix& other );

    /** Move constructor */

    StencilMatrix( StencilMatrix&& other ) noexcept;

    /** Override MatrixStorage<ValueType>::getLocalStorage with covariant return type. */

    virtual const StencilStorage<ValueType>& getLocalStorage() const;

    /** Override MatrixStorage<ValueType>::getHaloStorage with covariant return type. */

    virtual const CSRStorage<ValueType>& getHaloStorage() const;

    /* Implementation of pure method _Matrix::newMatrix with covariant return type */

    virtual StencilMatrix<ValueType>* newMatrix() const;

    /* Implementation of pure method _Matrix::copy with covariant return type */

    virtual StencilMatrix<ValueType>* copy() const;

    /* Implementation of pure method _Matrix::getFormat */

    virtual Format getFormat() const
    {
        return Format::CSR;
    }

    /** Implementation of _Matrix::buildLocalStorage, overrides SparseMatrix<ValueType>::buildLocalStorage 
     *
     *  This method gets an own implementatin as it fails for SparseMatrix and as it might be implemented
     *  more efficiently.
     */
    virtual void buildLocalStorage( _MatrixStorage& storage ) const;

    /* Query the stencil that has been used to define the stencil matrix */

    const common::Stencil<ValueType>& getStencil() const;

    /* Implementation of pure method of class _Matrix. */

    virtual const char* getTypeName() const;

    using _Matrix::getNumRows;
    using _Matrix::getNumColumns;
    using _Matrix::setIdentity;

    using _Matrix::getRowDistribution;
    using _Matrix::getRowDistributionPtr;
    using _Matrix::getColDistribution;
    using _Matrix::getColDistributionPtr;


    using Matrix<ValueType>::getValueType;
    using SparseMatrix<ValueType>::operator=;
    using SparseMatrix<ValueType>::operator-=;
    using SparseMatrix<ValueType>::operator+=;

    using SparseMatrix<ValueType>::setContextPtr;
    using SparseMatrix<ValueType>::redistribute;

protected:

    using SparseMatrix<ValueType>::mLocalData;
    using SparseMatrix<ValueType>::mHaloData;
    using SparseMatrix<ValueType>::mHalo;

private:

    static std::string initTypeName();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static void buildStencilHaloStorage(
        hmemo::HArray<IndexType>& haloIA,
        hmemo::HArray<IndexType>& haloJA,
        hmemo::HArray<ValueType>& haloValues,
        const dmemo::GridDistribution& gridDist,
        const common::Stencil<ValueType>& stencil );

};

} /* end namespace lama */

} /* end namespace scai */
