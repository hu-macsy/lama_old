/**
 * @file StencilMatrix.hpp
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
#include "StencilStorage.hpp"
#include "Stencil.hpp"

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

    typedef ValueType MatrixValueType; //!< This is the type of the matrix values.

    /** Static method that returns the name of the matrix class. */

    static const char* typeName();

    /** Constructor, creates a replicated zero-matrix of size numRows x numColums */

    StencilMatrix( const common::Grid& grid, const Stencil<ValueType>& stencil );

    /** Constructor of a stencil matrix by a distributed grid and a stencil 
     *
     *  @param dist must be a GridDistribution 
     */
    StencilMatrix( dmemo::DistributionPtr dist, const Stencil<ValueType>& stencil );

    // Expression constructors are not allowed

    /**
     * @brief Destructor. Releases all allocated resources.
     */
    ~StencilMatrix();

    /** Override MatrixStorage<ValueType>::getLocalStorage with covariant return type. */

    virtual const StencilStorage<ValueType>& getLocalStorage() const;

    /** Override MatrixStorage<ValueType>::getHaloStorage with covariant return type. */

    virtual const CSRStorage<ValueType>& getHaloStorage() const;

    /* Implementation of pure method Matrix::newMatrix with covariant return type */

    virtual StencilMatrix<ValueType>* newMatrix() const;

    /* Implementation of pure method Matrix::copy with covariant return type */

    virtual StencilMatrix<ValueType>* copy() const;

    /* Implementation of pure method Matrix::getFormat */

    virtual Format::MatrixStorageFormat getFormat() const
    {
        return Format::CSR;
    }

    /* Implementation of pure method of class Matrix. */

    virtual const char* getTypeName() const;

    using SparseMatrix<ValueType>::setContextPtr;

protected:

    using SparseMatrix<ValueType>::mLocalData;
    using SparseMatrix<ValueType>::mHaloData;
    using SparseMatrix<ValueType>::mHalo;

private:

    static std::string initTypeName();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

};

} /* end namespace lama */

} /* end namespace scai */
