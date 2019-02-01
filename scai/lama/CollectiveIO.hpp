/**
 * @file CollectioveIO.hpp
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
 * @brief Parallel I/O operations for LAMA matrices and vectors using CollectiveFile.
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#pragma once

#include <scai/dmemo/CollectiveFile.hpp>
#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
class Vector;

template<typename ValueType>
class DenseVector;

template<typename ValueType>
class SparseVector;

template<typename ValueType>
class Matrix;

template<typename ValueType>
class DenseMatrix;

template<typename ValueType>
class SparseMatrix;

template<typename ValueType>
class CSRSparseMatrix;

/**
 *  Class provides static methods for read and write operations of LAMA data structures 
 */
class CollectiveIO
{
public:

    /**
     *  @brief Write a vector into a collective file.
     *  
     *  @param[in] file is the collecitve file (must have been opened for write)
     *  @param[in] vector is the vector to be written
     *  @param[in] fileIndexType specifies the type for writing index type values
     *  @param[in] fileData specifies the type for writing data entries into the file
     *  
     */
    template<typename ValueType>
    static void write( 
        dmemo::CollectiveFile& file, 
        const Vector<ValueType>& vector,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileDataType );

    /**
     *  @brief Write a vector into a collective file.
     *  
     *  @param[in] file is the collecitve file (must have been opened for write)
     *  @param[in] vector is the vector to be written
     *  
     *  This method will write values without any conversion into the file.  
     */
    template<typename ValueType>
    static void write( dmemo::CollectiveFile& file, const Vector<ValueType>& vector );

    /**
     *  @brief Write a vector into a collective file.
     *  
     *  @param[in] file is the collecitve file (must have been opened for write)
     *  @param[in] vector is the vector to be written
     *  @param[in] fileIndexType specifies the type for writing index type values
     *  @param[in] fileData specifies the type for writing data entries into the file
     *  
     */
    template<typename ValueType>
    static void write(
        dmemo::CollectiveFile& file,
        const Matrix<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileDataType );

    /**
     *  @brief Write a matrix into a collective file.
     *  
     *  @param[in] file is the collecitve file (must have been opened for write)
     *  @param[in] matrix is the matrix to be written
     *  
     *  This method will write values without any conversion into the file.  
     */
    template<typename ValueType>
    static void write( dmemo::CollectiveFile& file, const Matrix<ValueType>& matrix );

    /**
     *  @brief Read a dense vector from a collective file.
     */
    template<typename ValueType>
    static void read( dmemo::CollectiveFile& file, Vector<ValueType>& vector );

    /**
     *  @brief Read a sparse matrix from a collective file.
     */
    template<typename ValueType>
    static void read( dmemo::CollectiveFile& file, Matrix<ValueType>& matrix );

    /**
     *  @brief Get the identification for a dense vector in a collective file.
     */
    static int getDenseVectorId();

    /**
     *  @brief Get the identification for a sparse matrix (CSR) in a collective file.
     */
    static int getSparseMatrixId();

private:

    template<typename ValueType>
    static void writeDenseMatrix( 
        dmemo::CollectiveFile& file, 
        const DenseMatrix<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileValueType );

    template<typename ValueType>
    static void writeIt( 
        dmemo::CollectiveFile& file, 
        const DenseMatrix<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileValueType );

    template<typename ValueType>
    static void writeDenseVector( 
        dmemo::CollectiveFile& file, 
        const DenseVector<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileValueType );

    template<typename ValueType>
    static void writeIt( 
        dmemo::CollectiveFile& file, 
        const DenseVector<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileValueType );

    template<typename ValueType>
    static void writeSparseMatrix( 
        dmemo::CollectiveFile& file, 
        const SparseMatrix<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileValueType );

    template<typename ValueType>
    static void writeIt( 
        dmemo::CollectiveFile& file, 
        const CSRSparseMatrix<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileValueType );

    template<typename ValueType>
    static void writeSparseVector( 
        dmemo::CollectiveFile& file, 
        const SparseVector<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileValueType );

    template<typename ValueType>
    static void writeIt( 
        dmemo::CollectiveFile& file, 
        const SparseVector<ValueType>& matrix,
        const common::ScalarType fileIndexType,
        const common::ScalarType fileValueType );

    template<typename ValueType>
    static void readIt( dmemo::CollectiveFile& file, CSRSparseMatrix<ValueType>& matrix );

    template<typename ValueType>
    static void readIt( dmemo::CollectiveFile& file, DenseMatrix<ValueType>& matrix );

    template<typename ValueType>
    static void readIt( dmemo::CollectiveFile& file, SparseVector<ValueType>& matrix );

    template<typename ValueType>
    static void readIt( dmemo::CollectiveFile& file, DenseVector<ValueType>& matrix );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

};

}  // namespace lama

}  // namespace scai
