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

namespace scai
{

namespace lama
{

template<typename ValueType>
class DenseVector;

/**
 *  Class provides static methods for read and write operations of LAMA data structures 
 */
class CollectiveIO
{
public:

    /**
     *  @brief Write a dense vector into a collective file.
     */
    template<typename ValueType>
    static void write( dmemo::CollectiveFile& file, const DenseVector<ValueType>& vector );

    /**
     *  @brief Read a dense vector from a collective file.
     */
    template<typename ValueType>
    static void read( dmemo::CollectiveFile& file, DenseVector<ValueType>& vector );

    /**
     *  @param get the identification for a dense vector in a collective file.
     */
    static int getDenseVectorId();
};

}  // namespace lama

}  // namespace scai
