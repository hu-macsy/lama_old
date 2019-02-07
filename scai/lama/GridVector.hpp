/**
 * @file GridVector.hpp
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
 * @brief Definition of a multi-dimensional vector.
 * @author Thomas Brandes
 * @date 07.02.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/GridSection.hpp>

// other SCAI libraries

#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/NoCommunicator.hpp>
#include <scai/common/Grid.hpp>
#include <scai/common/BinaryOp.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Multidimensioal vector is a dense vector that gives the data a grid shape
 * 
 * \code
 * DenseVector<double> x( nx * ny );
 * GridVector<double> x( nx, ny );
 * Grid grid( nx, ny );
 * DistributionPtr gridDist( new GridDistribution( grid, comm ) );
 * GridVector<double> x( gridDist );
 *
 * LGArray<double>& localX = x.getLocalData();
 * localX( 5, 3 ) = localX( 4, 2 );
 * \endcode
 *
 */

/**
 * @brief The GridVector represents a distributed n-dimensional grid with elements of type ValueType.
 *
 * @tparam ValueType the value type for the vector values.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT GridVector:

    public DenseVector<ValueType>

{
public:

    /** Default vector to enable just declaration of a grid vector. */

    explicit GridVector( hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    /**
     * @brief Create a replicated GridVector of a certain size initialized with same value for all elements.
     *
     * @param[in] grid  shape of the grid vector
     * @param[in] value the value to assign to all elements of the new GridVecto
     * @param[in] context   specifies optionally the context where grid vector should reside
     */
    GridVector( const common::Grid& grid, const ValueType value, hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    GridVector( dmemo::DistributionPtr dist, const ValueType value, hmemo::ContextPtr context = hmemo::Context::getContextPtr() );

    const common::Grid& globalGrid() const
    {
        const dmemo::GridDistribution* dist = dynamic_cast<const dmemo::GridDistribution*>( this->getDistributionPtr().get() );
        SCAI_ASSERT_ERROR( dist, "Grid vector has no grid distribution" )
        return dist->getGlobalGrid();
    }

    const common::Grid& localGrid() const
    {
        const dmemo::GridDistribution* dist = dynamic_cast<const dmemo::GridDistribution*>( this->getDistributionPtr().get() );
        SCAI_ASSERT_ERROR( dist, "Grid vector has no grid distribution" )
        return dist->getLocalGrid();
    }

    const IndexType* localLB() const
    {
        const dmemo::GridDistribution* dist = dynamic_cast<const dmemo::GridDistribution*>( this->getDistributionPtr().get() );
        SCAI_ASSERT_ERROR( dist, "Grid vector has no grid distribution" )
        return dist->localLB();
    }

    IndexType nDims() const
    {
        return globalGrid().nDims();
    }

    IndexType size( const IndexType dim ) const
    {
        return globalGrid().size( dim );
    }

    /** Make size() of Vector visible */

    using _Vector::size;

    /** Swap available array with data into a grid vector. */

    void swap( hmemo::HArray<ValueType>& data, const common::Grid& grid )
    {
        dmemo::CommunicatorPtr comm( new dmemo::NoCommunicator() );
        dmemo::DistributionPtr dist( new dmemo::GridDistribution( grid, comm ) );
   
        DenseVector<ValueType>::swap( data, dist );
    }

    /** Swap with another GridVector */

    void swap( GridVector<ValueType>& other )
    {
        DenseVector<ValueType>::swap( other );
    }

    using _Vector::swap;

    /** Set this grid vector by another grid vector that is reduced in one dimension */

    void reduce( const GridVector<ValueType>& other, IndexType dim, const common::BinaryOp redOp );

    /** _Matrix-matrix multiplication : this += alpha * v1 * v2
     *
     *  @param alpha is just an additional scaling factor to avoid temporaries
     *  @param v1 must be a two-dimensional grid-vector with v1.size( 0 ) == this->size( 0 )
     *  @param v2 must be a two-dimensional grid-vector with v2.size( 0 ) == v1.size( 1 ) and v2.size( 1 ) == this->size( 1 )
     */

    void gemm( const ValueType alpha, const GridVector<ValueType>& v1, const GridVector<ValueType>& v2 );

    /** Set a diagonal on a two-dimensional grid vector. 
     *
     *  @param[in] diagonal must be a one-dimensional section
     *  @param[in] diagonalNumber specifies the diagonal, 0 stands for main diagonal
     *
     *  This routine throws an exception if the dimensions are not correct.
     *  This routine throws an exception if the sizes do not fit.
     */
    void setDiagonal( const GridSection<ValueType>& diagonal, const int diagonalNumber );

    using DenseVector<ValueType>::operator=;

    void operator=( const GridSection<ValueType>& other );

    GridSection<ValueType> operator() ( const Range& r1 )
    {
        return GridSection<ValueType>( *this, r1 );
    }

    GridSection<ValueType> operator() ( const Range& r1, const Range& r2 )
    {
        return GridSection<ValueType>( *this, r1, r2 );
    }

    GridSection<ValueType> operator() ( const Range& r1, const Range& r2, const Range& r3 )
    {
        return GridSection<ValueType>( *this, r1, r2, r3 );
    }

    GridSection<ValueType> operator() ( const Range& r1, const Range& r2, const Range& r3, const Range& r4 )
    {
        return GridSection<ValueType>( *this, r1, r2, r3, r4 );
    }

    /** Implementation of _Vector::writeLocalToFile */

    virtual void writeLocalToFile( FileIO& file ) const;

    /** Implementation of _Vector::readFromFile */

    virtual void readFromFile( FileIO& file );

    using _Vector::readFromFile;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

};

} /* end namespace lama */

} /* end namespace scai */
