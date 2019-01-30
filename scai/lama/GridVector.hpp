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

    GridVector() : DenseVector<ValueType>()

    {
        // give it a default grid distribution for consistency

        dmemo::CommunicatorPtr comm( new dmemo::NoCommunicator() );
        dmemo::DistributionPtr dist( new dmemo::GridDistribution( common::Grid(), comm ) );
        this->allocate( dist );
    }

    /** Read in a grid from an input file. 
     *
     *  Note: suffix of filename decides about file format. 
     */
    GridVector( const std::string& inputFileName );

    GridVector( const common::Grid& grid ) : DenseVector<ValueType>()
    {
        dmemo::CommunicatorPtr comm( new dmemo::NoCommunicator() );
        dmemo::DistributionPtr dist( new dmemo::GridDistribution( grid, comm ) );
        this->allocate( dist );
    }

    GridVector( const common::Grid& grid, const ValueType initVal ) : DenseVector<ValueType>()
    {
        dmemo::CommunicatorPtr comm( new dmemo::NoCommunicator() );
        dmemo::DistributionPtr dist( new dmemo::GridDistribution( grid, comm ) );
        this->allocate( dist );
        this->setScalar( initVal );
    }

    GridVector( dmemo::DistributionPtr dist ) : DenseVector<ValueType>()
    {
        this->allocate( dist );
    }

    GridVector( dmemo::DistributionPtr dist, const ValueType initVal ) : DenseVector<ValueType>()
    {
        this->allocate( dist );
        this->setScalar( initVal );
    }

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

    /** Override of DenseVector::writeLocalToFile 
     *
     *  In contrary to DenseVector this method tries also to save the information about the grid shape
     */
    virtual void writeLocalToFile(
        const std::string& fileName,
        const std::string& fileType,
        const common::ScalarType dataType,
        const FileMode fileMode ) const;

    /** Override of _Vector::readLocalFromFile 
     *
     *  This method reads the local part of a grid vector and tries to read also the shape.
     */
    virtual IndexType readLocalFromFile( const std::string& fileName, const IndexType first = 0, const IndexType size = invalidIndex );
};

} /* end namespace lama */

} /* end namespace scai */
