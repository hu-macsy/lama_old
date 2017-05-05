/**
 * @file GridVector.hpp
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
 * @brief Definition of a multi-dimensional vector.
 * @author Thomas Brandes
 * @date 07.02.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/DenseVector.hpp>

// other SCAI libraries

#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/NoCommunicator.hpp>
#include <scai/common/Grid.hpp>

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
        dmemo::DistributionPtr dist( new dmemo::GridDistribution( common::Grid1D( 0 ), comm ) );
        this->allocate( dist );
    }

    GridVector( const common::Grid& grid, const ValueType initVal ) : DenseVector<ValueType>()

    {
        dmemo::CommunicatorPtr comm( new dmemo::NoCommunicator() );
        dmemo::DistributionPtr dist( new dmemo::GridDistribution( grid, comm ) );
        this->allocate( dist );
        this->assign( initVal );
    }

    const common::Grid& globalGrid() const
    {
        const dmemo::GridDistribution* dist = dynamic_cast<const dmemo::GridDistribution*>( this->getDistributionPtr().get() );
        SCAI_ASSERT_ERROR( dist, "Grid vector has no grid distribution" )
        return dist->getGlobalGrid();
    }

    IndexType nDims() const
    {
        return globalGrid().nDims();
    }

    /** Swap available array with data into a grid vector. */

    void swap( hmemo::HArray<ValueType>& data, const common::Grid& grid )
    {
        dmemo::CommunicatorPtr comm( new dmemo::NoCommunicator() );
        dmemo::DistributionPtr dist( new dmemo::GridDistribution( grid, comm ) );
   
        DenseVector<ValueType>::swap( data, dist );
    }

    /** Help class to observe the further use of operator[] in LArray */

    class GridElemProxy
    {
    public:

        /** Proxy constructed by ref to the grid and the index positions */

        GridElemProxy( GridVector& gridVector, const IndexType pos[] ) :

            mGridVector( gridVector )

        {
            // only copy relevant positions 

            for ( IndexType idim = 0; idim < gridVector.nDims(); ++idim )
            {
                mPos[idim] = pos[idim];
            }
        }

        /** indexed value proxy can be used to get its value */

        operator Scalar() const
        {
            return get();
        }

        /** indexed value proxy can be assigned a value */

        GridElemProxy& operator= ( Scalar val )
        {
            mGridVector.setValue( pos(), Scalar( val ) );
        }

        /** Override the default assignment operator to avoid ambiguous interpretation of a[i] = b[i] */

        GridElemProxy& operator= ( const GridElemProxy& other )
        {
            mGridVector.setValue( pos(), other.get() );
        }

        IndexType pos() const
        {
            return mGridVector.globalGrid().linearPos( mPos );
        }

        Scalar get() const
        {
            return mGridVector.getValue( pos() );
        }

    private:

        GridVector<ValueType>& mGridVector;        //!< keep reference to the grid vector
        IndexType mPos[SCAI_GRID_MAX_DIMENSION];   //!< holds the selected grid element
    };

    GridElemProxy operator() ( const IndexType i1, const IndexType i2 )
    {
        IndexType pos[] = { i1, i2, 0, 0, 0, 0 };
        return GridElemProxy( *this,  pos );
    }

    GridElemProxy operator() ( const IndexType i1, const IndexType i2, const IndexType i3 )
    {
        IndexType pos[] = { i1, i2, i3, 0, 0, 0 };
        return GridElemProxy( *this,  pos );
    }

    GridElemProxy operator() ( const IndexType i1, const IndexType i2, const IndexType i3, const IndexType i4 )
    {
        IndexType pos[] = { i1, i2, i3, i4, 0, 0 };
        return GridElemProxy( *this,  pos );
    }

    Scalar operator() ( const IndexType i1, const IndexType i2 ) const
    {
        const common::Grid2D& grid2D = static_cast<const common::Grid2D&>( globalGrid() );
        return DenseVector<ValueType>::getValue( grid2D.linearPos( i1, i2 ) );
    }

    Scalar operator() ( const IndexType i1, const IndexType i2, const IndexType i3 ) const
    {
        const common::Grid3D& grid3D = static_cast<const common::Grid3D&>( globalGrid() );
        return DenseVector<ValueType>::getValue( grid3D.linearPos( i1, i2, i3 ) );
    }

    Scalar operator() ( const IndexType i1, const IndexType i2, const IndexType i3, const IndexType i4 ) const
    {
        const common::Grid4D& grid4D = static_cast<const common::Grid4D&>( globalGrid() );
        return DenseVector<ValueType>::getValue( grid4D.linearPos( i1, i2, i3, i4 ) );
    }

    using DenseVector<ValueType>::operator=;

};

} /* end namespace lama */

} /* end namespace scai */
