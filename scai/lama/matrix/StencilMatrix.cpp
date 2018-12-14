/**
 * @file StencilMatrix.cpp
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
 * @brief Implementation of methods and constructors for template class StencilMatrix.
 * @author Thomas Brandes
 * @date 13.04.2017
 */

// hpp
#include "StencilMatrix.hpp"

#include <scai/lama/storage/StorageMethods.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/sparsekernel/openmp/OpenMPStencilKernel.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>

using std::shared_ptr;

namespace scai
{

using common::Stencil;
using namespace dmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, StencilMatrix<ValueType>::logger,
                              "Matrix.SparseMatrix.StencilMatrix" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>::StencilMatrix() : SparseMatrix<ValueType>()

{
    SCAI_LOG_INFO( logger, "create default stencil matrix" )

    mLocalData.reset( new StencilStorage<ValueType>() );
    mHaloData.reset( new CSRStorage<ValueType>() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StencilMatrix<ValueType>::define( const common::Grid& grid, const Stencil<ValueType>& stencil )
{
    const IndexType N = grid.size();   // shape of matrix is N x N

    _Matrix::setReplicatedMatrix( grid.size(), grid.size() );

    SCAI_LOG_INFO( logger, "create stencil matrix, grid = " << grid << ", stencil = " << stencil )

    mLocalData.reset( new StencilStorage<ValueType>( grid, stencil ) );
    mHaloData.reset( new CSRStorage<ValueType>() );
    mHaloData->allocate( N, 0 );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>::StencilMatrix( const common::Grid& grid, const Stencil<ValueType>& stencil ) : 

    SparseMatrix<ValueType>()

{
    // own method needed as halo storage will be of CSR format
    define( grid, stencil );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StencilMatrix<ValueType>::buildStencilHaloStorage(
    hmemo::HArray<IndexType>& haloIA,
    hmemo::HArray<IndexType>& haloJA,
    hmemo::HArray<ValueType>& haloValues,
    const dmemo::GridDistribution& gridDist,
    const Stencil<ValueType>& stencil )
{
    // Determine all non-local grid points of the boundaries

    const common::Grid& localGrid  = gridDist.getLocalGrid();
    const common::Grid& globalGrid = gridDist.getGlobalGrid();

    const IndexType n     = gridDist.getLocalSize();
    const IndexType nDims = stencil.nDims();

    IndexType localGridDistances[ SCAI_GRID_MAX_DIMENSION ];
    IndexType globalGridDistances[ SCAI_GRID_MAX_DIMENSION ];

    localGrid.getDistances( localGridDistances );
    globalGrid.getDistances( globalGridDistances );

    // get the offsets of the stencil points in the global grid 

    std::unique_ptr<int[]> stencilOffsets( new int[ stencil.nPoints() ] );
    stencil.getLinearOffsets( stencilOffsets.get(), globalGridDistances );

    {
        hmemo::WriteOnlyAccess<IndexType> wIA( haloIA, n + 1 );

        sparsekernel::OpenMPStencilKernel::stencilHaloSizes( 
            wIA.get(), nDims, localGrid.sizes(), localGridDistances, gridDist.localLB(),
            globalGrid.sizes(), globalGrid.borders(), stencil.nPoints(), stencil.positions() );

        wIA.resize( n );
    }

    // scan of IA 

    const IndexType nnz = utilskernel::HArrayUtils::scan1( haloIA );

    {
        hmemo::WriteOnlyAccess<IndexType> wJA( haloJA, nnz );
        hmemo::WriteOnlyAccess<ValueType> wValues( haloValues, nnz );
        hmemo::ReadAccess<IndexType> rIA( haloIA );

        sparsekernel::OpenMPStencilKernel::stencilHaloCSR(
            wJA.get(), wValues.get(), rIA.get(),
            nDims, localGrid.sizes(), localGridDistances, gridDist.localLB(),
            globalGrid.sizes(), globalGridDistances, globalGrid.borders(),
            stencil.nPoints(), stencil.positions(),
            stencil.values(), stencilOffsets.get() );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StencilMatrix<ValueType>::define( dmemo::DistributionPtr dist, const Stencil<ValueType>& stencil )
{
    SCAI_ASSERT_ERROR( dist.get(), "NULL dist" )

    using namespace dmemo;

    const GridDistribution* gridDist = dynamic_cast<const GridDistribution*>( dist.get() );
  
    SCAI_ASSERT_ERROR( gridDist, "not grid distribution: dist = " << *dist )

    _Matrix::setDistributedMatrix( dist, dist );

    // const common::Grid& globalGrid = gridDist->getGlobalGrid();
    const common::Grid& localGrid  = gridDist->getLocalGrid();

    mLocalData.reset( new StencilStorage<ValueType>( localGrid, stencil ) );

    CSRStorage<ValueType>* haloStorage = new CSRStorage<ValueType>();
    mHaloData.reset( haloStorage );   

    hmemo::HArray<IndexType> haloIA;
    hmemo::HArray<IndexType> haloJA;
    hmemo::HArray<ValueType> haloValues;

    // ToDo: build halo storage and halo 

    buildStencilHaloStorage( haloIA, haloJA, haloValues, *gridDist, stencil );
    
    SCAI_LOG_DEBUG( logger, "haloStorage, haloIA = " << haloIA )
    SCAI_LOG_DEBUG( logger, "haloStorage, haloJA = " << haloJA )

    IndexType haloNumColumns;

    _StorageMethods::buildHalo( mHaloPlan, haloJA, haloNumColumns, *gridDist );

    SCAI_LOG_DEBUG( logger, "haloPlan = " << mHaloPlan )

    haloStorage->setCompressThreshold( 0.5 );
    haloStorage->allocate( haloIA.size() - 1, haloNumColumns );
    haloStorage->swap( haloIA, haloJA, haloValues );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>::StencilMatrix( dmemo::DistributionPtr dist, const Stencil<ValueType>& stencil )

    : SparseMatrix<ValueType>()

{
    define( dist, stencil );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>::~StencilMatrix()
{
    SCAI_LOG_INFO( logger, "~StencilMatrix" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const StencilStorage<ValueType>& StencilMatrix<ValueType>::getLocalStorage() const
{
    // here we need a dynamic cast as for any stupid reasons somebody
    // has modified the underlying storage type
    const StencilStorage<ValueType>* local = dynamic_cast<const StencilStorage<ValueType>*>( this->mLocalData.get() );
    SCAI_ASSERT_ERROR( local, "StencilMatrix: local storage is no more Stencil: " << *this->mLocalData )
    return *local;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const CSRStorage<ValueType>& StencilMatrix<ValueType>::getHaloStorage() const
{
    const CSRStorage<ValueType>* halo = dynamic_cast<const CSRStorage<ValueType>*>( this->mHaloData.get() );
    SCAI_ASSERT_ERROR( halo, "StencilMatrix: local storage is no more CSR: " << *this->mHaloData )
    return *halo;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>* StencilMatrix<ValueType>::newMatrix() const
{
    return new StencilMatrix<ValueType>();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const Stencil<ValueType>& StencilMatrix<ValueType>::getStencil() const
{
    SCAI_ASSERT_ERROR( mLocalData.get(), "no local data available" )

    const StencilStorage<ValueType>* localStencilStorage = dynamic_cast<const StencilStorage<ValueType>*>( mLocalData.get() );

    SCAI_ASSERT_ERROR( localStencilStorage, "local data no stencil storage: " << mLocalData )

    return localStencilStorage->getStencil();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>::StencilMatrix( const StencilMatrix<ValueType>& other ) :

    SparseMatrix<ValueType>( other )

{
    // copy constructor of sparse matrix should have already done all 

    SCAI_LOG_INFO( logger, "copy constructed stencil matrix. " << *this )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>* StencilMatrix<ValueType>::copy() const
{
    return new StencilMatrix<ValueType>( *this );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>::StencilMatrix( StencilMatrix&& other ) noexcept :

    SparseMatrix<ValueType>()
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>&  StencilMatrix<ValueType>::operator=( const StencilMatrix& other )
{
    SparseMatrix<ValueType>::operator=( other );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
StencilMatrix<ValueType>&  StencilMatrix<ValueType>::operator=( StencilMatrix&& other )
{
    SparseMatrix<ValueType>::operator=( std::move( other ) );
    return *this;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StencilMatrix<ValueType>::buildLocalStorage( _MatrixStorage& storage ) const
{
    if ( getColDistribution().isReplicated() )
    {
        // copy local storage with format / value conversion
        storage = *mLocalData;
    }
    else
    {
        // temporary local storage with joined columns needed before

        if ( storage.getValueType() == getValueType() )
        {
            MatrixStorage<ValueType>& typedStorage = static_cast<MatrixStorage<ValueType>&>( storage );
            typedStorage.joinHalo( *mLocalData, *mHaloData, mHaloPlan, getColDistribution() );
        }     
        else
        {
            CSRStorage<ValueType> tmp; 
            tmp.joinHalo( *mLocalData, *mHaloData, mHaloPlan, getColDistribution() );
            storage = tmp;
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
const char* StencilMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
std::string StencilMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "StencilMatrix<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* StencilMatrix<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template specializations and nstantiations                          */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( StencilMatrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
