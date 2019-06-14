/** Example program for building a hybrid matrix. */

#include <scai/lama.hpp>

#include <scai/lama/matrix/StencilMatrix.hpp>
#include <scai/lama/matrix/HybridMatrix.hpp>
#include <scai/lama/matrix/MatrixAssembly.hpp>

#include <scai/common/Grid.hpp>
#include <scai/common/Stencil.hpp>

using namespace scai;
using namespace lama;

/** Function that returns a smooth matrix for Full Wavefield Migration.
 *
 *  It is very close to a stencil matrix but updates some elements. So it 
 *  returns CSR sparse matrix.
 */
template<typename ValueType>
CSRSparseMatrix<ValueType> createSmoothMatrix( 
    const common::Grid3D& grid,
    const ValueType strength,
    const ValueType factor[3] )
{
    const IndexType nx = grid.size( 0 );
    const IndexType ny = grid.size( 1 );
    const IndexType nz = grid.size( 2 );

    const ValueType factorx = factor[0];
    const ValueType factory = factor[1];
    const ValueType factorz = factor[2];

    MatrixAssembly<ValueType> assembly;

    auto gridDist = std::make_shared<dmemo::GridDistribution>( grid );

    IndexType central_index;
    IndexType entriesinline;
    ValueType denominator;
    ValueType norm;

    for ( IndexType i( 0 ); i < nx; i++ )
    {
        IndexType itmp( i * ny * nz );

        for ( IndexType j( 0 ); j < ny; j++ )
        {
            IndexType jtmp( j * nz );

            for ( IndexType k( 0 ); k < nz; k++ )
            {
                central_index = itmp + jtmp + k ;

                // we have to make sure that only one processor fills one line  

                if ( !gridDist->isLocal( central_index ) )
                {
                    continue;
                }

                // calculate denominator
                denominator = factorx + factory + factorz; //there is always a minor diagonal(z) and two far away diagonal(x,y)
                entriesinline = 4;                         //there are always at least 4 elements in a row

                if ( ( k > 0 ) && ( k < ( nz - 1 ) ) )
                {
                    denominator += factorz;
                    entriesinline++;
                }

                if ( ( j > 0 ) && ( j < ( ny - 1 ) ) )
                {
                    denominator += factory;
                    entriesinline++;
                }

                if ( ( i > 0 ) && ( i < ( nx - 1 ) ) )
                {
                    denominator += factorx;
                    entriesinline++;
                }

                norm = static_cast<ValueType>( entriesinline ) / static_cast<ValueType>( 7 );

                if ( k > 0 )
                {
                    assembly.push( central_index, central_index - 1, factorz * strength * norm / denominator );
                }

                if ( k < ( nz - 1 ) )
                {
                    assembly.push( central_index, central_index + 1, factorz * strength * norm / denominator );
                }

                if ( j > 0 )
                {
                    assembly.push( central_index, central_index - nz, factory * strength * norm / denominator );
                }

                if ( j < ( ny - 1 ) )
                {
                    assembly.push( central_index, central_index + nz, factory * strength * norm / denominator );
                }

                if ( i > 0 )
                {
                    assembly.push( central_index, central_index - nz * ny, factorx * strength * norm / denominator );
                }

                if ( i < ( nx - 1 ) )
                {
                    assembly.push( central_index, central_index + nz * ny, factorx * strength * norm / denominator );
                }

                assembly.push( central_index, central_index, - strength * norm );
            }
        }
    }

    auto smoothMatrix = zero<CSRSparseMatrix<ValueType>>( gridDist, gridDist );

    // filling will be very fast as there are no elements to exchange

    smoothMatrix.fillFromAssembly( assembly );

    return smoothMatrix;
}

/** This method creates a stencil as used for the smooth matrix. 
 *
 *  As we have no topology/grid we have always the same norm and denominator.
 */
template<typename ValueType>
common::Stencil3D<ValueType> createSmoothStencil( 
    const ValueType strength,
    const ValueType factor[3] )
{
    const ValueType factorx = factor[0];
    const ValueType factory = factor[1];
    const ValueType factorz = factor[2];

    common::Stencil3D<ValueType> stencil;

    ValueType denominator = factor[0] + factor[1] + factor[2] + factor[0] + factor[1] + factor[2];

    stencil.addPoint( -1,  0,  0, factorx * strength / denominator );
    stencil.addPoint(  1,  0,  0, factorx * strength / denominator );
    stencil.addPoint(  0, -1,  0, factory * strength / denominator );
    stencil.addPoint(  0,  1,  0, factory * strength / denominator );
    stencil.addPoint(  0,  0, -1, factorz * strength / denominator );
    stencil.addPoint(  0,  0,  1, factorz * strength / denominator );
    stencil.addPoint(  0,  0,  0, -strength );

    return stencil;
}

int main( )
{
    typedef double ValueType;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    const common::Grid3D grid( 100, 150, 120 );

    const ValueType factor[] = { 0.8, 1.2, 1.1 };

    const ValueType strength = 8.0;

    CSRSparseMatrix<ValueType> csrSmoothMatrix = createSmoothMatrix<ValueType>( grid, strength, factor );

    dmemo::DistributionPtr gridDistribution = csrSmoothMatrix.getRowDistributionPtr();

    // Define stencil as in createSmoothMatrix, but without norm and same denominator

    common::Stencil3D<ValueType> stencil = createSmoothStencil<ValueType>( strength, factor );

    StencilMatrix<ValueType> stencilMatrix( gridDistribution, stencil );

    auto diffMatrix = eval<DIASparseMatrix<ValueType>>( csrSmoothMatrix - stencilMatrix );

    // due to rounding errors the difference matrix might contain many nearly zero elements

    diffMatrix.compress( 1e-8, false );   // do not keep zero diagonal elements

    const IndexType nnzAll = csrSmoothMatrix.getNumValues();
    const IndexType nnz = diffMatrix.getNumValues();
    const double percentage = nnz * 100.0 / nnzAll;

    std::cout << *comm << ": Diff matrix contains " << nnz << " of " << nnzAll << " values, is " << percentage << "%" << std::endl;

    HybridMatrix<ValueType> hybridSmoothMatrix( stencilMatrix, diffMatrix );

    ValueType delta = ValueType( 1 ) / ValueType( grid.size() );

    // create an example vector

    DenseVector<ValueType> x = denseVectorLinear<ValueType>( gridDistribution, delta, delta );

    // Note: csrSmoothMatrix.purge() would free all used memory, but we keep it for verification of result

    //  matrix-vector multiplication should give 'nearly' same result for csr and hybrid matrix

    auto y1 = eval<DenseVector<ValueType>>( csrSmoothMatrix * x );
    auto y2 = eval<DenseVector<ValueType>>( hybridSmoothMatrix * x );

    std::cout << "norm diff = " << y1.maxDiffNorm( y2 ) << std::endl;
}
