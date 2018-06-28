/**
 * @file common/Stencil.hpp
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
 * @brief Definition of stencil classes
 * @author Thomas Brandes
 * @date 13.04.2017
 */

#pragma once

#include <scai/common/config.hpp>
#include <scai/common/Utils.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>

#include <memory>

namespace scai
{

namespace common
{

/** Common base class for all stencils. 
 *
 *  A stencil describes a linear mapping in an n-dimensional grid how a function value is computed
 *  by the values of neighbors in a certain range. The linear mapping is the same for each grid point
 *  as long as all neighbors are available.
 *
 *  If a neighbored element is not available, different strategies are possible:
 *
 *   - its value is assumed to be zero 
 *   - the value of the other side is taken( circular boundaries )
 *   - the value of other neighboris taken ( reflecting boundaries )
 *  
 *  How this is handled is not part of this class. It must be specified when a stencil
 *  is applied to a grid.
 * 
 *  A stencil does not contain any information about the grid size
 *  only about its dimenison.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Stencil
{
public:

    /** Construct an empty stencil for a n-dimensional grid */

    Stencil( const IndexType n );

    /** Return number of dimensions for the stencil */

    IndexType nDims() const;

    /** Return number of points for the stencil */

    IndexType nPoints() const;

    /** Reserve storage for a certain number of stencil points. */

    void reserve( const IndexType npoints );

    /** Add a stencil point  
     *
     * @param[in] relpos array with nDims values to specify stencil point distance
     * @param[in] val is the value to be used for scaling with the stencil point
     */
    void addPoint( const int relpos[], const ValueType val );

    /** Return the total number of stencil points. */

    IndexType size() const;

    /** Get stencil point offsets in linearized grid 
     *  
     *  @param[in] gridDistance contains for each dim the distance between two neighbored elements in that dim
     *  @param[out] offset contains fore each stencil point the relative offset in linearized grid
     *
     *  \code
     *      common::Stencil3D stencil( 7 );  
     *      common::Grid3D grid( N1, N2, N3 );
     *      IndexType distance[3];
     *      grid.getDistances( distance ); //   returns { N2 * N3, N3, 1 }
     *      int offset[ 7 ];
     *      stencil.getLinearOffsets( offset, distance ); 
     *      // now:  offsets = { 0, 1, -1, -N3, N3, -N2*N3, N2*N3 }
     *  \endcode
     */
    inline void getLinearOffsets( int offset[], const IndexType gridDistance[] ) const;

    /** Get valid positions */

    inline IndexType getValidPoints( bool valid[], const IndexType gridSizes[], const IndexType gridPos[] ) const;

    /** Return widths of a stencil in each dimension, both directions
     *
     *  @param[out] width result array with nDims * 2 entries
     *
     *  Biggest distance to the left side for a dimension i is stored in width[2*i]
     *  biggest distance to the right side is at width[2*i+1]
     */
    inline void getWidth( IndexType width[] ) const;

    /** Return number of entries required for the stencil matrix. 
     * 
     *  @returns ( ub[0] + lb[0] + 1 ) * ( ub[1] + lb[1] + 1 ) .... 
     */
    inline IndexType getMatrixSize() const;

    /** Set the values of the stencil matrix */

    inline void getMatrix( ValueType matrix[] ) const;

    /** Compares two stencils for equality. */

    bool operator==( const Stencil& other ) const;

    /** Compares two stencils for inequality. */

    bool operator!=( const Stencil& other ) const;

    /** Get pointer to the array with stencil positions, size is nPoints() * nDims() */

    const int* positions() const;

    /** Get pointer to the array with stencil values, size is nPoints() */

    const ValueType* values() const;

    /** Set this stencil by transpose of other stencil. */

    void transpose( const Stencil& other );

    /** Scale this stencily by a certain value. */

    void scale( const ValueType scaling );

protected:

    std::vector<ValueType> mValues;
    std::vector<int> mPositions;

private:

    Stencil();

    IndexType mNDims;
};

/*
 * Output of Stencil in stream by writing its extensions
 *
 * \code
 *    cout << stencil;   // e.g. prints Stencil3D( 27 )
 * \endcode
 */
template<typename ValueType>
COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const Stencil<ValueType>& stencil  )
{
    stream << "Stencil" << stencil.nDims() << "D<" << common::TypeTraits<ValueType>::id() << "> ( " << stencil.nPoints() << " )";

    return stream;
}

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline constructors                                             */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil<ValueType>::Stencil( const IndexType n )
{
    mNDims = n;
}

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline methods                                                  */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil<ValueType>::reserve( const IndexType npoints )
{
    mValues.reserve( npoints ); 
    mPositions.reserve( npoints * mNDims ); 
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType Stencil<ValueType>::nPoints() const
{
    SCAI_ASSERT_EQ_DEBUG( mNDims * mValues.size(), mPositions.size(), "serious size mismatch" )

    return mValues.size();
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType Stencil<ValueType>::nDims() const
{
    return mNDims;
}

template<typename ValueType>
const int* Stencil<ValueType>::positions() const
{
    return &mPositions[0];
}

template<typename ValueType>
const ValueType* Stencil<ValueType>::values() const
{
    return &mValues[0];
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil<ValueType>::addPoint( const int relpos[], const ValueType val )
{
    // ToDo: check if stencil point has already been defined 

    mValues.push_back( val );

    for ( IndexType i = 0; i < mNDims; ++i )
    {
        mValues.push_back( relpos[i] );
    }
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil<ValueType>::getLinearOffsets( int offset[], const IndexType gridDistance[] ) const
{
     // compute one offset for each stencil point

     for ( IndexType k = 0; k < nPoints(); ++k )
     {
         offset[k] = 0;
         
         for ( IndexType i = 0; i < mNDims; ++i )
         {
             offset[k] += mPositions[ k * mNDims + i ] * static_cast<int>( gridDistance[ i ] );
         }
     }
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil<ValueType>::getWidth( IndexType width[] ) const
{
    // initialize width array with 0

    for ( IndexType i = 0; i < 2 * mNDims; ++i )
    {
        width[i] = 0;
    }

    for ( IndexType k = 0; k < nPoints(); ++k )
    {
         for ( IndexType i = 0; i < mNDims; ++i )
         {
             int pos = mPositions[ k * mNDims + i ];

             if ( pos < 0 )
             {
                 IndexType lbi = static_cast<IndexType>( -pos );

                 if ( lbi > width[2 * i ] )
                 { 
                     width[2 * i] = lbi;    // new max value here
                 }
             }
          
             if ( pos > 0 )
             {
                 IndexType ubi = static_cast<IndexType>( pos );

                 if ( ubi > width[2 * i + 1] )
                 { 
                     width[2 * i + 1] = ubi;    // new max value here
                 }
             }
         }
    }
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType Stencil<ValueType>::getMatrixSize() const
{
    std::unique_ptr<IndexType[]> width( new IndexType[ 2 * mNDims ] );
  
    getWidth( width.get() );
 
    IndexType size = 1;
 
    for ( IndexType i = 0; i < mNDims; ++i )
    {
        size *= width[2 * i] + width[2 * i + 1] + 1;
    }

    return size;
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil<ValueType>::getMatrix( ValueType matrix[] ) const
{
    std::unique_ptr<IndexType[]> width( new IndexType[ 2 * mNDims ] );
  
    getWidth( width.get() );
 
    IndexType size = 1;
 
    for ( IndexType i = 0; i < mNDims; ++i )
    {
        size *= width[2 * i] + width[2 * i + 1] + 1;
    }

    // Initialize matrix with 0 

    for ( IndexType matrixPos = 0; matrixPos < size; ++matrixPos )
    {
        matrix[matrixPos] = ValueType( 0 );
    }

    for ( IndexType k = 0; k < nPoints(); ++k )
    {
        IndexType matrixPos = 0;

         for ( IndexType i = 0; i < mNDims; ++i )
         {
             int pos = mPositions[ k * mNDims + i ];

             IndexType kDim = width[2 * i];

             if ( pos < 0 )
             {
                 kDim -= static_cast<IndexType>( -pos );
             }
             else
             {
                 kDim += static_cast<IndexType>( pos );
             }

             if ( i == 0 )
             {
                 matrixPos = kDim;
             }
             else
             {
                 matrixPos = matrixPos * ( width[2 * i] + width[2 * i + 1] + 1 ) + kDim;
             }
         }

         matrix[ matrixPos ] = mValues[k];
    }
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
bool Stencil<ValueType>::operator==( const Stencil<ValueType>& other ) const
{
    // stencils are never the same if they have different dims

    if ( mNDims != other.mNDims )
    {
        return false;
    }

    // now check for equal width in all directions

    std::unique_ptr<IndexType[]> width( new IndexType[ 4 * mNDims ] );

    IndexType* width1 = &width[0];
    IndexType* width2 = &width[2 * mNDims];

    getWidth( width1 );
    other.getWidth( width2 );

    IndexType size = 1;  // determine also the size of stencil matrix for later

    for ( IndexType i = 0; i < mNDims; ++i ) 
    {
        if ( width1[2 * i] != width2[2 * i] ) return false;
        if ( width1[2 * i + 1] != width2[2 * i + 1] ) return false;

        size *= width1[2 * i] + width1[2 * i + 1] + 1;
    }

    std::unique_ptr<ValueType[]> matrix( new ValueType[ 2 * size ] );

    ValueType* matrix1 = &matrix[0];
    ValueType* matrix2 = &matrix[size];

    getMatrix( matrix1 );
    other.getMatrix( matrix2 );

    for ( IndexType k = 0; k < size; ++k )
    {
        if ( matrix1[k] != matrix2[k] ) 
        {
            // std::cout << "mismtach at k = " << k << ": " <<  matrix1[k] << ", " << matrix2[k] << std::endl;
            return false;
        }
    }

    return true;
}

template<typename ValueType>
bool Stencil<ValueType>::operator!=( const Stencil<ValueType>& other ) const
{
    return ! operator==( other );
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType Stencil<ValueType>::getValidPoints( bool valid[], const IndexType gridSizes[], const IndexType gridPos[] ) const
{
    IndexType cnt = 0;   // count the number of valid stencil neighbors

     for ( IndexType k = 0; k < nPoints(); ++k )
     {
         valid[k] = true;

         for ( IndexType i = 0; i < mNDims; ++i )
         {
             int pos = mPositions[ k * mNDims + i ];

             if ( pos < 0 )
             {
                 if ( gridPos[i] < static_cast<IndexType>( -pos ) )
                 {
                     valid[k] = false;
                 }
             }
             else 
             {
                 if ( gridPos[i] + static_cast<IndexType>( pos ) >= gridSizes[i] )
                 {
                     valid[k] = false;
                 }
             }
         }

         if ( valid[k] )
         {
            cnt++;
         }
     }

    return cnt;
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil<ValueType>::scale( const ValueType scaling )
{
    for ( IndexType p = 0; p < nPoints(); ++p )
    {
        mValues[p] *= scaling;
    }
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil<ValueType>::transpose( const Stencil<ValueType>& other )
{
    mNDims = other.mNDims;
    mValues = other.mValues;
    mPositions = other.mPositions;

    for ( size_t k = 0; k < mPositions.size(); k++ )
    {
        mPositions[ k ] = -mPositions[k];
    }
}

/* ==================================================================================== */
/*   Stencil1D  one-dimensional stencil                                                 */
/* ==================================================================================== */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Stencil1D : public Stencil<ValueType>
{
public:

    /** Create an empty stencil */

    Stencil1D();

    /** Create a default stencil with a certain number of points */

    Stencil1D( const IndexType nPoints );

    /** Construct the stencil by a matrix that covers all neighbors involved.
     *
     *  @param[in] n is the size of the matrix
     *  @param[in] matrix contains the values  
     *
     *  If n is odd, the range -n2 .. n2 is covered with n2 = ( n - 1 ) / 2
     *  If n is even, the range -n2+1 .. n2 is covered with n2 = n / 2
     */

    template<typename OtherValueType>
    Stencil1D( const IndexType n, const OtherValueType matrix[] );

    /** More convenient interface for addPoint on 1-dimensional stencil */

    void addPoint( int pos, ValueType val );

    /** More convenient interface for getPoint on 1-dimensional stencil */

    void getPoint( int& pos, ValueType& val, IndexType k ) const;

protected:

    using Stencil<ValueType>::mPositions;
    using Stencil<ValueType>::mValues;

};

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline constructors                                             */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil1D<ValueType>::Stencil1D() : Stencil<ValueType>( 1 )
{
}

template<typename ValueType>
Stencil1D<ValueType>::Stencil1D( const IndexType nPoints ) : Stencil<ValueType>( 1 )
{
    Stencil<ValueType>::reserve( nPoints );

    ValueType minusOne( -1 );
 
    switch( nPoints ) 
    {
        case 7 : 
            addPoint( -3, minusOne );
            addPoint( 3, minusOne );
            // fall through
        case 5 : 
            addPoint( -2, minusOne );
            addPoint( 2, minusOne );
            // fall through
        case 3 : 
            addPoint( -1, minusOne );
            addPoint( 1, minusOne );
            // fall through
        case 1 :
            addPoint( ( int ) 0, ValueType( nPoints - 1 ) );
            break;
        default:
            COMMON_THROWEXCEPTION( "Unsupported type of Stencil1D, #points = " << nPoints )
    }
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
Stencil1D<ValueType>::Stencil1D( const IndexType n, const OtherValueType matrix[] ) :

    Stencil<ValueType>( 1 )

{
    // if n is even, we use one element more on the right, e.g. n = 8 covers -3 : 4, 9 covers -4 : 4

    int ub = static_cast<int>( n / 2 );
    int lb = ub + 1 - static_cast<int>( n );
 
    IndexType cnt = 0;

    for ( int i = lb; i <= ub; ++i )
    {
        ValueType val = static_cast<ValueType>( matrix[ cnt++ ] );

        if ( val != 0 )
        {
            addPoint( i, val );
        }
    }

    SCAI_ASSERT_EQ_ERROR( cnt, n, "serious mismatch" )
}

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline methods                                                  */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil1D<ValueType>::addPoint( int pos, ValueType val )
{
    mPositions.push_back( pos );
    mValues.push_back( val );
}

template<typename ValueType>
void Stencil1D<ValueType>::getPoint( int& pos, ValueType& val, IndexType k ) const
{
    pos = mPositions[ k ];
    val = mValues[ k ];
}

/* ==================================================================================== */
/*   Stencil2D  two-dimensional stencil                                                 */
/* ==================================================================================== */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Stencil2D : public Stencil<ValueType>
{
public:

    /** Create an empty stencil */

    Stencil2D();

    /** Create a default stencil with a certain number of points 
     * 
     *  @param[in] nPoints is the number of stencil points, can be 5 or 9
     */
    Stencil2D( const IndexType nPoints );

    /** Create a stencil by two one-dimensional stencils */

    Stencil2D( const Stencil1D<ValueType>& stencilX, const Stencil1D<ValueType>& stencilY );

    /** More convenient interface for addPoint on Stencil2D */

    void addPoint( int posX, int posY, ValueType val );

    /** More convenient interface for getPoint on 2-dimensional stencil */

    void getPoint( int& posX, int& posY, ValueType& val, IndexType k ) const;

    /** Construct the stencil by a n1 x n2 matrix that covers all neighbors involved.
     *
     *  @param[in] n1 is the size of the matrix in the first dimension, must be odd
     *  @param[in] n2 is the size of the matrix in the second dimension, must be odd
     *  @param[in] matrix contains the values  
     */

    template<typename OtherValueType>
    Stencil2D( const IndexType n1, const IndexType n2, const OtherValueType matrix[] );

protected:

    using Stencil<ValueType>::mPositions;
    using Stencil<ValueType>::mValues;
};

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline constructors                                             */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil2D<ValueType>::Stencil2D() : Stencil<ValueType>( 2 )
{
}

template<typename ValueType>
Stencil2D<ValueType>::Stencil2D( const IndexType nPoints ) : Stencil<ValueType>( 2 )
{
    Stencil<ValueType>::reserve( nPoints );

    ValueType minusOne( -1 );
 
    switch( nPoints ) 
    {
        case 9 : 
            addPoint( -1, -1, minusOne );
            addPoint( -1,  1, minusOne );
            addPoint(  1, -1, minusOne );
            addPoint(  1,  1, minusOne );
            // fall through
        case 5 : 
            addPoint(  0, -1, minusOne );
            addPoint(  0,  1, minusOne );
            addPoint( -1,  0, minusOne );
            addPoint(  1,  0, minusOne );
            // fall through
        case 1 :
            addPoint( 0, 0, ValueType( nPoints - 1 ) );
            break;
        default:

            COMMON_THROWEXCEPTION( "Unsupported type of Stencil2D, #points = " << nPoints )
    }
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
Stencil2D<ValueType>::Stencil2D( const IndexType n1, const IndexType n2, const OtherValueType matrix[] ) :

    Stencil<ValueType>( 2 )

{
    int ub1 = static_cast<int>( n1 / 2 );
    int lb1 = ub1 + 1 - static_cast<int>( n1 );
    int ub2 = static_cast<int>( n2 / 2 );
    int lb2 = ub2 + 1 - static_cast<int>( n2 );

    IndexType cnt = 0;

    for ( int i1 = lb1; i1 <= ub1; ++i1 )
    {
        for ( int i2 = lb2; i2 <= ub2; ++i2 )
        {
            ValueType val = static_cast<ValueType>( matrix[ cnt++ ] );

            if ( val != 0 )
            {
                addPoint( i1, i2, val );
            }
        }
    }

    SCAI_ASSERT_EQ_ERROR( n1 * n2, cnt, "serious mismatch" )
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil2D<ValueType>::Stencil2D( 
    const Stencil1D<ValueType>& stencilX, 
    const Stencil1D<ValueType>& stencilY ) :

    Stencil<ValueType>( 2 )

{
    IndexType nPoints = stencilX.nPoints() + stencilY.nPoints() -  2 + 1 ;

    ValueType diagValue = 0;

    Stencil<ValueType>::reserve( nPoints );

    int pos;        // is relative, so must be signed
    ValueType val;

    for ( IndexType i = 0; i < stencilX.nPoints(); ++i )
    {
        stencilX.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( pos, 0, val );
        }
    }

    for ( IndexType i = 0; i < stencilY.nPoints(); ++i )
    {
        stencilY.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( 0, pos, val );
        }
    }

    addPoint( 0, 0, diagValue );
}

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline methods                                                  */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
void Stencil2D<ValueType>::addPoint( int posX, int posY, ValueType val )
{
    mPositions.push_back( posX );
    mPositions.push_back( posY );
    mValues.push_back( val );
}

template<typename ValueType>
void Stencil2D<ValueType>::getPoint( int& posX, int& posY, ValueType& val, IndexType k ) const
{
    posX = mPositions[ k * 2 ];
    posY = mPositions[ k * 2 + 1 ];
    val = mValues[ k ];
}

/* ==================================================================================== */
/*   Stencil3D  three-dimensional stencil                                               */
/* ==================================================================================== */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Stencil3D : public Stencil<ValueType>
{
public:

    /** Create an empty stencil */

    Stencil3D();

    /** Create a default stencil with a certain number of points */

    Stencil3D( const IndexType nPoints );

    /** Create a stencil by two one-dimensional stencils */

    Stencil3D( const Stencil1D<ValueType>& stencilX, const Stencil1D<ValueType>& stencilY, const Stencil1D<ValueType>& stencilZ );

    /** Construct the stencil by a n1 x n2 x n3 matrix that covers all neighbors involved.
     *
     *  @param[in] n1 is the size of the matrix in the first dimension, must be odd
     *  @param[in] n2 is the size of the matrix in the second dimension, must be odd
     *  @param[in] n3 is the size of the matrix in the third dimension, must be odd
     *  @param[in] matrix contains the values  
     */

    template<typename OtherValueType>
    Stencil3D( const IndexType n1, const IndexType n2, const IndexType n3, const OtherValueType matrix[] );

    /** More convenient interface for addPoint */

    void addPoint( int posX, int posY, int posZ, ValueType val );

    /** More convenient interface for getPoint on 3-dimensional stencil */

    void getPoint( int& posX, int& posY, int& posZ, ValueType& val, IndexType k ) const;

protected:

    using Stencil<ValueType>::mPositions;
    using Stencil<ValueType>::mValues;
};

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline constructors                                             */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil3D<ValueType>::Stencil3D() : Stencil<ValueType>( 3 )
{
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil3D<ValueType>::Stencil3D( const IndexType nPoints ) : Stencil<ValueType>( 3 )
{
    Stencil<ValueType>::reserve( nPoints );

    ValueType minusOne( -1 );
 
    switch( nPoints ) 
    {
        case 27 : 

            addPoint( -1, -1, -1, minusOne );
            addPoint( -1, -1,  1,  minusOne );
            addPoint( -1,  1, -1, minusOne );
            addPoint( -1,  1,  1,  minusOne );
            addPoint(  1, -1, -1, minusOne );
            addPoint(  1, -1,  1,  minusOne );
            addPoint(  1,  1, -1, minusOne );
            addPoint(  1,  1,  1,  minusOne );

            // fall through

        case 19 : 

            addPoint( -1,  0, -1, minusOne );
            addPoint( -1,  0,  1, minusOne );
            addPoint( -1, -1,  0, minusOne );
            addPoint( -1,  1,  0, minusOne );
            addPoint(  0, -1, -1, minusOne );
            addPoint(  0,  1, -1, minusOne );
            addPoint(  0, -1,  1, minusOne );
            addPoint(  0,  1,  1, minusOne );
            addPoint(  1,  0, -1, minusOne );
            addPoint(  1,  0,  1, minusOne );
            addPoint(  1, -1,  0, minusOne );
            addPoint(  1,  1,  0, minusOne );

            // fall through

        case 7 :

            addPoint( -1,  0,  0, minusOne );
            addPoint(  1,  0,  0, minusOne );
            addPoint(  0, -1,  0, minusOne );
            addPoint(  0,  1,  0, minusOne );
            addPoint(  0,  0, -1, minusOne );
            addPoint(  0,  0,  1, minusOne );

            // fall through

        case 1 :

            addPoint( 0, 0, 0, ValueType( nPoints - 1 ) );
            break;

        default:

            COMMON_THROWEXCEPTION( "Unsupported type of Stencil3D, #points = " << nPoints )
    }

    SCAI_ASSERT_EQ_DEBUG( Stencil<ValueType>::nPoints(), nPoints, "serious mismatch" )
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil3D<ValueType>::Stencil3D( 
    const Stencil1D<ValueType>& stencilX, 
    const Stencil1D<ValueType>& stencilY, 
    const Stencil1D<ValueType>& stencilZ ) : 

    Stencil<ValueType>( 3 )

{
    IndexType nPoints = stencilX.nPoints() + stencilY.nPoints() + stencilZ.nPoints() -  3 + 1 ;

    ValueType diagValue = 0;

    Stencil<ValueType>::reserve( nPoints );

    int pos;
    ValueType val;

    for ( IndexType i = 0; i < stencilX.nPoints(); ++i )
    {
        stencilX.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( pos, 0, 0, val );
        }
    }

    for ( IndexType i = 0; i < stencilY.nPoints(); ++i )
    {
        stencilY.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( 0, pos, 0, val );
        }
    }

    for ( IndexType i = 0; i < stencilZ.nPoints(); ++i )
    {
        stencilZ.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( 0, 0, pos, val );
        }
    }

    addPoint( 0, 0, 0, diagValue );
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
Stencil3D<ValueType>::Stencil3D( const IndexType n1, const IndexType n2, const IndexType n3, const OtherValueType matrix[] ) :

    Stencil<ValueType>( 3 )

{
    int ub1 = static_cast<int>( n1 / 2 );
    int lb1 = ub1 + 1 - static_cast<int>( n1 );
    int ub2 = static_cast<int>( n2 / 2 );
    int lb2 = ub2 + 1 - static_cast<int>( n2 );
    int ub3 = static_cast<int>( n3 / 2 );
    int lb3 = ub3 + 1 - static_cast<int>( n3 );

    IndexType cnt = 0;

    for ( int i1 = lb1; i1 <= ub1; ++i1 )
    {
        for ( int i2 = lb2; i2 <= ub2; ++i2 )
        {
            for ( int i3 = lb3; i3 <= ub3; ++i3 )
            {
                ValueType val = static_cast<ValueType>( matrix[ cnt++ ] );

                if ( val != 0 )
                {
                    addPoint( i1, i2, i3, val );
                }
            }
        }
    }

    SCAI_ASSERT_EQ_ERROR ( n1 * n2 * n3, cnt, "serious mismatch" )
}

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline methods for Stencil3D                                    */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void Stencil3D<ValueType>::addPoint( int posX, int posY, int posZ, ValueType val )
{
    mPositions.push_back( posX );
    mPositions.push_back( posY );
    mPositions.push_back( posZ );
    mValues.push_back( val );
}

template<typename ValueType>
void Stencil3D<ValueType>::getPoint( int& posX, int& posY, int& posZ, ValueType& val, IndexType k ) const
{
    posX = mPositions[ k * 3 ];
    posY = mPositions[ k * 3 + 1 ];
    posZ = mPositions[ k * 3 + 2 ];
    val = mValues[ k ];
}

/* ==================================================================================== */
/*   Stencil4D  four-dimensional stencil                                                */
/* ==================================================================================== */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Stencil4D : public Stencil<ValueType>
{
public:

    /** Create an empty stencil */

    Stencil4D();

    /** Create a default stencil with a certain number of points */

    Stencil4D( const IndexType nPoints );

    /** Create a stencil by two one-dimensional stencils */

    Stencil4D( 
        const Stencil1D<ValueType>& stencilX, 
        const Stencil1D<ValueType>& stencilY, 
        const Stencil1D<ValueType>& stencilZ,
        const Stencil1D<ValueType>& stencilT );

    /** More convenient interface for addPoint */

    void addPoint( int posX, int posY, int posZ, int posT, ValueType val );

    /** More convenient interface for getPoint on 4-dimensional stencil */

    void getPoint( int& posX, int& posY, int& posZ, int& posT, ValueType& val, IndexType k ) const;

protected:

    using Stencil<ValueType>::mPositions;
    using Stencil<ValueType>::mValues;
};

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline constructors                                             */
/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil4D<ValueType>::Stencil4D() : Stencil<ValueType>( 4 )
{
}

template<typename ValueType>
Stencil4D<ValueType>::Stencil4D( const IndexType nPoints ) : Stencil<ValueType>( 4 )
{
    Stencil<ValueType>::reserve( nPoints );

    ValueType minusOne( -1 );
 
    switch( nPoints ) 
    {
        case 9 :

            addPoint( -1,  0,  0,  0, minusOne );
            addPoint(  1,  0,  0,  0, minusOne );
            addPoint(  0, -1,  0,  0, minusOne );
            addPoint(  0,  1,  0,  0, minusOne );
            addPoint(  0,  0, -1,  0, minusOne );
            addPoint(  0,  0,  1,  0, minusOne );
            addPoint(  0,  0,  0, -1, minusOne );
            addPoint(  0,  0,  0,  1, minusOne );

            // fall through

        case 1 :

            addPoint( 0, 0, 0, 0, ValueType( nPoints - 1 ) );
            break;

        default:

            COMMON_THROWEXCEPTION( "Unsupported type of Stencil4D, #points = " << nPoints )
    }

    SCAI_ASSERT_EQ_DEBUG( Stencil<ValueType>::nPoints(), nPoints, "serious mismatch" )
}

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
Stencil4D<ValueType>::Stencil4D( 
    const Stencil1D<ValueType>& stencilX, 
    const Stencil1D<ValueType>& stencilY, 
    const Stencil1D<ValueType>& stencilZ,
    const Stencil1D<ValueType>& stencilT ) :

    Stencil<ValueType>( 4 )

{
    IndexType nPoints = stencilX.nPoints() + stencilY.nPoints() + stencilZ.nPoints() + stencilT.nPoints() - 4 + 1 ;

    ValueType diagValue = 0;

    Stencil<ValueType>::reserve( nPoints );

    int pos;         // must be signed 
    ValueType val;

    for ( IndexType i = 0; i < stencilX.nPoints(); ++i )
    {
        stencilX.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( pos, 0, 0, 0, val );
        }
    }

    for ( IndexType i = 0; i < stencilY.nPoints(); ++i )
    {
        stencilY.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( 0, pos, 0, 0, val );
        }
    }

    for ( IndexType i = 0; i < stencilZ.nPoints(); ++i )
    {
        stencilZ.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( 0, 0, pos, 0, val );
        }
    }

    for ( IndexType i = 0; i < stencilT.nPoints(); ++i )
    {
        stencilT.getPoint( pos, val, i );

        if ( pos == 0 )
        {
            diagValue += val;
        }
        else
        {
            addPoint( 0, 0, 0, pos, val );
        }
    }

    addPoint( 0, 0, 0, 0, diagValue );
}

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline methods for Stencil4D                                    */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void Stencil4D<ValueType>::addPoint( int posX, int posY, int posZ, int posT, ValueType val )
{
    mPositions.push_back( posX );
    mPositions.push_back( posY );
    mPositions.push_back( posZ );
    mPositions.push_back( posT );
    mValues.push_back( val );
}

template<typename ValueType>
void Stencil4D<ValueType>::getPoint( int& posX, int& posY, int& posZ, int& posT, ValueType& val, IndexType k ) const
{
    posX = mPositions[ k * 4 ];
    posY = mPositions[ k * 4 + 1 ];
    posZ = mPositions[ k * 4 + 2 ];
    posT = mPositions[ k * 4 + 3 ];
    val = mValues[ k ];
}

/* ---------------------------------------------------------------------------------- */

}  // namespace common

}  // namespace scai

