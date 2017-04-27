/**
 * @file lama/Stencil.hpp
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
 * @brief Definition of stencil classes
 * @author Thomas Brandes
 * @date 13.04.2017
 */

#pragma once

// for dll_import

#include <scai/common/config.hpp>
#include <scai/common/Utils.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>

// std
#include <stdint.h>

#define SCAI_STENCIL_MAX_POINTS 128

namespace scai
{

namespace lama
{

/** Common base class for all stencils. 
 *
 *  A stencil describes in an n-dimensional grid how one element is connected
 *  with a certain number of neighbors.
 *
 *  Please note that a stencil does not contain any information about the grid size
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

    /** Add a stencil point */

    void addPoint( const int relpos[], const ValueType val );

    /** Return the total number of stencil points. */

    IndexType size() const;

    /** Get stencil positions in linearized grid 
     *  
     *  @param[in] gridDistances contains for each dim the distance between two neighbored elements in that dim
     *  @param[out] pos with entry fore each stencil point the relative position in linearized grid
     */
    inline void getLinearPositions( int pos[], IndexType gridDistances[] ) const;

    /** Get valid positions */

    inline IndexType getValidPoints( bool valid[], const IndexType gridSizes[], const IndexType gridPos[] ) const;

    /** Return widths of a stencil in each dimension, both directions
     *
     *  @param[out] lb is array with width to the left for each dimension
     *  @param[out] ub is array with width to the right for each dimension
     */
    inline void getWidths( IndexType lb[], IndexType ub[] ) const;

    /** Compares two stencils for equality. */

    bool operator==( const Stencil& other ) const;

    /** Compares two stencils for inequality. */

    bool operator!=( const Stencil& other ) const;

    /** Get pointer to the array with stencil positions, size is nPoints() * nDims() */

    const int* positions() const;

    /** Get pointer to the array with stencil values, size is nPoints() */

    const ValueType* values() const;

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
void Stencil<ValueType>::getLinearPositions( int pos[], IndexType gridDistances[] ) const
{
     for ( IndexType k = 0; k < nPoints(); ++k )
     {
         pos[k] = 0;
         
         for ( IndexType i = 0; i < mNDims; ++i )
         {
             pos[k] += mPositions[ k * mNDims + i ] * static_cast<int>( gridDistances[ i ] );
         }
     }
}

template<typename ValueType>
void Stencil<ValueType>::getWidths( IndexType lb[], IndexType ub[] ) const
{
    // initialize lb and ub arrays with 0

    for ( IndexType i = 0; i < mNDims; ++i )
    {
        lb[i] = 0;
        ub[i] = 0;
    }

    for ( IndexType k = 0; k < nPoints(); ++k )
    {
         for ( IndexType i = 0; i < mNDims; ++i )
         {
             int pos = mPositions[ k * mNDims + i ];

             if ( pos < 0 )
             {
                 IndexType lbi = static_cast<IndexType>( -pos );

                 if ( lbi > lb[ i ] )
                 { 
                     lb[i] = lbi;    // new max value here
                 }
             }
          
             if ( pos > 0 )
             {
                 IndexType ubi = static_cast<IndexType>( pos );

                 if ( ubi > ub[ i ] )
                 { 
                     ub[i] = ubi;    // new max value here
                 }
             }
         }
    }
}

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
class COMMON_DLL_IMPORTEXPORT Stencil1D : public Stencil<ValueType>
{
public:

    /** Create an empty stencil */

    Stencil1D();

    /** Create a default stencil with a certain number of points */

    Stencil1D( const IndexType nPoints );

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
        {
            addPoint( -3, minusOne );
            addPoint( 3, minusOne );
        }
        case 5 : 
        {
            addPoint( -2, minusOne );
            addPoint( 2, minusOne );
        }
        case 3 : 
        {
            addPoint( -1, minusOne );
            addPoint( 1, minusOne );
        }
        case 1 :
        {
            addPoint( ( int ) 0, ValueType( nPoints - 1 ) );
            break;
        }
        default:
            COMMON_THROWEXCEPTION( "Unsupported type of Stencil1D, #points = " << nPoints )
    }
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

/* ------------------------------------------------------------------------------------ */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT Stencil2D : public Stencil<ValueType>
{
public:

    /** Create an empty stencil */

    Stencil2D();

    /** Create a default stencil with a certain number of points 
     * 
     *  @param[in] nPoints, can be 5 or 9
     */
    Stencil2D( const IndexType nPoints );

    /** Create a stencil by two one-dimensional stencils */

    Stencil2D( const Stencil1D<ValueType>& stencilX, const Stencil1D<ValueType>& stencilY );

    /** More convenient interface for addPoint on Stencil2D */

    void addPoint( int posX, int posY, ValueType val );

    /** More convenient interface for getPoint on 2-dimensional stencil */

    void getPoint( int& posX, int& posY, ValueType& val, IndexType k ) const;

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
        {
            addPoint( -1, -1, minusOne );
            addPoint( -1,  1, minusOne );
            addPoint(  1, -1, minusOne );
            addPoint(  1,  1, minusOne );
        }
        case 5 : 
        {
            addPoint(  0, -1, minusOne );
            addPoint(  0,  1, minusOne );
            addPoint( -1,  0, minusOne );
            addPoint(  1,  0, minusOne );
        }
        case 1 :
        {
            addPoint( 0, 0, ValueType( nPoints - 1 ) );
            break;
        }
        default:

            COMMON_THROWEXCEPTION( "Unsupported type of Stencil2D, #points = " << nPoints )
    }
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

    IndexType pos;
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

/* ------------------------------------------------------------------------------------ */

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
        {
            addPoint( -1, -1, -1, minusOne );
            addPoint( -1, -1,  1,  minusOne );
            addPoint( -1,  1, -1, minusOne );
            addPoint( -1,  1,  1,  minusOne );
            addPoint(  1, -1, -1, minusOne );
            addPoint(  1, -1,  1,  minusOne );
            addPoint(  1,  1, -1, minusOne );
            addPoint(  1,  1,  1,  minusOne );

            // no break here, continue with points for Stencil3D( 19 )
        }
        case 19 : 
        {
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

            // no break here, continue with points for Stencil3D( 7 )
        }
        case 7 :
        {
            addPoint( -1,  0,  0, minusOne );
            addPoint(  1,  0,  0, minusOne );
            addPoint(  0, -1,  0, minusOne );
            addPoint(  0,  1,  0, minusOne );
            addPoint(  0,  0, -1, minusOne );
            addPoint(  0,  0,  1, minusOne );

            // no break here, continue with points for Stencil3D( 1 )
        }
        case 1 :
        {
            addPoint( 0, 0, 0, ValueType( nPoints - 1 ) );
            break;
        }
        default:

            COMMON_THROWEXCEPTION( "Unsupported type of Stencil3D, #points = " << nPoints )

        SCAI_ASSERT_EQ_DEBUG( Stencil<ValueType>::nPoints(), nPoints, "serious mismatch" )
    }
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

    IndexType pos;
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

/* ---------------------------------------------------------------------------------- */

}  // namespace lama

}  // namespace scai

