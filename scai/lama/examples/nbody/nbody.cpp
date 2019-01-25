/**
 * @file nbody.cpp
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
 * @brief nbody.cpp is an example for an nbody simulation - here one fixed big body in the middle of the system,
 *        all others movable (random positions, masses and velocities in the beginning)
 * @author Vanessa Wolff
 * @date 26.09.2016
 */

#include <stdio.h>

#ifdef __APPLE__
#include <glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include <scai/lama.hpp>

#include <scai/lama/matrix/_Matrix.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <scai/dmemo/Distributed.hpp>

#include <scai/common/Math.hpp>

#include <math.h>

using namespace scai;
using namespace lama;
using hmemo::ReadAccess;
using hmemo::WriteAccess;
using dmemo::Distributed;
using common::Math;
using common::UnaryOp;


typedef DefaultReal ValueType;


#define DIMx 900                // dimension of displayed window
#define DIMy 900
#define maxIteration 500        // max. number of iteration
#define nBodies 200             // number of bodies in the system

ValueType maxMass  = 1.98892e30;  // max. mass of biggest body
ValueType dt = 1e11;              // time step
ValueType softening = 3e4f;       // softening parameter (used for force computation to avoid infinities)
ValueType G = 6.673e-11;            // gravity constant
ValueType radius = 1e18;           // radius of the universe

IndexType iter = 0;               // global iteration count

const ValueType ZERO = 0;               // typed 0, used to fill vectors

auto x = denseVector( nBodies, ZERO );
auto y = denseVector( nBodies, ZERO );

auto vx = denseVector( nBodies, ZERO );
auto vy = denseVector( nBodies, ZERO );

auto fx = denseVector( nBodies, ZERO );
auto fy = denseVector( nBodies, ZERO );

auto mass        = denseVector( nBodies, ZERO );
auto inversemass = denseVector( nBodies, ZERO );

void randomBodies( )
{
    // Random Positions (relative to the radius of the universe at the beginning)

    x.setRandom( nBodies, 2 );
    x -= 1;
    x *= radius * ( ValueType )Math::exp( -1.8 ) ;

    y.setRandom( nBodies, 2 );
    y -= 1;
    y *= radius * ( ValueType )Math::exp( -1.8 ) ;


    // circular velocities for each particle depenent on particle position
    for ( IndexType i = 0; i < nBodies; i++ )
    {
        ValueType px = x.getValue( i );
        ValueType py = y.getValue( i );

        ValueType denum = Math::sqrt( px * px + py * py );
        ValueType num = ( 6.67e-11 ) * 1e6 * maxMass;
        ValueType magv = Math::sqrt( num / denum );

        ValueType thetav = ( ValueType )M_PI / 2.0 - ( ValueType )atan( Math::abs( py / px ) );

        vx.setValue( i, -1.0 * copysign( 1.0, py ) * ( ValueType )cos( thetav ) * magv );
        vy.setValue( i, copysign( 1.0, px ) * ( ValueType )sin( thetav ) * magv );
    }


    // Random mass of Particles
    mass.setRandom( mass.getDistributionPtr(), 1 );  // random numbers between 0 and 1
    mass *= maxMass ;

    // put a heavy body in the center
    x.setValue( 0, 0.0 );
    y.setValue( 0, 0.0 );
    vx.setValue( 0, 0.0 );
    vy.setValue( 0, 0.0 );
    mass.setValue( 0, maxMass );

    inversemass = 1 / mass;
}

void updateParticles( )
{
    DenseVector<ValueType> help;

    // update velocities  v  = v + (dt * fx)/mass

    help = dt * fx;
    help *= inversemass;  // elementwise multiplication
    vx += help;

    help = 0.0;

    help = dt * fy;
    help *= inversemass;  
    vy += help;

    // update coordinates new-pos = old-pos + dt * velocity
    x += dt * vx;
    y += dt * vy;
}

void resetForce( )
{
    fx = 0.0;
    fy = 0.0;
}

void computeForce( IndexType i, IndexType j )
{
    ValueType dx = x.getValue( j ) - x.getValue( i );
    ValueType dy = y.getValue( j ) - y.getValue( i );

    // distance between the bodies
    ValueType dist = Math::sqrt( dx * dx + dy * dy );

    // force between the bodies = ( gravity * mass(a) * mass(b) ) / ( distance^2 + softening^2)         -> softening used to avoid infinities
    ValueType force = ( G * mass.getValue( i ) * mass.getValue( j ) ) / ( dist * dist + softening * softening );

    // force update for bodie a : force(a) = sum of the forces between a and all other bodies
    // force(a) += dx * force between a and b / distance between a and b
    ValueType f1 = fx.getValue( i ) + dx * force / dist;
    ValueType f2 = fy.getValue( i ) + dy * force / dist;

    // update global Force
    fx.setValue( i, f1 );
    fy.setValue( i, f2 );

}

void Nbody()
{
    for ( IndexType i = 0; i < nBodies; i++ )
    {
        for ( IndexType j = 0; j < nBodies; j++ )
        {
            if ( i != j )
                // compute new force between particle i and j
            {
                computeForce( i, j );
            }
        }
    }

    // update global particles coordinates and velocities with respect to new interacting forces
    updateParticles( );

}


//draw particles
void onDisplay()
{
    glClear( GL_COLOR_BUFFER_BIT );

    for ( IndexType i = 0; i < nBodies; i++ )
    {
        // get position and mass of each body
        ValueType xpos = x.getValue( i ) * ( DIMx - 5.0 ) / radius ;
        ValueType ypos = y.getValue( i ) * ( DIMy - 5.0 ) / radius;
        ValueType m = mass.getValue( i );
        m *= 1.0 / ( ( ValueType )maxMass );

        glColor3f( 0.90, 0.91, 0.98 );

        if ( i == 0 )   // large body in the middle is colored red
        {
            glColor3f( 1.0 , 0.0, 0.0 );
        }

        glBegin( GL_POLYGON );

        for ( float a = 0; a < 2 * M_PI; a += 0.1 )
        {
            glVertex2f( ( GLfloat )( 5.0 * m * cos( a ) + xpos ), ( GLfloat )( 5.0 * m * sin( a ) + ypos ) );
        }

        glEnd();
    }

    glFlush();
    glutSwapBuffers();
}



// timer function to redisplay changes
void timer( int /*t*/ )
{
    onDisplay();

    // resett all forces between particles to 0 to compute them again
    resetForce( );

    // calculate new placement of the bodies -> movement because of own velocity and force directet
    Nbody();

    // cout iteration recompute until max iteraion
    iter++;

    if ( iter <= maxIteration )
    {
        glutTimerFunc( 1, timer, 0 );
    }
}



int main( int argc, char** argv )
{
    // Here basic Opengl initialization.
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB );
    glutInitWindowSize( DIMx, DIMy );
    glutInitWindowPosition( 50, 50 );
    glutCreateWindow( "N-Body Simmulation" );

    glClearColor( 0, 0, 0, 1 );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho( -DIMx, DIMx, DIMy, -DIMy, 0, 1 );

    // create Random Bodies' placement and velocities
    randomBodies( );

    glutDisplayFunc( onDisplay );

    timer( 0 );

    glutMainLoop();

    return 0;
}


