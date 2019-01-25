/**
 * @file mandelbrot.cpp
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
 * @brief mandelbrot.cpp is an example for Julia Set with GPU.
 * @author Vanessa Wolff
 * @date 18.08.2016
 */

#include <stdio.h>

#ifdef __APPLE__
#include <glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

using namespace scai;
using namespace lama;
using hmemo::ReadAccess;
using hmemo::WriteAccess;

typedef DefaultReal ValueType;

/*
*defining a RGB struct to color the pixel
*/

#define DIMx 1440
#define DIMy 841
#define maxIteration 200

struct Type_rgb
{
    ValueType r;
    ValueType g;
    ValueType b;
};
/*
* pixel variable contain the value of the color pixel in
* the picture.
*/
struct Type_rgb pixels[DIMy* DIMx];

/*
* function mandelbrotset find where the number is in
* mandelbrotset or not and also assign a color to that
* coordinate with that iteration pattern.
*/

void Julia()
{
    auto ind   = denseVectorZero<ValueType>( DIMy * DIMx );
    auto cReal = denseVectorZero<ValueType>( DIMx * DIMy );
    auto cImag = denseVectorZero<ValueType>( DIMx * DIMy );
    auto mag   = denseVectorZero<ValueType>( DIMx * DIMy );
    auto help  = denseVectorZero<ValueType>( DIMx * DIMy );
    auto help2 = denseVectorZero<ValueType>( DIMx * DIMy );

    auto aReal = denseVector<ValueType>( DIMx * DIMy, -0.8 );
    auto aImag = denseVector<ValueType>( DIMx * DIMy, 0.156 );

    auto yreal = denseVectorLinear<ValueType>( DIMy, 0, 1 );
    auto xreal = denseVectorLinear<ValueType>( DIMx, 0, 1 );

    auto eye1 = denseVector( DIMx, ValueType( 1 ) );
    auto eye2 = denseVector( DIMy, ValueType( 1 ) );

    ValueType w = 0.5 * DIMx;
    ValueType h = 0.5 * DIMy;

    xreal = eye1 * w - xreal;
    xreal = xreal / w;
    xreal *= 1.5;
    yreal = eye2 * h - yreal;
    yreal = yreal / h;
    yreal *= 1.5;

    {
        ReadAccess<ValueType> x1( xreal.getLocalValues() );
        ReadAccess<ValueType> y1( yreal.getLocalValues() );
        WriteAccess<ValueType> c1( cReal.getLocalValues() );
        WriteAccess<ValueType> c2( cImag.getLocalValues() );

        for ( IndexType i = 0; i < DIMy; i++ )
        {
            for ( IndexType j = 0; j < DIMx; j++ )
            {
                c1[i * DIMx + j] = x1[j];
                c2[i * DIMx + j] = y1[i];
            }
        }
    }
    int iteration = 0;

    while ( iteration < maxIteration )
    {
        help = cReal * cReal;
        help2 = cImag * cImag;
        cImag = ValueType( 2 ) * cReal * cImag;
        cReal = help - help2;
        cReal += aReal;
        cImag += aImag;
        help = 0.0;
        mag = cReal * cReal;
        help = cImag * cImag;
        mag += help;
        {
            ReadAccess<ValueType> m( mag.getLocalValues() );
            WriteAccess<ValueType> in( ind.getLocalValues() );
            #pragma omp parallel for

            for ( IndexType i = 0; i < DIMx * DIMy; i++ )
            {
                if ( m[i] < 4 )
                {
                    in[i] ++;
                }
            }
        }
        iteration++;
    }

    ReadAccess<ValueType> in( ind.getLocalValues() );
    #pragma omp parallel for

    for ( IndexType i = 0; i < DIMx * DIMy; i++ )
    {
        if ( in[i] >= maxIteration )
        {
            // setting color pixel to dark blue
            pixels[i].r = 0;
            pixels[i].g = 207;
            pixels[i].b = 207;
        }
        else
        {
            // setting color pixel to light blue
            pixels[i].r = 0;
            pixels[i].g = 0;
            pixels[i].b = 73;
        }
    }
}

void Init( )
{
    /*
    * Basic Opengl initialization.
    */
    glViewport( 0, 0, DIMx, DIMy );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity( );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity( );
    gluOrtho2D( 0, DIMx, 0, DIMy );
    /*
    * Initializing all the pixels to white.
    */
    #pragma omp parallel for

    for ( int i = 0; i < DIMy * DIMx; i++ )
    {
        pixels[i].r = 1;
        pixels[i].g = 1;
        pixels[i].b = 1;
    }

    // call of Julia Set Routine
    Julia();
}

void onDisplay()
{
    /*
    * Clearing the initial buffer
    */
    glClearColor( 1, 1, 1, 0 );
    glClear( GL_COLOR_BUFFER_BIT );
    /*
    * Draw the complete Mandelbrot set picture.
    */
    glDrawPixels( DIMx, DIMy, GL_RGB, GL_FLOAT, pixels );
    glutSwapBuffers();
}

int main( int argc, char** argv )
{
    /*
    * Here basic Opengl initialization.
    */
    glutInit( &argc, argv );
    glutInitWindowSize ( DIMx, DIMy );
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
    glutInitWindowPosition ( 100, 100 );
    glutCreateWindow ( "Mandelbrotset" );
    Init ();
    /*
    * connecting the Display function
    */
    glutDisplayFunc( onDisplay );
    /*
    * starting the activities
    */
    glutMainLoop();
    return 0;
}
