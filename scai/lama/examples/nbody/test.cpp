/**
 * @file lama/examples/nbody/test.cpp
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
 * @brief ToDo: Missing description in ./lama/examples/nbody/test.cpp
 * @author Eric Schricker
 * @date 05.10.2016
 */
#ifdef __APPLE__
#include <glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

GLfloat xRotated, yRotated, zRotated;
GLdouble radius = 1;

#define DIMx 1440
#define DIMy 841


void display( void );
void reshape( int x, int y );


int main ( int argc, char** argv )
{
    glutInit( &argc, argv );
    glutInitWindowSize( DIMx, DIMy );
    glutCreateWindow( "Solid Sphere" );
    xRotated = yRotated = zRotated = 30.0;
    xRotated = 43;
    yRotated = 50;

    glutDisplayFunc( display );
    glutReshapeFunc( reshape );
    glutMainLoop();
    return 0;
}

void display( void )
{

    glMatrixMode( GL_MODELVIEW );
    // clear the drawing buffer.
    glClear( GL_COLOR_BUFFER_BIT );
    // clear the identity matrix.
    glLoadIdentity();
    // traslate the draw by z = -4.0
    // Note this when you decrease z like -8.0 the drawing will looks far , or smaller.
    glTranslatef( 0.0, 0.0, -5.0 );
    // Red color used to draw.
    glColor3f( 0.9, 0.3, 0.2 );
    // changing in transformation matrix.
    // rotation about X axis
    glRotatef( xRotated, 1.0, 0.0, 0.0 );
    // rotation about Y axis
    glRotatef( yRotated, 0.0, 1.0, 0.0 );
    // rotation about Z axis
    glRotatef( zRotated, 0.0, 0.0, 1.0 );
    // scaling transfomation
    glScalef( 1.0, 1.0, 1.0 );
    // built-in (glut library) function , draw you a sphere.
    glutSolidSphere( radius / 10 , 20, 20 );
    // Flush buffers to screen

    //glFlush();
    // sawp buffers called because we are using double buffering
    glutSwapBuffers();
}

void reshape( int x, int y )
{
    if ( y == 0 || x == 0 )
    {
        return;
    }

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 39.0, ( GLdouble )x / ( GLdouble )y, 0.6, 21.0 );
    glMatrixMode( GL_MODELVIEW );
    glViewport( 0, 0, x, y ); //Use the whole window for rendering
}
