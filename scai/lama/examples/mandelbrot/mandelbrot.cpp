#include <stdio.h>

#include <GL/glut.h>

int  main( int argc, char **argv )
{
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA );
    glutInitWindowSize( 100, 100 );
    glutCreateWindow( "Mandelbrot" );
    // glutDisplayFunc( renderFunc );
    // glutMouseFunc( mouseFunc );
    // glutReshapeFunc( reshapeFunc );
    // glutKeyboardFunc( keyboardFunc );
    // glutMotionFunc( motionFunc );

    glutMainLoop();


	return 0;
}
