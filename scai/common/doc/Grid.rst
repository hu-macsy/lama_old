.. _Grid:

Grid
====

Grid is a simple class to define an n-dimensional grid. A grid itself
does not contain any data. Its main purpose
is to support linear addressing of the grid points that are assumed to
be ordered in row-major order.

.. code-block:: c++

    Grid grid1( 10 );
    Grid grid2( 2, 5 );
    Grid grid3( 2, 3, 2 ); 

    IndexType p[] = { 1, 4 };
    grid2.linearPos( p ) -> returns 9
    IndexType q[] = { 1, 1, 1 };
    grid3.linearPos( q ) -> returns 9 

    grid3.gridPos( 10, q );   // q = { 1, 2, 0 }
