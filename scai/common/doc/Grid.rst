.. _Grid:

Grid
====

Grid is a simple base class to define an n-dimensional grid. A grid itself
does not contain any data. Its main purpose
is to support linear addressing of the grid points that are assumed to
be ordered in row-major order.

For each dimension (up to 4) there is an own derived class.

.. code-block:: c++

    #include <scai/common/Grid.hpp>

    common::Grid1D grid1( 10 );
    common::Grid2D grid2( 2, 5 );
    common::Grid3D grid3( 2, 3, 2 ); 

    grid2.linearPos( 1, 4 ) -> returns 9
    grid3.linearPos( 1, 1, 1 ) -> returns 9 

    for ( IndexType i = 0; i < grid3.size(); ++i )
    {
        IndexType q[3];
        grid3.gridPos( q, i );   
        std::cout << "Element " << i << " in grid " << grid3 << " is " 
                  << q[0] << ", " << q[1] << ", " << q[2] << std::endl;
    }

The output of the loop over the grid elements shows that the grid is traversed in
row major order.

.. code-block:: c++

   Element 0 in grid Grid3D( 2 x 3 x 2 ) is 0, 0, 0
   Element 1 in grid Grid3D( 2 x 3 x 2 ) is 0, 0, 1
   Element 2 in grid Grid3D( 2 x 3 x 2 ) is 0, 1, 0
   Element 3 in grid Grid3D( 2 x 3 x 2 ) is 0, 1, 1
   Element 4 in grid Grid3D( 2 x 3 x 2 ) is 0, 2, 0
   Element 5 in grid Grid3D( 2 x 3 x 2 ) is 0, 2, 1
   Element 6 in grid Grid3D( 2 x 3 x 2 ) is 1, 0, 0
   Element 7 in grid Grid3D( 2 x 3 x 2 ) is 1, 0, 1
   Element 8 in grid Grid3D( 2 x 3 x 2 ) is 1, 1, 0
   Element 9 in grid Grid3D( 2 x 3 x 2 ) is 1, 1, 1
   Element 10 in grid Grid3D( 2 x 3 x 2 ) is 1, 2, 0
   Element 11 in grid Grid3D( 2 x 3 x 2 ) is 1, 2, 1

Like for all data structures in LAMA indexing starts always with index 0.
