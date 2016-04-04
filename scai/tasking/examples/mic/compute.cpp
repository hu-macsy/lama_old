
#include <cstdlib>

__declspec( target(mic, host) )
void compute( size_t NSIZE )
{
    double* a = new double[ NSIZE * NSIZE ];
    double* b = new double[ NSIZE * NSIZE ];
    double* c = new double[ NSIZE * NSIZE ];

    for ( int i = 0; i < NSIZE; ++i )
    {
        for ( int j = 0; j < NSIZE; ++j )
        {
            a[ i * NSIZE + j ] = 1;
            b[ i * NSIZE + j ] = 1;
            c[ i * NSIZE + j ] = 0;
        }
    }

    #pragma omp parallel for

    for ( int i = 0; i < NSIZE; ++i )
    {
        for ( int j = 0; j < NSIZE; ++j )
        {
            double tmp = 0;
            for ( int k = 0; k < NSIZE; ++k )
            {
                tmp += a[ i * NSIZE + k ] * b[ k * NSIZE + j ];
            }
            c[ i * NSIZE + j ] = tmp;
        }
    }

    delete[] a;
    delete[] b;
    delete[] c;
}

