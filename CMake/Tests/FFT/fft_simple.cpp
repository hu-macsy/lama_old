#include<vector>
#include<fftw3.h>

int main( void )
{
    int N = 4;
    fftw_complex* in;
    fftw_complex* out;
    fftw_plan my_plan;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N );
    out= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N );
    my_plan = fftw_plan_dft_1d( N, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute( my_plan );
    fftw_destroy_plan( my_plan );
    fftw_free( in );
    fftw_free( out );
}
