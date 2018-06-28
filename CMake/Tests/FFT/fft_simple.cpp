/**
 * @file FFT/fft_simple.cpp
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
 * @brief ToDo: Missing description in ./FFT/fft_simple.cpp
 * @author Thomas Brandes
 * @date 17.04.2018
 */
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
