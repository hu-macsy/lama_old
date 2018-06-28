/**
 * @file OpenMP/example.cpp
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
 * @brief Simple example program that runs only with OpenMP enabled
 * @author Thomas Brandes
 * @date 18.04.2018
 */

#include <iostream>

#include <omp.h>

int main( int, char* [] )
{
    std::cout << "#Threads = " << omp_get_max_threads() << std::endl;
    std::cout << "Version " << _OPENMP << std::endl;

    #pragma omp parallel
    {
        #pragma omp critical
        {
            std::cout << "Thread " << omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
        }
    }
}
