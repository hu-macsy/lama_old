/**
 * @file GPI/hello_world.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief ToDo: Missing description in ./GPI/hello_world.cpp
 * @author Thomas Brandes
 * @date 30.03.2016
 */
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <GASPI.h>

#define GASPI_CALL( call )                                                \
    {                                                                     \
        gaspi_return_t status = call;                                     \
                                                                          \
        gaspi_printf( "status = %d\n", status) ;               \
        if ( status != GASPI_SUCCESS )                                    \
        {                                                                 \
            std::ostringstream errorStr;                                  \
            errorStr << "GPI error in line " << __LINE__ ;                \
            errorStr << " of file " << __FILE__ << ": ";                  \
            errorStr << #call << "\n";                                    \
            gaspi_printf( "%s\n", errorStr.str().c_str() );               \
            exit(-1);                                                     \
        }                                                                 \
    }

int main(int argc, char *argv[])
{

  gaspi_printf( "Main\n" );

  GASPI_CALL( gaspi_proc_init(GASPI_BLOCK) );

  gaspi_printf( "Started\n" );

  gaspi_rank_t rank, num;

  gaspi_proc_rank(&rank);
  gaspi_proc_num(&num);

  gaspi_printf( "Hello from rank %d of %d\n", rank, num);

  gaspi_proc_term(GASPI_BLOCK);
}

