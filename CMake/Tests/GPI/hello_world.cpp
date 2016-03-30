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

